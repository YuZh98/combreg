#include <RcppArmadillo.h>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

// Hit-and-run sampler for the dual-certificate polytope
//   D(y, zeta) = { u >= 0 : u supported on active constraints,
//                  (A^T u)_j <= zeta_j if y_j = 1, >= zeta_j if y_j = 0 }
// targeting either the exponential kernel exp(-rho * sum(u)) or the
// half-Gaussian kernel exp(-0.5 * sum(u^2)).
//
// RNG: R's generator is not thread-safe, so each chain (row) gets its own
// std::mt19937_64 stream, seeded from R's RNG before the parallel region.
// Results are therefore invariant to the number of threads and reproducible
// under set.seed().
//
// Optimizations relative to a naive implementation:
//   - nonnegativity handled per-coordinate instead of via -I augmentation
//   - single gemv per iteration for the constraint bounds
//   - A * x maintained incrementally (axpy)
//   - direction sampled over active coordinates only
//   - exponential kernel: sign-flip instead of rejection when sum(v) < 0
//   - buffers reused across iterations

namespace {

struct Rng {
  std::mt19937_64 eng;
  std::normal_distribution<double> norm;
  std::uniform_real_distribution<double> unif;

  explicit Rng(uint64_t seed) : eng(seed), norm(0.0, 1.0), unif(0.0, 1.0) {}
  double rnorm() { return norm(eng); }
  double runif() { return unif(eng); }
};

// Rmath pnorm/qnorm are pure functions (no RNG state): thread-safe.
double sample_truncnorm_icdf(Rng& rng, double mean, double sd,
                             double lower, double upper) {
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  a = std::max(-37.0, std::min(37.0, a));
  b = std::max(-37.0, std::min(37.0, b));

  double Phi_a = R::pnorm(a, 0.0, 1.0, 1, 0);
  double Phi_b = R::pnorm(b, 0.0, 1.0, 1, 0);

  if (Phi_b - Phi_a < 1e-15) {
    return mean + sd * (a + b) / 2.0;
  }

  double u = rng.runif();
  double Phi_x = Phi_a + u * (Phi_b - Phi_a);
  double x_std = R::qnorm(Phi_x, 0.0, 1.0, 1, 0);
  return mean + sd * x_std;
}

arma::vec hit_and_run_chain(
    const arma::mat& A,
    const arma::vec& z,
    const arma::vec& greater_equal,
    arma::vec x0,
    int n_iter,
    const arma::vec& kappa,
    double rho,
    const std::string& kernel,
    int max_dir_tries,
    double bound_truncation,
    Rng& rng
) {
  const int d = x0.n_elem;
  const int m = A.n_rows;
  const bool exponential = (kernel == "exponential");

  // Transform constraints to A_trans * x <= z_trans form.
  arma::mat A_trans(m, d);
  arma::vec z_trans(m);
  for (int j = 0; j < m; ++j) {
    if (greater_equal[j] > 0.5) {
      A_trans.row(j) = -A.row(j);
      z_trans[j]     = -z[j];
    } else {
      A_trans.row(j) = A.row(j);
      z_trans[j]     = z[j];
    }
  }

  arma::uvec active_idx = arma::find(kappa > 0.5);
  const int n_active    = active_idx.n_elem;

  arma::vec direction(d, arma::fill::zeros);
  arma::vec Ad(m);
  arma::vec Ax = A_trans * x0;
  arma::vec rand_active(n_active);

  for (int iter = 1; iter < n_iter; ++iter) {

    // ----- STEP 1: sample direction (active coords only) -----
    bool direction_valid = false;
    double direction_sum = 0.0;

    for (int tries = 0; tries < max_dir_tries && !direction_valid; ++tries) {

      for (int i = 0; i < n_active; ++i) {
        rand_active[i] = rng.rnorm();
      }

      double norm2 = 0.0;
      for (int i = 0; i < n_active; ++i) norm2 += rand_active[i] * rand_active[i];
      if (norm2 < 1e-24) continue;
      double inv_norm = 1.0 / std::sqrt(norm2);

      direction.zeros();
      for (int i = 0; i < n_active; ++i) {
        direction[active_idx[i]] = rand_active[i] * inv_norm;
      }

      direction_sum = 0.0;
      for (int i = 0; i < n_active; ++i) {
        direction_sum += direction[active_idx[i]];
      }

      if (exponential) {
        // Kernel requires sum(v) > 0; direction law is symmetric, so flip
        // sign instead of rejecting.
        if (direction_sum < -1e-12) {
          direction      = -direction;
          direction_sum  = -direction_sum;
        }
        if (direction_sum > 1e-12) direction_valid = true;
      } else {
        direction_valid = true;
      }
    }

    if (!direction_valid)
      throw std::runtime_error("Failed to sample valid direction after max attempts");

    // ----- STEP 2: feasible interval via single gemv + per-coord nonneg -----
    Ad = A_trans * direction;

    double t_min = -std::numeric_limits<double>::infinity();
    double t_max =  std::numeric_limits<double>::infinity();

    for (int j = 0; j < m; ++j) {
      const double Ad_j = Ad[j];
      const double Ax_j = Ax[j];
      if (std::abs(Ad_j) > 1e-12) {
        const double t_cand = (z_trans[j] - Ax_j) / Ad_j;
        if (Ad_j > 0) t_max = std::min(t_max, t_cand);
        else          t_min = std::max(t_min, t_cand);
      } else if (Ax_j > z_trans[j] + 1e-12) {
        throw std::runtime_error("Starting point not in feasible region");
      }
    }

    for (int i = 0; i < n_active; ++i) {
      const arma::uword k  = active_idx[i];
      const double v_k     = direction[k];
      const double x_k     = x0[k];
      if (std::abs(v_k) > 1e-12) {
        const double t_cand = -x_k / v_k;
        if (v_k < 0) t_max = std::min(t_max, t_cand);
        else         t_min = std::max(t_min, t_cand);
      } else if (x_k < -1e-12) {
        throw std::runtime_error("Starting point not in feasible region (nonneg)");
      }
    }

    if (t_min > t_max) continue;

    const double t_min_f = std::max(t_min, -bound_truncation);
    const double t_max_f = std::min(t_max,  bound_truncation);

    // ----- STEP 3: sample step size (kernel-specific) -----
    double alpha = 0.0;

    if (exponential) {
      const double lambda_rate = rho * direction_sum;
      if (lambda_rate <= 0)
        throw std::runtime_error("Internal error: lambda_rate <= 0 for exponential kernel");

      const double c = -lambda_rate;
      const double a = c * t_min_f;
      const double b = c * t_max_f;

      double u_rand = rng.runif();
      u_rand = std::max(1e-15, std::min(1.0 - 1e-15, u_rand));

      if (std::abs(b - a) < 1e-12) {
        alpha = t_min_f + u_rand * (t_max_f - t_min_f);
      } else if (b - a > 700) {
        alpha = (1.0 / c) * (b + std::log(u_rand));
      } else if (b - a < -700) {
        alpha = (1.0 / c) * (a + std::log(1.0 - u_rand));
      } else {
        alpha = (1.0 / c) * (a + std::log(1.0 - u_rand * (1.0 - std::exp(b - a))));
      }
    } else {
      double a_coeff = 0.0, b_coeff = 0.0;
      for (int i = 0; i < n_active; ++i) {
        const arma::uword k = active_idx[i];
        const double v_k = direction[k];
        const double w_k = x0[k];
        a_coeff += v_k * v_k;
        b_coeff += w_k * v_k;
      }
      a_coeff *= 0.5;

      if (a_coeff < 1e-15) {
        alpha = t_min_f + rng.runif() * (t_max_f - t_min_f);
      } else {
        const double mu_g = -b_coeff / (2.0 * a_coeff);
        const double sd_g = std::sqrt(1.0 / (2.0 * a_coeff));
        alpha = sample_truncnorm_icdf(rng, mu_g, sd_g, t_min_f, t_max_f);
      }
    }

    // ----- STEP 4: update position + Ax incrementally -----
    x0 += alpha * direction;
    Ax += alpha * Ad;
  }

  return x0;
}

void validate_common(const arma::mat& A, const arma::vec& greater_equal,
                     int d, const arma::vec& kappa, const std::string& kernel) {
  if ((int)kappa.n_elem != d)
    Rcpp::stop("Length of kappa must equal dimension of x0");
  if (greater_equal.n_elem != A.n_rows)
    Rcpp::stop("Length of greater_equal must equal number of rows of A");
  if (kernel != "exponential" && kernel != "half_gaussian")
    Rcpp::stop("kernel must be 'exponential' or 'half_gaussian'");
}

} // namespace

// [[Rcpp::export]]
arma::vec hit_and_run_cpp(
    const arma::mat& A,
    const arma::vec& z,
    const arma::vec& greater_equal,
    arma::vec x0,
    int n_iter,
    const arma::vec& kappa,
    double rho = 1.0,
    std::string kernel = "exponential",
    int max_dir_tries = 1000,
    double bound_truncation = 1e5
) {
  validate_common(A, greater_equal, x0.n_elem, kappa, kernel);
  if (arma::accu(kappa) < 1)
    Rcpp::stop("kappa must have at least one non-zero element");

  uint64_t seed = (uint64_t)(R::runif(0, 1) * 9007199254740992.0);
  Rng rng(seed);
  return hit_and_run_chain(A, z, greater_equal, x0, n_iter, kappa,
                           rho, kernel, max_dir_tries, bound_truncation, rng);
}

// [[Rcpp::export]]
arma::mat loop_hit_and_run_cpp(
    const arma::mat& A,
    const arma::mat& Z,
    const arma::mat& Greater_equal,
    const arma::mat& X0,
    const arma::mat& Kappa,
    int n_iter,
    double rho = 1.0,
    std::string kernel = "exponential",
    int max_dir_tries = 1000,
    double bound_truncation = 1e5,
    int n_threads = 1
) {
  const int n = X0.n_rows;
  const int d = X0.n_cols;

  if ((int)Z.n_rows != n || (int)Greater_equal.n_rows != n || (int)Kappa.n_rows != n)
    Rcpp::stop("Z, Greater_equal, and Kappa must have same number of rows as X0");
  if ((int)Kappa.n_cols != d)
    Rcpp::stop("Kappa must have same number of columns as X0");
  if ((int)A.n_cols != d)
    Rcpp::stop("ncol(A) must equal ncol(X0)");
  if (Z.n_cols != A.n_rows || Greater_equal.n_cols != A.n_rows)
    Rcpp::stop("Z and Greater_equal must have one column per row of A");
  if (kernel != "exponential" && kernel != "half_gaussian")
    Rcpp::stop("kernel must be 'exponential' or 'half_gaussian'");

  // One RNG stream per chain, seeded from R's RNG on the main thread, so
  // results do not depend on n_threads.
  std::vector<uint64_t> seeds(n);
  for (int i = 0; i < n; ++i) {
    seeds[i] = (uint64_t)(R::runif(0, 1) * 9007199254740992.0);
  }

  arma::mat final_samples(n, d);
  std::string err_msg;
  bool failed = false;

  if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n; ++i) {
    if (failed) continue;
    try {
      arma::vec z_i     = Z.row(i).t();
      arma::vec ge_i    = Greater_equal.row(i).t();
      arma::vec x0_i    = X0.row(i).t();
      arma::vec kappa_i = Kappa.row(i).t();

      if (arma::accu(kappa_i) == 0) {
        final_samples.row(i) = x0_i.t();
      } else {
        Rng rng(seeds[i]);
        arma::vec sample_i = hit_and_run_chain(
          A, z_i, ge_i, x0_i, n_iter, kappa_i,
          rho, kernel, max_dir_tries, bound_truncation, rng
        );
        final_samples.row(i) = sample_i.t();
      }
    } catch (const std::exception& e) {
#ifdef _OPENMP
      #pragma omp critical
#endif
      {
        failed = true;
        err_msg = e.what();
      }
    }
  }

  if (failed) Rcpp::stop(err_msg);

  return final_samples;
}

// [[Rcpp::export]]
Rcpp::IntegerVector check_feasible_dual_cpp(
    const arma::mat& A,
    const arma::mat& Z,
    const arma::mat& Greater_equal,
    const arma::mat& X0,
    int n_threads = 1
) {
  const int n = X0.n_rows;
  const int d = X0.n_cols;
  const int m = A.n_rows;

  Rcpp::IntegerVector feas(n);
  int* feas_ptr = INTEGER(feas);

  if (n_threads < 1) n_threads = 1;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int i = 0; i < n; ++i) {
    bool ok = true;
    arma::rowvec x0_i = X0.row(i);

    for (int k = 0; k < d; ++k) {
      if (x0_i[k] < 0) { ok = false; break; }
    }

    if (ok) {
      for (int j = 0; j < m; ++j) {
        const double dotVal = arma::dot(A.row(j), x0_i);
        const double z_val  = Z(i, j);
        const double ge_val = Greater_equal(i, j);
        if (ge_val > 0.5) {
          if (dotVal < z_val - 1e-12) { ok = false; break; }
        } else {
          if (dotVal > z_val + 1e-12) { ok = false; break; }
        }
      }
    }
    feas_ptr[i] = ok ? 1 : 0;
  }
  return feas;
}
