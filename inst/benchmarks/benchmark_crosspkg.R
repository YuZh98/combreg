#!/usr/bin/env Rscript
# Cross-package benchmark: combreg vs. competing CRAN Bayesian samplers on
# combreg's problem setup (binary response vector on a totally-unimodular
# integral polytope, y in {0,1}^d s.t. A y <= b).
#
# This is a standalone research/benchmark asset, NOT package code. Run it with
#   Rscript inst/benchmarks/benchmark_crosspkg.R
# It loads the package with pkgload::load_all(), fits every available method on
# the SAME simulated data across a few seeds, and reports RMSE, 90% credible
# interval coverage, posterior-predictive feasibility, and sampling efficiency
# (min bulk-ESS / wall-clock second). Absent competitor packages are skipped.
#
# Methods:
#   1. combreg mhwg           - the proposed constrained sampler
#   2. combreg unconstrained  - independent Albert-Chib probit (bias baseline)
#   3. bayesm::rmvpGibbs      - multivariate probit that IGNORES the constraint
#   4. MNP::mnp               - multinomial probit, valid ONLY on the simplex
#
# Scenario A: general matching polytope (2x2 bipartite assignment) - methods 1-3
# Scenario B: simplex special case sum(y) <= 1 ("at most one") - methods 1-4

# load_all() walks up from cwd to the package root (run from anywhere in the repo).
suppressWarnings(suppressMessages(pkgload::load_all(quiet = TRUE)))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SMOKE <- TRUE   # TRUE: fast self-test (2 seeds, tiny iters, small n). FALSE: full run.

if (SMOKE) {
  SEEDS   <- 1:2
  N_OBS   <- 80
  N_ITER  <- 300
  WARMUP  <- 150
  S_PRED  <- 30    # posterior-predictive draws for the feasibility metric
} else {
  SEEDS   <- 1:5
  N_OBS   <- 300
  N_ITER  <- 2000
  WARMUP  <- 1000
  S_PRED  <- 100
}
P_COV <- 3        # number of covariates (no intercept: matches simulate_crr)

# ---------------------------------------------------------------------------
# Generic metric helpers (operate on a flat R x (p*d) posterior draw matrix
# whose columns are already mapped to combreg's beta parameterization, i.e.
# column order = as.vector(beta) = [beta[,1]; beta[,2]; ...]).
# ---------------------------------------------------------------------------

ess_bulk_min <- function(flat) {
  if (requireNamespace("posterior", quietly = TRUE)) {
    e <- apply(flat, 2, function(col) posterior::ess_bulk(as.numeric(col)))
  } else {
    e <- coda::effectiveSize(coda::mcmc(flat))
  }
  min(e[is.finite(e)])
}

metrics_from_draws <- function(flat, truth, time_sec) {
  pm <- colMeans(flat)
  q  <- apply(flat, 2, stats::quantile, probs = c(0.05, 0.95))
  list(
    rmse            = sqrt(mean((pm - truth)^2)),
    coverage90      = mean(q[1, ] <= truth & truth <= q[2, ]),
    min_ess_per_sec = ess_bulk_min(flat) / time_sec
  )
}

# Feasibility of posterior-predictive responses: draw beta from the posterior,
# simulate utilities, map to a response by `rule`, and report the fraction of
# response rows satisfying A y <= b. Constraint-ignoring rules expose infeasible
# posterior mass; combreg's constrained argmax is feasible by construction.
feasibility_rate <- function(flat, X, con, rule, S) {
  d <- con$d; pcov <- ncol(X); n <- nrow(X)
  idx <- sample.int(nrow(flat), min(S, nrow(flat)))
  ok <- 0; tot <- 0
  for (r in idx) {
    beta <- matrix(flat[r, ], pcov, d)
    zeta <- X %*% beta + matrix(stats::rnorm(n * d), n, d)
    Yr <- rule(zeta)
    ok  <- ok + sum(is_feasible(con, Yr))
    tot <- tot + n
  }
  ok / tot
}

rule_threshold <- function(zeta) (zeta > 0) * 1          # coordinatewise probit
rule_argmax_con <- function(con) function(zeta) {         # constrained maximizer
  t(apply(zeta, 1, ilp_argmax, constraints = con))
}

# ---------------------------------------------------------------------------
# Method runners. Each returns a one-row data.frame of metrics (or NULL to skip).
# ---------------------------------------------------------------------------

combreg_flat <- function(fit) {
  keep <- seq(fit$warmup + 1, fit$n_iter)
  dr <- fit$draws[keep, , , drop = FALSE]
  matrix(dr, prod(dim(dr)[1:2]), dim(dr)[3])
}

run_combreg <- function(method, Y, X, con, beta_true, seed, feasible_rule) {
  ctrl <- crr_control(n_iter_hit_and_run = if (SMOKE) 10 else 50)
  fit <- if (method == "mhwg") {
    crr(Y, X, con, method = "mhwg", n_iter = N_ITER, warmup = WARMUP,
        seed = seed, control = ctrl)
  } else {
    suppressWarnings(
      crr(Y, X, method = "unconstrained", n_iter = N_ITER, warmup = WARMUP,
          seed = seed))
  }
  flat <- combreg_flat(fit)
  truth <- as.vector(beta_true)
  m <- metrics_from_draws(flat, truth, sum(fit$timing))
  feas <- feasibility_rate(flat, X, con, feasible_rule, S_PRED)
  data.frame(method = paste0("combreg ", method),
             rmse = m$rmse, coverage90 = m$coverage90,
             feasibility = feas, min_ess_per_sec = m$min_ess_per_sec,
             time_sec = sum(fit$timing), note = "", stringsAsFactors = FALSE)
}

# bayesm multivariate probit. PARAMETERIZATION MAPPING (the crux):
#   bayesm model: w_ij = x_i' beta_j + e_i, e ~ N(0, Sigma) with UNRESTRICTED
#   Sigma, so (beta_j, Sigma) are identified only up to a per-equation scale.
#   The sign pattern y = 1{w>0} is invariant to scaling row j by 1/sqrt(Sigma_jj),
#   so the identified coefficient is beta_j / sqrt(Sigma_jj) -- exactly combreg's
#   unit-variance-noise scale (eps ~ N(0, I_d) => unit variance per coordinate).
#   Data layout: y and X are stacked observation-major, p=d rows per subject;
#   row (i-1)*d+j carries x_i in coordinate-j's coefficient block, so the mapped
#   coefficient vector already matches combreg's as.vector(beta) column order.
run_bayesm <- function(Y, X, con, beta_true, seed) {
  if (!requireNamespace("bayesm", quietly = TRUE)) {
    message("  [skip] bayesm not installed")
    return(NULL)
  }
  d <- con$d; pcov <- ncol(X); n <- nrow(X); k <- d * pcov
  Xb <- matrix(0, n * d, k)
  for (i in seq_len(n)) for (j in seq_len(d)) {
    Xb[(i - 1) * d + j, (j - 1) * pcov + seq_len(pcov)] <- X[i, ]
  }
  yb <- as.vector(t(Y))                       # y[(i-1)*d + j] = Y[i,j]
  set.seed(seed)
  fit <- NULL
  t0 <- proc.time()[["elapsed"]]
  # capture.output evaluates in this frame, so `fit` is assigned here; it also
  # swallows rmvpGibbs's progress printing.
  utils::capture.output(
    fit <- tryCatch(
      bayesm::rmvpGibbs(Data = list(p = d, y = yb, X = Xb),
                        Mcmc = list(R = N_ITER, keep = 1L, nprint = 0)),
      error = function(e) e))
  time_sec <- proc.time()[["elapsed"]] - t0
  if (inherits(fit, "error")) {
    message("  [skip] bayesm::rmvpGibbs failed: ", conditionMessage(fit))
    return(NULL)
  }
  bd <- fit$betadraw; sd <- fit$sigmadraw
  drop <- seq_len(floor(nrow(bd) * WARMUP / N_ITER))
  bd <- bd[-drop, , drop = FALSE]; sd <- sd[-drop, , drop = FALSE]
  # Rescale each coordinate block j by 1/sqrt(Sigma_jj) -> combreg scale.
  flat <- t(vapply(seq_len(nrow(bd)), function(r) {
    S <- matrix(sd[r, ], d, d)
    v <- bd[r, ]
    for (j in seq_len(d)) {
      idx <- (j - 1) * pcov + seq_len(pcov)
      v[idx] <- v[idx] / sqrt(S[j, j])
    }
    v
  }, numeric(k)))
  truth <- as.vector(beta_true)
  m <- metrics_from_draws(flat, truth, time_sec)
  # Feasibility: identified model has unit-variance (correlated) noise; iid unit
  # noise on the mapped beta is faithful for the marginal feasibility check.
  feas <- feasibility_rate(flat, X, con, rule_threshold, S_PRED)
  data.frame(method = "bayesm::rmvpGibbs",
             rmse = m$rmse, coverage90 = m$coverage90,
             feasibility = feas, min_ess_per_sec = m$min_ess_per_sec,
             time_sec = time_sec,
             note = "MVP ignores constraint; beta/sqrt(Sigma_jj)",
             stringsAsFactors = FALSE)
}

# MNP multinomial probit on the simplex ("at most one" => choose among the d
# coordinates plus an outside good = "none", whose utility is the fixed base 0).
#   IDENTIFICATION MISMATCH: MNP works with base-relative utility DIFFERENCES and
#   fixes the differenced-error scale, whereas combreg uses absolute utilities
#   with unit-variance iid noise. There is no exact scale match, so we compare on
#   the identified STRUCTURE: the best least-squares scalar between estimated and
#   true coefficients is removed (reported as est_scale), and RMSE/coverage are
#   computed on the rescaled draws. base="none" has true beta 0, so MNP's
#   alternative-j coefficient targets beta[,j] up to that scale.
run_mnp <- function(Y, X, con, beta_true, seed) {
  if (!requireNamespace("MNP", quietly = TRUE)) {
    message("  [skip] MNP not installed")
    return(NULL)
  }
  d <- con$d; pcov <- ncol(X); n <- nrow(X)
  choice <- apply(Y, 1, function(row) if (sum(row) == 0) 0L else which(row == 1)[1])
  df <- data.frame(choice = factor(choice, levels = c(0, seq_len(d))), X)
  names(df) <- c("choice", paste0("x", seq_len(pcov)))
  form <- stats::as.formula(paste("choice ~", paste0("x", seq_len(pcov),
                                                      collapse = " + ")))
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    MNP::mnp(form, data = df, base = "0", n.draws = N_ITER, burnin = WARMUP,
             verbose = FALSE),
    error = function(e) e)
  if (inherits(fit, "error")) {
    message("  [skip] MNP::mnp failed: ", conditionMessage(fit))
    return(NULL)
  }
  time_sec <- proc.time()[["elapsed"]] - t0
  param <- fit$param
  # Coefficient columns: those not belonging to the covariance block (MNP names
  # covariance draws "cov(.)"). Keep only slope coefficients (drop intercepts).
  is_cov <- grepl("^cov", colnames(param))
  bcols <- colnames(param)[!is_cov]
  slope <- bcols[grepl(paste0("^(", paste0("x", seq_len(pcov), collapse = "|"),
                              "):"), bcols)]
  if (length(slope) == 0L) {
    message("  [skip] MNP: could not locate slope-coefficient columns")
    return(NULL)
  }
  draws <- param[, slope, drop = FALSE]
  # Map MNP draws to a (pcov*d)-vector in combreg's as.vector(beta) order. Column
  # names look like "x<k>:<alt>"; place coefficient (k, alt) into block `alt`.
  key <- do.call(rbind, lapply(strsplit(slope, ":"), function(s) {
    c(as.integer(sub("x", "", s[1])), as.integer(s[2]))
  }))
  target <- (key[, 2] - 1L) * pcov + key[, 1]   # index into as.vector(beta)
  flat <- matrix(0, nrow(draws), pcov * d)
  flat[, target] <- as.matrix(draws)
  truth <- as.vector(beta_true)
  # Remove the unknown identification scale via least squares on posterior means.
  pm <- colMeans(flat)
  est_scale <- sum(pm * truth) / sum(pm * pm)
  flat_s <- flat * est_scale
  m <- metrics_from_draws(flat_s, truth, time_sec)
  data.frame(method = "MNP::mnp",
             rmse = m$rmse, coverage90 = m$coverage90,
             feasibility = 1,           # choice model => exactly one 1 (sum<=1)
             min_ess_per_sec = m$min_ess_per_sec, time_sec = time_sec,
             note = sprintf("contrast-scale; est_scale=%.2f", est_scale),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Scenario drivers
# ---------------------------------------------------------------------------

run_scenario <- function(name, con, methods) {
  message("== Scenario ", name, " (d=", con$d, ") ==")
  feas_rule <- rule_argmax_con(con)
  rows <- list()
  for (seed in SEEDS) {
    sim <- simulate_crr(N_OBS, P_COV, con, seed = seed)
    per <- list()
    if ("mhwg" %in% methods)
      per[["mhwg"]] <- run_combreg("mhwg", sim$Y, sim$X, con, sim$beta, seed,
                                   feas_rule)
    if ("unconstrained" %in% methods)
      per[["unc"]] <- run_combreg("unconstrained", sim$Y, sim$X, con, sim$beta,
                                  seed, rule_threshold)
    if ("bayesm" %in% methods)
      per[["bayesm"]] <- run_bayesm(sim$Y, sim$X, con, sim$beta, seed)
    if ("mnp" %in% methods)
      per[["mnp"]] <- run_mnp(sim$Y, sim$X, con, sim$beta, seed)
    per <- Filter(Negate(is.null), per)
    if (length(per)) {
      dfp <- do.call(rbind, per)
      dfp$seed <- seed
      rows[[length(rows) + 1L]] <- dfp
    }
  }
  if (!length(rows)) return(NULL)
  long <- do.call(rbind, rows)
  long$scenario <- name
  long
}

# Scenario A: 2x2 bipartite assignment (matching-type). Variables y_{rc},
# ordered (11,12,21,22); each row and each column selected at most once.
A_assign <- rbind(c(1, 1, 0, 0),   # row 1
                  c(0, 0, 1, 1),   # row 2
                  c(1, 0, 1, 0),   # col 1
                  c(0, 1, 0, 1))   # col 2
con_A <- crr_constraints(A_assign, b = rep(1, 4))

# Scenario B: simplex special case, sum(y) <= 1 ("at most one").
con_B <- crr_constraints(matrix(1, 1, 4), b = 1)

resA <- run_scenario("A: matching polytope", con_A,
                     c("mhwg", "unconstrained", "bayesm"))
resB <- run_scenario("B: simplex (sum<=1)", con_B,
                     c("mhwg", "unconstrained", "bayesm", "mnp"))

long <- do.call(rbind, Filter(Negate(is.null), list(resA, resB)))

# ---------------------------------------------------------------------------
# Aggregate over seeds and report
# ---------------------------------------------------------------------------
agg <- stats::aggregate(
  cbind(rmse, coverage90, feasibility, min_ess_per_sec, time_sec) ~
    scenario + method, data = long, FUN = mean)
notes <- stats::aggregate(note ~ scenario + method, data = long,
                          FUN = function(x) x[nzchar(x)][1])
tab <- merge(agg, notes, by = c("scenario", "method"), all.x = TRUE)
tab$note[is.na(tab$note)] <- ""
tab <- tab[order(tab$scenario, tab$method), ]
num <- vapply(tab, is.numeric, logical(1))
tab[num] <- lapply(tab[num], function(v) round(v, 3))

cat("\n")
cat("Cross-package benchmark: combreg vs. CRAN Bayesian samplers\n")
cat("Averaged over ", length(SEEDS), " seed(s); n=", N_OBS, ", p=", P_COV,
    ", n_iter=", N_ITER, " (warmup ", WARMUP, ")",
    if (SMOKE) "  [SMOKE]" else "", "\n", sep = "")
cat("RMSE / coverage vs true beta; feasibility = P(pred. y satisfies A y<=b);",
    "eff = min bulk-ESS/sec\n")
print(tab, row.names = FALSE)
cat("\n")

out_csv <- file.path("inst", "benchmarks", "benchmark_crosspkg_results.csv")
if (!dir.exists(dirname(out_csv))) out_csv <- "benchmark_crosspkg_results.csv"
utils::write.csv(tab, out_csv, row.names = FALSE)
cat("Wrote ", normalizePath(out_csv), "\n", sep = "")
