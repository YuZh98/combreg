# Layer-1 sampling primitives. These are the building blocks used by crr()
# and are exported so that models with non-linear mean structure (e.g. the
# waterfowl matching analysis) can compose them with custom coefficient
# updates without duplicating any constraint logic.

#' Initialize dual certificates
#'
#' Computes the active-constraint mask `active` (\eqn{1\{A y_i = b\}}) and an
#' initial dual matrix `U` with independent Exp(1) entries on the active
#' support.
#'
#' @param constraints A [crr_constraints] object.
#' @param Y Response matrix (`n` x `d`).
#'
#' @return List with `U` (`n` x `m`) and `active` (`n` x `m` binary mask).
#'
#' @keywords internal
#' @export
init_dual <- function(constraints, Y) {
  stopifnot(inherits(constraints, "crr_constraints"))
  Y <- as.matrix(Y)
  n <- nrow(Y)
  m <- constraints$m
  active <- (Y %*% t(constraints$A) == rep(constraints$b, each = n)) * 1
  U <- matrix(rexp(n * m), n, m) * active
  list(U = U, active = active)
}

#' Update dual certificates given latent utilities
#'
#' One Gibbs update of \eqn{U \mid \zeta, y}: for each observation, runs
#' hit-and-run over the dual-certificate polytope
#' \eqn{\mathcal{U}(y_i, \zeta_i) \cap \{u \ge 0\}} targeting the chosen
#' kernel, and returns the endpoint.
#'
#' @param constraints A [crr_constraints] object.
#' @param zeta Latent utility matrix (`n` x `d`).
#' @param Y Response matrix (`n` x `d`).
#' @param U Current dual matrix (`n` x `m`); must be feasible.
#' @param active Active-constraint mask from [init_dual()].
#' @param kernel `"exponential"` or `"half_gaussian"`.
#' @param control A [crr_control] object.
#'
#' @return Updated dual matrix (`n` x `m`).
#'
#' @keywords internal
#' @export
sample_dual <- function(constraints, zeta, Y, U, active,
                        kernel = "exponential",
                        control = crr_control()) {
  loop_hit_and_run_cpp(
    t(constraints$A), zeta, 1 - Y, U, active,
    n_iter = control$n_iter_hit_and_run,
    rho = control$rho,
    kernel = kernel,
    max_dir_tries = control$max_dir_tries,
    bound_truncation = control$bound_truncation,
    n_threads = control$n_threads
  )
}

#' Draw latent utilities from their truncated-normal full conditional
#'
#' For each coordinate `j` in `subset`, draws
#' \eqn{\zeta_{\cdot j} \sim TN(\mu_{\cdot j}, 1)} truncated below at
#' \eqn{(UA)_{\cdot j}} where \eqn{y = 1} and above where \eqn{y = 0}.
#'
#' @param zeta Current latent utility matrix (`n` x `d`); returned unchanged
#'   outside `subset`.
#' @param Mu Mean matrix (`n` x `d`), e.g. `X %*% beta`.
#' @param Y Response matrix (`n` x `d`).
#' @param UA Dual bounds `U %*% A` (`n` x `d`).
#' @param subset Integer vector of coordinates to update.
#'
#' @return Matrix (`n` x `d`).
#'
#' @keywords internal
#' @export
draw_utility <- function(zeta, Mu, Y, UA, subset = seq_len(ncol(zeta))) {
  n <- nrow(zeta)
  for (j in subset) {
    lower <- ifelse(Y[, j] == 1, UA[, j], -Inf)
    upper <- ifelse(Y[, j] == 1, Inf, UA[, j])
    zeta[, j] <- truncnorm::rtruncnorm(n, a = lower, b = upper,
                                       mean = Mu[, j], sd = 1)
  }
  zeta
}

#' Check dual-certificate feasibility
#'
#' Verifies, per observation, that a dual matrix lies in
#' \eqn{\mathcal{U}(y_i, \zeta_i) \cap \{u \ge 0\}}.
#'
#' @inheritParams sample_dual
#'
#' @return Integer 0/1 vector of length `n`.
#'
#' @keywords internal
#' @export
dual_feasible <- function(constraints, zeta, Y, U,
                          control = crr_control()) {
  check_feasible_dual_cpp(t(constraints$A), zeta, 1 - Y, U,
                          n_threads = control$n_threads)
}

#' Metropolis-Hastings update of the latent utilities
#'
#' One MH step of \eqn{\zeta \mid U, \beta, y}: proposes new utilities for a
#' random block of coordinates from the truncated-normal proposal, forms the
#' coordinatewise envelope \eqn{\tilde\zeta} (max where \eqn{y=1}, min where
#' \eqn{y=0}), runs hit-and-run under \eqn{\tilde\zeta} to obtain a dual
#' certificate \eqn{U^\star}, and accepts each observation's proposal iff its
#' \eqn{u_i^\star} is feasible for the current \eqn{\zeta_i}.
#'
#' @inheritParams sample_dual
#' @param Mu Mean matrix (`n` x `d`), e.g. `X %*% beta`.
#' @param block Number of coordinates to update this step. When `NULL`
#'   (default), it is resolved from `control$zeta_block`: the fixed integer if
#'   set, otherwise `min(d, 100)` (the adaptive controller in [crr()] supplies
#'   an explicit value each step, so this fallback only applies to direct
#'   callers under an adaptive control).
#'
#' @return List with `zeta` (updated matrix) and `accept` (0/1 vector,
#'   per-observation acceptance).
#'
#' @keywords internal
#' @export
sample_utility <- function(constraints, zeta, Y, U, active, Mu,
                           kernel = "exponential",
                           control = crr_control(), block = NULL) {
  d <- ncol(zeta)
  UA <- U %*% constraints$A

  if (is.null(block)) {
    block <- if (is.numeric(control$zeta_block)) control$zeta_block else 100
  }
  block <- max(1L, min(d, as.integer(block)))
  subset <- sample.int(d, block, replace = FALSE)
  zeta_new <- draw_utility(zeta, Mu, Y, UA, subset)

  zeta_tilde <- (Y > 0.5) * pmax(zeta, zeta_new) + (Y < 0.5) * pmin(zeta, zeta_new)
  # The certificate governs the acceptance probability, so it gets its own
  # sweep count; falls back to n_iter_hit_and_run for controls built elsewhere.
  control_mh <- control
  if (!is.null(control$n_iter_hit_and_run_mh)) {
    control_mh$n_iter_hit_and_run <- control$n_iter_hit_and_run_mh
  }
  U_star <- sample_dual(constraints, zeta_tilde, Y, U, active, kernel,
                        control_mh)

  accept <- dual_feasible(constraints, zeta, Y, U_star, control)
  list(
    zeta = zeta_new * accept + zeta * (1 - accept),
    accept = accept
  )
}

#' Precompute quantities for the conjugate coefficient update
#'
#' @param X Covariate matrix (`n` x `p`).
#' @param prior A [crr_prior] object.
#'
#' @return List with posterior covariance `V` and its lower Cholesky factor
#'   `L`, for use in [update_coef()].
#'
#' @keywords internal
#' @export
coef_precompute <- function(X, prior = crr_prior()) {
  stopifnot(inherits(prior, "crr_prior"))
  p <- ncol(X)
  V <- solve(diag(1 / prior$sd^2, p) + crossprod(X))
  list(V = V, L = t(chol(V)))
}

#' Conjugate Gaussian update of the regression coefficients
#'
#' Draws \eqn{\beta_{\cdot j} \sim N(V X^\top \zeta_{\cdot j}, V)} jointly for
#' all response coordinates.
#'
#' @param X Covariate matrix (`n` x `p`).
#' @param zeta Latent utility matrix (`n` x `d`).
#' @param precomp Output of [coef_precompute()].
#'
#' @return Coefficient matrix (`p` x `d`).
#'
#' @keywords internal
#' @export
update_coef <- function(X, zeta, precomp) {
  p <- ncol(X)
  d <- ncol(zeta)
  mean_mat <- precomp$V %*% crossprod(X, zeta)
  mean_mat + precomp$L %*% matrix(rnorm(p * d), p, d)
}
