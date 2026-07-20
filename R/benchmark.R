# Benchmark comparison across samplers on a common data set.

#' Benchmark samplers on a common data set
#'
#' Fits several methods to the same data and collects a comparison table:
#' wall-clock time, MH acceptance, effective sample sizes and sampling
#' efficiency (min ESS per second), in-sample fit, and --- when the true
#' coefficients are supplied, e.g. from [simulate_crr()] --- estimation error
#' (RMSE) and coverage of the central posterior intervals. This is the
#' package's canonical way to compare the paper's constrained sampler with
#' the unconstrained probit baseline (or any subset of registered methods)
#' without writing bespoke benchmarking code.
#'
#' @param Y Response matrix (`n` x `d`).
#' @param X Covariate matrix (`n` x `p`).
#' @param constraints A [crr_constraints] object (used by constrained
#'   methods only).
#' @param methods Character vector of methods accepted by [crr()].
#' @param beta True coefficient matrix (`p` x `d`), or `NULL`. When supplied,
#'   RMSE and interval coverage are reported.
#' @param prob Central posterior interval probability for coverage.
#' @param ... Passed on to [crr()] (e.g. `n_iter`, `warmup`, `chains`,
#'   `kernel`, `seed`, `control`). The same arguments --- including any seed ---
#'   are used for every method.
#'
#' @return An object of class `crr_benchmark`: a list with `table` (one row
#'   per method) and `fits` (the underlying `crr_fit` objects, named by
#'   method).
#'
#' @examples
#' if (requireNamespace("lpSolve", quietly = TRUE)) {
#'   con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
#'   sim <- simulate_crr(n = 50, p = 2, constraints = con, seed = 1)
#'   bm <- crr_benchmark(sim$Y, sim$X, con, beta = sim$beta,
#'                       n_iter = 200, warmup = 100, seed = 1)
#'   bm
#' }
#'
#' @export
crr_benchmark <- function(Y, X, constraints,
                          methods = c("mhwg", "unconstrained"),
                          beta = NULL, prob = 0.95, ...) {
  methods <- unique(methods)
  if (!is.null(beta)) beta <- as.matrix(beta)

  fits <- lapply(methods, function(mth) {
    if (mth == "mhwg") {
      crr(Y, X, constraints, method = mth, ...)
    } else {
      crr(Y, X, method = mth, ...)
    }
  })
  names(fits) <- methods

  rows <- lapply(methods, function(mth) {
    fit <- fits[[mth]]
    diag <- crr_diagnostics(fit, prob = prob, n_rep = 0)
    tab <- diag$table
    row <- data.frame(
      method = mth,
      time_sec = sum(fit$timing),
      accept_rate = mean(fit$accept_rate),
      min_ess = min(tab$ess),
      median_ess = stats::median(tab$ess),
      min_ess_per_sec = min(tab$ess_per_sec),
      exact_match = if (is.null(diag$fit_stats)) NA_real_ else
        diag$fit_stats$exact_match
    )
    if (!is.null(beta)) {
      truth <- as.vector(beta)
      row$rmse <- sqrt(mean((as.vector(coef(fit)) - truth)^2))
      row$coverage <- mean(tab$lower <= truth & truth <= tab$upper)
    }
    row
  })

  structure(
    list(table = do.call(rbind, rows), fits = fits, prob = prob),
    class = "crr_benchmark"
  )
}

#' @export
print.crr_benchmark <- function(x, digits = 3, ...) {
  cat("<crr_benchmark>: ", nrow(x$table), " method(s), ",
      100 * x$prob, "% intervals\n", sep = "")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], function(v) signif(v, digits))
  print(tab, row.names = FALSE)
  invisible(x)
}
