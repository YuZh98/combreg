# MCMC and regression diagnostics for crr_fit objects.

#' Effective sample sizes
#'
#' Effective sample size per coefficient, computed from post-warmup thinned
#' draws and summed across chains (via [coda::effectiveSize()]).
#'
#' @param object A `crr_fit` object.
#'
#' @return Named numeric vector of length `p * d`.
#'
#' @export
crr_ess <- function(object) {
  stopifnot(inherits(object, "crr_fit"))
  coda::effectiveSize(as.mcmc(object))
}

# Split-Rhat (Gelman et al., BDA3): each chain is halved, so the statistic is
# defined even for a single chain.
split_rhat <- function(mat) {
  half <- floor(nrow(mat) / 2)
  if (half < 2) return(NA_real_)
  splits <- cbind(mat[seq_len(half), , drop = FALSE],
                  mat[seq_len(half) + half, , drop = FALSE])
  means <- colMeans(splits)
  vars <- apply(splits, 2, stats::var)
  W <- mean(vars)
  B <- half * stats::var(means)
  if (W == 0) return(NA_real_)
  sqrt(((half - 1) / half * W + B / half) / W)
}

#' Split-Rhat convergence diagnostics
#'
#' Potential scale reduction factor per coefficient, computed on split
#' half-chains so it is defined even for `chains = 1`.
#'
#' @param object A `crr_fit` object.
#'
#' @return Named numeric vector of length `p * d`.
#'
#' @export
crr_rhat <- function(object) {
  stopifnot(inherits(object, "crr_fit"))
  dr <- kept_draws(object)
  out <- apply(dr, 3, split_rhat)
  names(out) <- dimnames(dr)[[3]]
  out
}

#' MCMC and regression diagnostics report
#'
#' Assembles a structured diagnostics report for a fitted model: per-parameter
#' posterior summaries with effective sample sizes, sampling efficiency
#' (ESS per second), split-Rhat convergence diagnostics, MH acceptance rates,
#' and in-sample regression fit (exact-match and per-coordinate accuracy of
#' the fitted responses). Potential problems (low ESS, high Rhat, low
#' acceptance) are collected as warnings and flagged by `print()`.
#'
#' @param object A `crr_fit` object.
#' @param prob Central posterior interval probability for the summary table.
#'
#' @return An object of class `crr_diagnostics`: a list with elements
#'   `table` (per-parameter data frame with columns `mean`, `sd`, interval
#'   bounds, `ess`, `ess_per_sec`, `rhat`), `fit_stats` (exact-match rate and
#'   per-coordinate accuracy, or `NULL` if fitted responses are unavailable),
#'   `warnings` (character vector), and the sampler configuration.
#'
#' @examples
#' con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
#' sim <- simulate_crr(n = 50, p = 2, constraints = con, seed = 1)
#' fit <- crr(sim$Y, sim$X, con, n_iter = 200, warmup = 100, seed = 1)
#' crr_diagnostics(fit)
#'
#' @export
crr_diagnostics <- function(object, prob = 0.95) {
  stopifnot(inherits(object, "crr_fit"), prob > 0, prob < 1)

  dr <- kept_draws(object)
  flat <- matrix(dr, prod(dim(dr)[1:2]), dim(dr)[3],
                 dimnames = list(NULL, dimnames(dr)[[3]]))
  alpha <- (1 - prob) / 2
  total_time <- sum(object$timing)
  ess <- crr_ess(object)
  rhat <- crr_rhat(object)

  tab <- data.frame(
    parameter = colnames(flat),
    mean = colMeans(flat),
    sd = apply(flat, 2, stats::sd),
    lower = apply(flat, 2, stats::quantile, alpha),
    upper = apply(flat, 2, stats::quantile, 1 - alpha),
    ess = ess,
    ess_per_sec = if (total_time > 0) ess / total_time else NA_real_,
    rhat = rhat,
    row.names = NULL
  )

  fit_stats <- tryCatch({
    Y_hat <- stats::fitted(object)
    list(
      exact_match = mean(rowSums(object$Y != Y_hat) == 0),
      coordinate_accuracy = colMeans(object$Y == Y_hat)
    )
  }, error = function(e) NULL)

  warnings <- character()
  bad_rhat <- tab$parameter[!is.na(tab$rhat) & tab$rhat > 1.05]
  if (length(bad_rhat)) {
    warnings <- c(warnings, paste0(
      "split-Rhat > 1.05 for ", length(bad_rhat), " parameter(s): ",
      paste(utils::head(bad_rhat, 5), collapse = ", "),
      if (length(bad_rhat) > 5) ", ..."
    ))
  }
  low_ess <- tab$parameter[tab$ess < 100]
  if (length(low_ess)) {
    warnings <- c(warnings, paste0(
      "ESS < 100 for ", length(low_ess), " parameter(s): ",
      paste(utils::head(low_ess, 5), collapse = ", "),
      if (length(low_ess) > 5) ", ..."
    ))
  }
  if (object$method == "mhwg" && mean(object$accept_rate) < 0.05) {
    warnings <- c(warnings, paste0(
      "mean MH acceptance rate is low (",
      format(mean(object$accept_rate), digits = 3),
      "); consider more iterations or a smaller zeta_block"
    ))
  }

  structure(
    list(
      table = tab,
      fit_stats = fit_stats,
      warnings = warnings,
      prob = prob,
      method = object$method, kernel = object$kernel,
      n = object$n, p = object$p, d = object$d,
      chains = object$chains, n_iter = object$n_iter,
      warmup = object$warmup, thin = object$thin,
      n_kept = prod(dim(dr)[1:2]),
      accept_rate = object$accept_rate,
      timing = object$timing
    ),
    class = "crr_diagnostics"
  )
}

#' @export
print.crr_diagnostics <- function(x, digits = 3, max_rows = 20, ...) {
  rule <- function(title) {
    cat(title, " ", strrep("-", max(1, 72 - nchar(title) - 1)), "\n", sep = "")
  }

  rule("Model")
  cat("  method ", x$method,
      if (!is.na(x$kernel)) paste0(" (", x$kernel, " dual kernel)"),
      ": n = ", x$n, ", p = ", x$p, ", d = ", x$d, "\n", sep = "")

  rule("Sampling")
  cat("  ", x$chains, " chain(s), ", x$n_iter, " iterations (warmup ",
      x$warmup, ", thin ", x$thin, "); ", x$n_kept,
      " retained draws\n", sep = "")
  cat("  total time ", format(sum(x$timing), digits = digits), "s", sep = "")
  if (x$method == "mhwg") {
    cat(", mean MH acceptance ",
        format(mean(x$accept_rate), digits = digits), sep = "")
  }
  cat("\n")

  rule("Coefficients (posterior summary + MCMC diagnostics)")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], function(v) signif(v, digits))
  names(tab)[names(tab) == "lower"] <- paste0("q", 100 * (1 - x$prob) / 2)
  names(tab)[names(tab) == "upper"] <- paste0("q", 100 * (1 + x$prob) / 2)
  if (nrow(tab) > max_rows) {
    print(utils::head(tab, max_rows), row.names = FALSE)
    cat("  ... ", nrow(tab) - max_rows, " more row(s); see $table\n", sep = "")
  } else {
    print(tab, row.names = FALSE)
  }
  cat("  worst split-Rhat ", format(max(x$table$rhat, na.rm = TRUE),
                                    digits = digits),
      "; min ESS ", format(min(x$table$ess), digits = digits),
      "; min ESS/sec ", format(min(x$table$ess_per_sec), digits = digits),
      "\n", sep = "")

  rule("Regression fit (training data)")
  if (is.null(x$fit_stats)) {
    cat("  fitted responses unavailable (requires the 'lpSolve' package",
        " for constrained fits)\n", sep = "")
  } else {
    cat("  exact-match rate: ",
        format(x$fit_stats$exact_match, digits = digits),
        "  (fraction of rows with y_i == y_hat_i)\n", sep = "")
    cat("  per-coordinate accuracy: ",
        paste(format(x$fit_stats$coordinate_accuracy, digits = digits),
              collapse = " "),
        "\n", sep = "")
  }

  rule("Warnings")
  if (length(x$warnings)) {
    for (w in x$warnings) cat("  - ", w, "\n", sep = "")
  } else {
    cat("  none\n")
  }
  invisible(x)
}

#' Fitted responses
#'
#' Posterior-point-estimate responses on the training covariates:
#' `predict(object, type = "response")`.
#'
#' @param object A `crr_fit` object.
#' @param ... Unused.
#'
#' @return Binary matrix (`n` x `d`).
#'
#' @export
fitted.crr_fit <- function(object, ...) {
  predict(object, type = "response")
}

#' Response residuals
#'
#' Training residuals `Y - fitted(object)`, entries in \{-1, 0, 1\}.
#'
#' @param object A `crr_fit` object.
#' @param ... Unused.
#'
#' @return Integer matrix (`n` x `d`).
#'
#' @export
residuals.crr_fit <- function(object, ...) {
  object$Y - fitted(object)
}
