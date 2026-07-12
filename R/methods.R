# S3 methods for crr_fit.

kept_draws <- function(object) {
  keep <- seq(object$warmup + 1, object$n_iter, by = object$thin)
  object$draws[keep, , , drop = FALSE]
}

#' @export
print.crr_fit <- function(x, ...) {
  cat("<crr_fit> method=", x$method,
      if (!is.na(x$kernel)) paste0(" kernel=", x$kernel),
      ": n=", x$n, ", p=", x$p, ", d=", x$d,
      ", chains=", x$chains, ", iterations=", x$n_iter,
      " (warmup ", x$warmup, ", thin ", x$thin, ")\n", sep = "")
  if (x$method == "mhwg") {
    cat("Mean MH acceptance rate: ",
        format(mean(x$accept_rate), digits = 3), "\n", sep = "")
  }
  cat("Total sampling time: ", format(sum(x$timing), digits = 3), "s\n",
      sep = "")
  invisible(x)
}

#' Posterior summary of a crr fit
#'
#' @param object A `crr_fit` object.
#' @param ... Unused.
#'
#' @return A data frame with posterior mean, sd, and 2.5/50/97.5 percent
#'   quantiles per coefficient, computed from post-warmup thinned draws
#'   pooled across chains.
#'
#' @export
summary.crr_fit <- function(object, ...) {
  dr <- kept_draws(object)
  flat <- matrix(dr, prod(dim(dr)[1:2]), dim(dr)[3],
                 dimnames = list(NULL, dimnames(dr)[[3]]))
  out <- data.frame(
    parameter = colnames(flat),
    mean = colMeans(flat),
    sd = apply(flat, 2, stats::sd),
    q2.5 = apply(flat, 2, stats::quantile, 0.025),
    q50 = apply(flat, 2, stats::quantile, 0.5),
    q97.5 = apply(flat, 2, stats::quantile, 0.975),
    row.names = NULL
  )
  class(out) <- c("summary.crr_fit", "data.frame")
  out
}

#' Posterior mean coefficients
#'
#' @param object A `crr_fit` object.
#' @param ... Unused.
#'
#' @return Posterior mean coefficient matrix (`p` x `d`).
#'
#' @export
coef.crr_fit <- function(object, ...) {
  dr <- kept_draws(object)
  flat <- matrix(dr, prod(dim(dr)[1:2]), dim(dr)[3])
  matrix(colMeans(flat), object$p, object$d)
}

#' Trace and autocorrelation plots
#'
#' Plots post-warmup traces (and optionally autocorrelations) for a subset of
#' coefficients.
#'
#' @param x A `crr_fit` object.
#' @param pars Character vector of parameter names (e.g. `"beta[1,2]"`), or
#'   `NULL` for the first four.
#' @param type `"trace"` or `"acf"`.
#' @param ... Passed to the base plotting functions.
#'
#' @return `x`, invisibly.
#'
#' @export
plot.crr_fit <- function(x, pars = NULL, type = c("trace", "acf"), ...) {
  type <- type[1]
  dr <- kept_draws(x)
  all_pars <- dimnames(dr)[[3]]
  if (is.null(pars)) pars <- all_pars[seq_len(min(4, length(all_pars)))]
  bad <- setdiff(pars, all_pars)
  if (length(bad)) stop("Unknown parameters: ", paste(bad, collapse = ", "),
                        call. = FALSE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  nr <- ceiling(length(pars) / 2)
  graphics::par(mfrow = c(nr, min(2, length(pars))))

  for (pn in pars) {
    idx <- match(pn, all_pars)
    if (type == "trace") {
      graphics::matplot(dr[, , idx], type = "l", lty = 1,
                        main = pn, xlab = "iteration (post-warmup)",
                        ylab = pn, ...)
    } else {
      stats::acf(as.vector(dr[, , idx]), main = paste("ACF:", pn), ...)
    }
  }
  invisible(x)
}

#' Posterior predictions
#'
#' @param object A `crr_fit` object.
#' @param newdata Covariate matrix (`n_new` x `p`); defaults to the training
#'   covariates.
#' @param type `"utility"` returns posterior-mean latent utilities
#'   \eqn{X \hat\beta}; `"response"` additionally maps each utility row to
#'   its constrained maximizer via integer programming (requires the
#'   'lpSolve' package and a fit with constraints).
#' @param ... Unused.
#'
#' @return Matrix (`n_new` x `d`) of utilities or feasible responses.
#'
#' @export
predict.crr_fit <- function(object, newdata = NULL,
                            type = c("utility", "response"), ...) {
  type <- match.arg(type)
  X_new <- if (is.null(newdata)) object$X else as.matrix(newdata)
  if (ncol(X_new) != object$p) {
    stop("newdata must have ", object$p, " columns", call. = FALSE)
  }
  utility <- X_new %*% coef(object)
  if (type == "utility") return(utility)

  if (is.null(object$constraints)) {
    stop("type = \"response\" requires a fit with constraints", call. = FALSE)
  }
  t(apply(utility, 1, ilp_argmax, constraints = object$constraints))
}

#' Convert to coda mcmc.list
#'
#' Post-warmup, thinned draws per chain.
#'
#' @param x A `crr_fit` object.
#' @param ... Unused.
#'
#' @return A [coda::mcmc.list].
#'
#' @importFrom coda as.mcmc
#' @export
as.mcmc.crr_fit <- function(x, ...) {
  dr <- kept_draws(x)
  coda::mcmc.list(lapply(seq_len(x$chains), function(ch) {
    coda::mcmc(dr[, ch, ], start = x$warmup + 1, thin = x$thin)
  }))
}

#' Convert to posterior draws_array
#'
#' @param x A `crr_fit` object.
#' @param ... Unused.
#'
#' @return A `posterior::draws_array` (requires the 'posterior' package).
#'
#' @exportS3Method posterior::as_draws
as_draws.crr_fit <- function(x, ...) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for as_draws()", call. = FALSE)
  }
  posterior::as_draws_array(kept_draws(x))
}
