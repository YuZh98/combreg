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
    if (identical(x$control$zeta_block, "adaptive") &&
        !is.null(x$zeta_block_tuned)) {
      tuned <- x$zeta_block_tuned
      cat("Adaptive zeta_block (target accept ",
          format(x$control$zeta_target_accept, digits = 2), "): ",
          if (length(unique(tuned)) == 1L) tuned[1]
          else paste(tuned, collapse = ", "), "\n", sep = "")
    }
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

#' Diagnostic plots for a crr fit
#'
#' All plots use post-warmup thinned draws.
#'
#' - `"trace"`: per-parameter trace plots, chains overlaid.
#' - `"acf"`: per-parameter autocorrelation functions (chains pooled).
#' - `"violin"`: violin plot of the marginal posterior of each coefficient,
#'   with median and central 95% interval.
#' - `"ess"`: effective sample size per coefficient, with a dashed reference
#'   line at 100.
#' - `"ess_time"`: effective sample size per second of sampling time ---
#'   the efficiency measure used for method comparison in the paper.
#' - `"residual"`: heat map of the training residuals `Y - fitted(x)`
#'   (entries in \{-1, 0, 1\}); rows are observations, columns response
#'   coordinates.
#' - `"beta_diff"`: heat map of the coefficient estimation error
#'   `coef(x) - beta` against a known true coefficient matrix (supplied via
#'   `beta`, e.g. from [simulate_crr()]), with the error printed in each
#'   cell.
#'
#' @param x A `crr_fit` object.
#' @param pars Character vector of parameter names (e.g. `"beta[1,2]"`).
#'   Defaults to the first four parameters for `"trace"`/`"acf"` and all
#'   parameters otherwise. Ignored for `"residual"` and `"beta_diff"`.
#' @param type Plot type; see Details.
#' @param beta True coefficient matrix (`p` x `d`); required for
#'   `type = "beta_diff"`.
#' @param ... Passed to the underlying base plotting functions.
#'
#' @return `x`, invisibly.
#'
#' @export
plot.crr_fit <- function(x, pars = NULL,
                         type = c("trace", "acf", "violin", "ess",
                                  "ess_time", "residual", "beta_diff"),
                         beta = NULL,
                         ...) {
  type <- match.arg(type)
  dr <- kept_draws(x)
  all_pars <- dimnames(dr)[[3]]
  if (is.null(pars)) {
    pars <- if (type %in% c("trace", "acf")) {
      all_pars[seq_len(min(4, length(all_pars)))]
    } else {
      all_pars
    }
  }
  bad <- setdiff(pars, all_pars)
  if (length(bad)) stop("Unknown parameters: ", paste(bad, collapse = ", "),
                        call. = FALSE)

  switch(type,
    trace = ,
    acf = plot_panels(dr, pars, all_pars, type, ...),
    violin = plot_violin(dr, pars, all_pars, ...),
    ess = plot_ess(x, pars, per_second = FALSE, ...),
    ess_time = plot_ess(x, pars, per_second = TRUE, ...),
    residual = plot_residual(x, ...),
    beta_diff = plot_beta_diff(x, beta, ...)
  )
  invisible(x)
}

plot_panels <- function(dr, pars, all_pars, type, ...) {
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
}

plot_violin <- function(dr, pars, all_pars, ...) {
  k <- length(pars)
  idx <- match(pars, all_pars)
  vals <- lapply(idx, function(i) as.vector(dr[, , i]))

  graphics::plot(NA, xlim = c(0.5, k + 0.5), ylim = range(unlist(vals)),
                 xaxt = "n", xlab = "", ylab = "coefficient",
                 main = "Posterior distributions", ...)
  graphics::axis(1, at = seq_len(k), labels = pars, las = 2, cex.axis = 0.8)
  graphics::abline(h = 0, lty = 3, col = "grey60")

  for (i in seq_len(k)) {
    v <- vals[[i]]
    den <- stats::density(v)
    w <- 0.4 * den$y / max(den$y)
    graphics::polygon(c(i - w, rev(i + w)), c(den$x, rev(den$x)),
                      col = "grey85", border = "grey40")
    q <- stats::quantile(v, c(0.025, 0.5, 0.975))
    graphics::segments(i, q[1], i, q[3], lwd = 2)
    graphics::points(i, q[2], pch = 19)
  }
}

plot_ess <- function(x, pars, per_second, ...) {
  ess <- crr_ess(x)[pars]
  ylab <- "effective sample size"
  if (per_second) {
    total_time <- sum(x$timing)
    if (total_time <= 0) stop("No recorded sampling time", call. = FALSE)
    ess <- ess / total_time
    ylab <- "effective sample size per second"
  }
  mid <- graphics::barplot(ess, las = 2, cex.names = 0.8, ylab = ylab,
                           main = if (per_second) "Sampling efficiency"
                                  else "Effective sample sizes",
                           ...)
  if (!per_second) graphics::abline(h = 100, lty = 2, col = "red")
  invisible(mid)
}

plot_beta_diff <- function(x, beta, ...) {
  if (is.null(beta)) {
    stop("type = \"beta_diff\" requires the true coefficient matrix ",
         "via the beta argument", call. = FALSE)
  }
  beta <- as.matrix(beta)
  if (nrow(beta) != x$p || ncol(beta) != x$d) {
    stop("beta must be a ", x$p, " x ", x$d, " matrix", call. = FALSE)
  }
  err <- coef(x) - beta
  lim <- max(abs(err), 1e-12)
  pal <- grDevices::colorRampPalette(c("#2166AC", "grey95", "#B2182B"))(64)

  graphics::image(seq_len(x$d), seq_len(x$p), t(err),
                  zlim = c(-lim, lim), col = pal,
                  xlab = "response coordinate j", ylab = "covariate k",
                  main = "Coefficient error  beta_hat - beta",
                  axes = FALSE, ...)
  graphics::axis(1, at = seq_len(x$d))
  graphics::axis(2, at = seq_len(x$p))
  graphics::box()
  for (k in seq_len(x$p)) {
    for (j in seq_len(x$d)) {
      graphics::text(j, k, signif(err[k, j], 2), cex = 0.8)
    }
  }
}

plot_residual <- function(x, ...) {
  res <- stats::residuals(x)
  n <- nrow(res)
  d <- ncol(res)
  graphics::image(seq_len(d), seq_len(n), t(res),
                  breaks = c(-1.5, -0.5, 0.5, 1.5),
                  col = c("#2166AC", "grey95", "#B2182B"),
                  xlab = "response coordinate", ylab = "observation",
                  main = "Residuals  y - y_hat", axes = FALSE, ...)
  graphics::axis(1, at = seq_len(d))
  graphics::axis(2)
  graphics::box()
  graphics::legend("topright", legend = c("-1", "0", "+1"),
                   fill = c("#2166AC", "grey95", "#B2182B"),
                   bg = "white", cex = 0.8)
}

#' Posterior predictions
#'
#' @param object A `crr_fit` object.
#' @param newdata Covariate matrix (`n_new` x `p`); defaults to the training
#'   covariates.
#' @param type `"utility"` returns posterior-mean latent utilities
#'   \eqn{X \hat\beta}; `"response"` additionally maps each utility row to a
#'   feasible response --- the constrained maximizer via integer programming
#'   for constrained fits (requires the 'lpSolve' package), or the
#'   coordinatewise sign indicator \eqn{1\{x^\top\hat\beta > 0\}} for
#'   `method = "unconstrained"` fits.
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
    return((utility > 0) * 1)
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
