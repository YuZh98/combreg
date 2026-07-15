#' Bayesian combinatorial response regression
#'
#' Fits the latent-utility regression model
#' \deqn{\zeta_i = \beta^\top x_i + \varepsilon_i, \quad
#'       \varepsilon_i \sim N(0, I_d), \quad
#'       y_i = \arg\max_{y \in \mathcal{Y}} \zeta_i^\top y,}
#' where \eqn{\mathcal{Y} = \{y \in \{0,1\}^d : Ay \le b\}} with totally
#' unimodular `A`, by MH-within-Gibbs sampling on the dual-certificate
#' augmented posterior (`method = "mhwg"`). `method = "unconstrained"` fits
#' independent Albert-Chib probit regressions per coordinate, ignoring the
#' constraints — useful as a baseline demonstrating constraint-ignoring bias.
#'
#' @param Y Response matrix (`n` x `d`), rows in \eqn{\mathcal{Y}}.
#' @param X Covariate matrix (`n` x `p`).
#' @param constraints A [crr_constraints] object. Ignored (with a warning)
#'   for `method = "unconstrained"`.
#' @param method Sampler: `"mhwg"` (the paper's method) or
#'   `"unconstrained"` (baseline).
#' @param kernel Dual kernel for `"mhwg"`: `"exponential"` or
#'   `"half_gaussian"`.
#' @param prior A [crr_prior] object.
#' @param n_iter Total MCMC iterations per chain (including warmup).
#' @param warmup Warmup iterations discarded by `summary()`, `plot()`, and
#'   the `as.mcmc`/`as_draws` converters. All draws are stored.
#' @param thin Thinning applied by the converters (draws are stored unthinned).
#' @param chains Number of independent chains, run serially.
#' @param seed Optional integer seed (`set.seed()` is called if supplied).
#' @param control A [crr_control] object.
#' @param verbose Print progress every 1000 iterations.
#'
#' @return An object of class `crr_fit`; see [summary.crr_fit()].
#'
#' @examples
#' con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
#' sim <- simulate_crr(n = 50, p = 2, constraints = con, seed = 1)
#' fit <- crr(sim$Y, sim$X, con, n_iter = 200, warmup = 100,
#'            chains = 1, seed = 1)
#' summary(fit)
#'
#' @export
crr <- function(Y, X, constraints,
                method = c("mhwg", "unconstrained"),
                kernel = c("exponential", "half_gaussian"),
                prior = crr_prior(),
                n_iter = 5000, warmup = floor(n_iter / 2), thin = 1,
                chains = 1, seed = NULL,
                control = crr_control(),
                verbose = FALSE) {
  method <- match.arg(method)
  kernel <- match.arg(kernel)
  stopifnot(
    inherits(prior, "crr_prior"),
    inherits(control, "crr_control"),
    n_iter >= 2, warmup >= 0, warmup < n_iter, thin >= 1, chains >= 1
  )

  Y <- as.matrix(Y)
  X <- as.matrix(X)
  n <- nrow(Y)
  d <- ncol(Y)
  p <- ncol(X)
  if (nrow(X) != n) stop("nrow(X) must equal nrow(Y)", call. = FALSE)

  if (method == "mhwg") {
    stopifnot(inherits(constraints, "crr_constraints"))
    if (constraints$d != d) {
      stop("constraints$d must equal ncol(Y)", call. = FALSE)
    }
    if (!all(is_feasible(constraints, Y))) {
      stop("All rows of Y must be feasible under the constraints",
           call. = FALSE)
    }
  } else if (!missing(constraints) && !is.null(constraints)) {
    warning("constraints are ignored for method = \"unconstrained\"",
            call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  sampler <- switch(method,
    mhwg = fit_mhwg,
    unconstrained = fit_unconstrained
  )

  draws <- array(
    NA_real_, dim = c(n_iter, chains, p * d),
    dimnames = list(
      NULL, paste0("chain:", seq_len(chains)),
      paste0("beta[", rep(seq_len(p), d), ",", rep(seq_len(d), each = p), "]")
    )
  )
  accept_rate <- numeric(chains)
  timing <- numeric(chains)

  zeta_block_tuned <- integer(chains)
  for (ch in seq_len(chains)) {
    t0 <- proc.time()[["elapsed"]]
    res <- sampler(Y, X, constraints, kernel, prior, n_iter, warmup, control,
                   verbose)
    timing[ch] <- proc.time()[["elapsed"]] - t0
    draws[, ch, ] <- res$draws
    accept_rate[ch] <- res$accept_rate
    zeta_block_tuned[ch] <- res$zeta_block
  }

  structure(
    list(
      draws = draws,
      n = n, p = p, d = d,
      method = method, kernel = if (method == "mhwg") kernel else NA_character_,
      prior = prior, control = control,
      n_iter = n_iter, warmup = warmup, thin = thin, chains = chains,
      accept_rate = accept_rate, timing = timing,
      zeta_block_tuned = zeta_block_tuned,
      constraints = if (method == "mhwg") constraints else NULL,
      X = X, Y = Y
    ),
    class = "crr_fit"
  )
}

# Dual-averaging controller (Hoffman & Gelman 2014, NUTS Alg. 5) on the knob
# x = log(block size), targeting a fixed MH acceptance rate. Acceptance is
# monotone decreasing in the block size, matching the NUTS step-size sign, so
# the update applies unchanged. The caller adapts during warmup then freezes,
# leaving the sampling phase a fixed proposal (valid MCMC).
.zeta_block_controller <- function(d, k0, target,
                                   gamma = 0.05, t0 = 10, kappa = 0.75) {
  mu <- log(k0); x <- log(k0); x_bar <- 0; H_bar <- 0; t <- 0L
  frozen <- FALSE; x_use <- log(k0)
  # Cap x at log(d) before exp: when even k = d clears the target acceptance,
  # x grows without bound and exp(x) would overflow to NA. d is the ceiling.
  to_k <- function(xx) max(1L, min(d, as.integer(round(exp(min(xx, log(d)))))))
  list(
    block = function() to_k(x_use),
    update = function(acc) {
      if (frozen) return(invisible(NULL))
      t <<- t + 1L
      eta <- 1 / (t + t0)
      H_bar <<- (1 - eta) * H_bar + eta * (target - acc)
      x <<- mu - sqrt(t) / gamma * H_bar
      w <- t^(-kappa)
      x_bar <<- w * x + (1 - w) * x_bar
      x_use <<- x
      invisible(NULL)
    },
    # Switch to the averaged knob, unless no update ran (x_bar is still 0).
    freeze = function() { frozen <<- TRUE; if (t > 0L) x_use <<- x_bar
                          invisible(NULL) }
  )
}

# Number of iterations between adaptive block-size updates during warmup.
.zeta_adapt_window <- 25L

fit_mhwg <- function(Y, X, constraints, kernel, prior, n_iter, warmup, control,
                     verbose) {
  n <- nrow(Y); d <- ncol(Y); p <- ncol(X)

  precomp <- coef_precompute(X, prior)
  dual <- init_dual(constraints, Y)
  U <- dual$U
  active <- dual$active

  beta <- matrix(0, p, d)
  Mu <- X %*% beta
  UA <- U %*% constraints$A
  zeta <- draw_utility(matrix(NA_real_, n, d), Mu, Y, UA)

  adaptive <- identical(control$zeta_block, "adaptive")
  ctrl <- if (adaptive)
    .zeta_block_controller(d, min(d, 100L), control$zeta_target_accept) else NULL
  # Shrink the adaptation window on short warmups so the controller still gets
  # several updates instead of silently never adapting.
  window <- if (warmup > 0L)
    max(1L, min(.zeta_adapt_window, warmup %/% 8L)) else .zeta_adapt_window
  win_acc <- numeric(0)

  draws <- matrix(NA_real_, n_iter, p * d)
  # Acceptance is accumulated over the post-warmup sampling phase only, so it
  # reflects the frozen proposal rather than warmup adaptation transients.
  accept_sum <- 0
  n_sample <- n_iter - warmup

  for (iter in seq_len(n_iter)) {
    U <- sample_dual(constraints, zeta, Y, U, active, kernel, control)

    block <- if (adaptive) ctrl$block() else NULL
    step <- sample_utility(constraints, zeta, Y, U, active, Mu, kernel,
                           control, block = block)
    zeta <- step$zeta
    acc <- mean(step$accept)
    if (iter > warmup) accept_sum <- accept_sum + acc

    beta <- update_coef(X, zeta, precomp)
    Mu <- X %*% beta

    draws[iter, ] <- as.vector(beta)

    if (adaptive && iter <= warmup) {
      win_acc <- c(win_acc, acc)
      if (length(win_acc) >= window) {
        ctrl$update(mean(win_acc)); win_acc <- numeric(0)
      }
      if (iter == warmup) ctrl$freeze()
    }
    if (verbose && iter %% 1000 == 0) {
      message("iteration ", iter, "/", n_iter,
              " (acceptance ", format(accept_sum / max(1L, iter - warmup),
                                      digits = 3),
              if (adaptive) paste0(", zeta_block ", ctrl$block()) else "", ")")
    }
  }

  tuned <- if (adaptive) ctrl$block() else min(d, control$zeta_block)
  list(draws = draws, accept_rate = accept_sum / n_sample, zeta_block = tuned)
}

fit_unconstrained <- function(Y, X, constraints, kernel, prior, n_iter, warmup,
                              control, verbose) {
  # Albert-Chib probit per coordinate: y_ij = 1{zeta_ij > 0}.
  n <- nrow(Y); d <- ncol(Y); p <- ncol(X)

  precomp <- coef_precompute(X, prior)
  beta <- matrix(0, p, d)
  Mu <- X %*% beta
  zero_bound <- matrix(0, n, d)
  zeta <- draw_utility(matrix(NA_real_, n, d), Mu, Y, zero_bound)

  draws <- matrix(NA_real_, n_iter, p * d)

  for (iter in seq_len(n_iter)) {
    zeta <- draw_utility(zeta, Mu, Y, zero_bound)
    beta <- update_coef(X, zeta, precomp)
    Mu <- X %*% beta
    draws[iter, ] <- as.vector(beta)
    if (verbose && iter %% 1000 == 0) {
      message("iteration ", iter, "/", n_iter)
    }
  }

  list(draws = draws, accept_rate = 1, zeta_block = NA_integer_)
}
