# Data simulation utilities (require the 'lpSolve' package).

require_lpsolve <- function() {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required for this function", call. = FALSE)
  }
}

ilp_argmax <- function(utility, constraints) {
  require_lpsolve()
  lpSolve::lp(
    direction = "max",
    objective.in = utility,
    const.mat = constraints$A,
    const.dir = rep("<=", constraints$m),
    const.rhs = constraints$b,
    all.bin = TRUE
  )$solution
}

#' Generate a random totally unimodular constraint system
#'
#' For small problems (`d * m <= 50` and `m >= 2` fails), draws random
#' \{-1, 0, 1\} matrices until one is totally unimodular, with `b = 1`. For
#' larger problems, generates a random network (incidence-type) matrix with
#' one `+1` and one `-1` per row — totally unimodular by construction — and
#' Bernoulli right-hand sides. This mirrors the simulation design of the
#' paper.
#'
#' @param d Response dimension.
#' @param m Number of constraints.
#'
#' @return A [crr_constraints] object.
#'
#' @export
random_constraints <- function(d, m) {
  if ((d * m <= 50) || (m < 2)) {
    repeat {
      A <- round(matrix(stats::runif(m * d, -1, 1), m, d))
      if (any(rowSums(abs(A)) == 0)) next
      if (is_tum(A)) break
    }
    b <- rep(1, m)
    crr_constraints(A, b, check_tum = FALSE)
  } else {
    A <- matrix(0, m, d)
    for (i in seq_len(m)) {
      ind <- sample.int(d, 2, replace = FALSE)
      A[i, ind[1]] <- 1
      A[i, ind[2]] <- -1
    }
    b <- sample(0:1, m, replace = TRUE)
    crr_constraints(A, b, check_tum = FALSE)
  }
}

#' Simulate combinatorial response regression data
#'
#' Draws \eqn{X_{ik} \sim N(0,1)}, coefficients
#' \eqn{\beta_{kj} \sim N(0,1)} (unless supplied), latent utilities
#' \eqn{\zeta_i = \beta^\top x_i + \varepsilon_i} with standard normal noise,
#' and responses \eqn{y_i = \arg\max_{y \in \mathcal{Y}} \zeta_i^\top y}
#' solved by integer programming.
#'
#' @param n Number of observations.
#' @param p Number of covariates.
#' @param constraints A [crr_constraints] object defining \eqn{\mathcal{Y}}.
#' @param beta Optional true coefficient matrix (`p` x `d`).
#' @param seed Optional integer seed.
#'
#' @return List with `Y` (`n` x `d`), `X` (`n` x `p`), `beta` (`p` x `d`,
#'   the truth), `zeta` (`n` x `d`), and `constraints`.
#'
#' @examples
#' con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
#' sim <- simulate_crr(n = 20, p = 2, constraints = con, seed = 1)
#'
#' @export
simulate_crr <- function(n, p, constraints, beta = NULL, seed = NULL) {
  stopifnot(inherits(constraints, "crr_constraints"))
  require_lpsolve()
  if (!is.null(seed)) set.seed(seed)

  d <- constraints$d
  X <- matrix(stats::rnorm(n * p), n, p)
  if (is.null(beta)) {
    beta <- matrix(stats::rnorm(p * d), p, d)
  } else {
    beta <- as.matrix(beta)
    stopifnot(nrow(beta) == p, ncol(beta) == d)
  }

  zeta <- X %*% beta + matrix(stats::rnorm(n * d), n, d)
  Y <- t(apply(zeta, 1, ilp_argmax, constraints = constraints))

  list(Y = Y, X = X, beta = beta, zeta = zeta, constraints = constraints)
}
