# Distributional correctness: on an enumerable problem with a constant
# covariate, the likelihood is multinomial over Y = {(0,0), (1,0), (0,1)}
# with category probabilities computable by 1-D quadrature, so the exact
# posterior is available on a grid. The sampler must reproduce it.

exact_posterior <- function(counts, grid = seq(-1.5, 1.5, by = 0.02)) {
  z <- seq(0, 8, by = 0.02)
  phi <- outer(z, grid, function(z, b) dnorm(z - b))
  Phi <- outer(z, grid, function(z, b) pnorm(z - b))
  P10 <- crossprod(phi, Phi) * 0.02  # P10[i, j] = P(y = (1,0) | b1_i, b2_j)
  P00 <- outer(pnorm(-grid), pnorm(-grid))

  lp <- counts[1] * log(P00) + counts[2] * log(P10) +
    counts[3] * log(t(P10)) +
    outer(dnorm(grid, log = TRUE), dnorm(grid, log = TRUE), `+`)
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  list(
    mean = c(sum(rowSums(w) * grid), sum(colSums(w) * grid)),
    sd = c(sqrt(sum(rowSums(w) * grid^2) - sum(rowSums(w) * grid)^2),
           sqrt(sum(colSums(w) * grid^2) - sum(colSums(w) * grid)^2))
  )
}

test_that("sampler reproduces the exact posterior on an enumerable problem", {
  set.seed(123)
  n <- 400
  beta_true <- c(0.5, -0.3)
  zeta <- cbind(beta_true[1] + rnorm(n), beta_true[2] + rnorm(n))
  Y <- t(apply(zeta, 1, function(z) {
    switch(which.max(c(0, z[1], z[2])), c(0, 0), c(1, 0), c(0, 1))
  }))
  counts <- c(sum(Y[, 1] == 0 & Y[, 2] == 0), sum(Y[, 1] == 1),
              sum(Y[, 2] == 1))

  exact <- exact_posterior(counts)

  con <- crr_constraints(matrix(c(1, 1), 1, 2), b = 1)
  fit <- crr(Y, matrix(1, n, 1), con,
             n_iter = 12000, warmup = 2000, thin = 5, seed = 1)
  s <- summary(fit)

  expect_lt(max(abs(s$mean - exact$mean)), 0.03)
  expect_lt(max(abs(s$sd - exact$sd)), 0.015)
})
