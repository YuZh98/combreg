make_setup <- function(n = 30, seed = 42) {
  set.seed(seed)
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = n, p = 2, constraints = con)
  list(con = con, sim = sim)
}

test_that("init_dual respects active-constraint support", {
  s <- make_setup()
  dual <- init_dual(s$con, s$sim$Y)
  expect_equal(dim(dual$U), c(30, 2))
  expect_true(all(dual$U[dual$active == 0] == 0))
  expect_true(all(dual$U >= 0))
})

test_that("sample_dual keeps dual certificates feasible", {
  skip_if_not_installed("lpSolve")
  s <- make_setup()
  Y <- s$sim$Y
  dual <- init_dual(s$con, Y)
  Mu <- matrix(0, nrow(Y), ncol(Y))
  UA <- dual$U %*% s$con$A
  zeta <- draw_utility(matrix(NA_real_, nrow(Y), ncol(Y)), Mu, Y, UA)

  U_new <- sample_dual(s$con, zeta, Y, dual$U, dual$active)
  feas <- dual_feasible(s$con, zeta, Y, U_new)
  expect_true(all(feas == 1))
  expect_true(all(U_new >= -1e-9))
  expect_true(all(U_new[dual$active == 0] == 0))
})

test_that("draw_utility respects truncation bounds", {
  s <- make_setup()
  Y <- s$sim$Y
  dual <- init_dual(s$con, Y)
  UA <- dual$U %*% s$con$A
  Mu <- matrix(0, nrow(Y), ncol(Y))
  zeta <- draw_utility(matrix(NA_real_, nrow(Y), ncol(Y)), Mu, Y, UA)

  expect_true(all(zeta[Y == 1] >= UA[Y == 1]))
  expect_true(all(zeta[Y == 0] <= UA[Y == 0]))
})

test_that("update_coef targets the conjugate posterior", {
  set.seed(7)
  n <- 2000; p <- 2; d <- 1
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- matrix(c(1, -2), p, d)
  zeta <- X %*% beta_true + matrix(rnorm(n * d), n, d)

  pc <- coef_precompute(X, crr_prior(sd = 10))
  draws <- replicate(500, update_coef(X, zeta, pc))
  post_mean <- apply(draws, c(1, 2), mean)
  expect_equal(as.vector(post_mean), as.vector(beta_true),
               tolerance = 0.15)
})
