test_that("formula path matches equivalent matrix path bit-for-bit", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 150, p = 2, constraints = con, seed = 7)

  dat <- data.frame(x1 = sim$X[, 1], x2 = sim$X[, 2])
  X_mat <- model.matrix(~ x1 + x2, dat)

  f_mat <- crr(sim$Y, X_mat, con, n_iter = 300, warmup = 150,
               chains = 1, seed = 42)
  f_form <- crr(sim$Y, ~ x1 + x2, con, data = dat, n_iter = 300,
                warmup = 150, chains = 1, seed = 42)

  expect_identical(f_form$draws, f_mat$draws)
  expect_identical(coef(f_form), coef(f_mat))
})

test_that("a factor predictor expands to the right number of columns", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 150, p = 2, constraints = con, seed = 8)

  dat <- data.frame(
    x1 = sim$X[, 1],
    g = factor(rep(c("a", "b", "c"), length.out = 150))
  )
  fit <- crr(sim$Y, ~ x1 + g, con, data = dat, n_iter = 300,
             warmup = 150, chains = 1, seed = 8)

  # intercept + x1 + two dummy columns for the 3-level factor
  expect_equal(fit$p, 4L)
  expect_equal(nrow(coef(fit)), 4L)
})

test_that("formula with data = NULL errors clearly", {
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  Y <- rbind(c(1, 0, 0), c(0, 1, 0))
  expect_error(
    crr(Y, ~ x1, con, data = NULL, n_iter = 10, warmup = 5),
    "data"
  )
})
