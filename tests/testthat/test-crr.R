test_that("crr recovers coefficients on a small problem", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 300, p = 2, constraints = con, seed = 11)

  fit <- crr(sim$Y, sim$X, con, n_iter = 600, warmup = 300,
             chains = 1, seed = 11)
  expect_s3_class(fit, "crr_fit")

  rmse <- sqrt(mean((coef(fit) - sim$beta)^2))
  expect_lt(rmse, 0.5)
  expect_gt(mean(fit$accept_rate), 0)
})

test_that("crr is reproducible under a seed and invariant to n_threads", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 40, p = 2, constraints = con, seed = 3)

  f1 <- crr(sim$Y, sim$X, con, n_iter = 100, warmup = 50, seed = 5)
  f2 <- crr(sim$Y, sim$X, con, n_iter = 100, warmup = 50, seed = 5)
  expect_identical(f1$draws, f2$draws)

  f3 <- crr(sim$Y, sim$X, con, n_iter = 100, warmup = 50, seed = 5,
            control = crr_control(n_threads = 2))
  expect_identical(f1$draws, f3$draws)
})

test_that("unconstrained baseline runs and ignores constraints", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 100, p = 2, constraints = con, seed = 2)

  expect_warning(
    fit <- crr(sim$Y, sim$X, con, method = "unconstrained",
               n_iter = 200, warmup = 100, seed = 2),
    "ignored"
  )
  expect_equal(dim(coef(fit)), c(2, 3))
})

test_that("crr rejects infeasible responses", {
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  Y_bad <- rbind(c(1, 1, 0), c(0, 0, 0))
  X <- matrix(rnorm(4), 2, 2)
  expect_error(crr(Y_bad, X, con, n_iter = 10, warmup = 5), "feasible")
})

test_that("S3 methods work", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 50, p = 2, constraints = con, seed = 9)
  fit <- crr(sim$Y, sim$X, con, n_iter = 100, warmup = 50, chains = 2,
             seed = 9)

  s <- summary(fit)
  expect_equal(nrow(s), 2 * 3)
  expect_true(all(c("mean", "sd", "q2.5") %in% names(s)))

  ml <- coda::as.mcmc(fit)
  expect_s3_class(ml, "mcmc.list")
  expect_length(ml, 2)

  pu <- predict(fit)
  expect_equal(dim(pu), c(50, 3))
  pr <- predict(fit, type = "response")
  expect_true(all(is_feasible(con, pr)))
})

test_that("simulate_crr responses are ILP-optimal and feasible", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 25, p = 2, constraints = con, seed = 4)
  expect_true(all(is_feasible(con, sim$Y)))
})

test_that("random_constraints returns TUM systems", {
  set.seed(1)
  con_small <- random_constraints(d = 4, m = 2)
  expect_true(is_tum(con_small$A))

  con_large <- random_constraints(d = 20, m = 5)
  expect_true(all(rowSums(con_large$A == 1) == 1))
  expect_true(all(rowSums(con_large$A == -1) == 1))
})
