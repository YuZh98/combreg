test_that("crr_benchmark compares methods on common data", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 80, p = 2, constraints = con, seed = 21)

  bm <- crr_benchmark(sim$Y, sim$X, con, beta = sim$beta,
                      n_iter = 300, warmup = 150, seed = 21)

  expect_s3_class(bm, "crr_benchmark")
  expect_equal(bm$table$method, c("mhwg", "unconstrained"))
  expect_true(all(c("time_sec", "min_ess", "min_ess_per_sec", "rmse",
                    "coverage", "exact_match") %in% names(bm$table)))
  expect_true(all(bm$table$time_sec >= 0))
  expect_true(all(bm$table$coverage >= 0 & bm$table$coverage <= 1))
  expect_s3_class(bm$fits$mhwg, "crr_fit")

  expect_output(print(bm), "crr_benchmark")
})

test_that("crr_benchmark omits truth-based columns without beta", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 40, p = 2, constraints = con, seed = 22)

  bm <- crr_benchmark(sim$Y, sim$X, con, methods = "mhwg",
                      n_iter = 100, warmup = 50, seed = 22)
  expect_false("rmse" %in% names(bm$table))
  expect_equal(nrow(bm$table), 1)
})
