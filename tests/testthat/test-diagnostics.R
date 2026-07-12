small_fit <- function(chains = 2, method = "mhwg", n = 60, n_iter = 200) {
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = n, p = 2, constraints = con, seed = 7)
  if (method == "mhwg") {
    crr(sim$Y, sim$X, con, n_iter = n_iter, warmup = n_iter / 2,
        chains = chains, seed = 7)
  } else {
    crr(sim$Y, sim$X, method = method, n_iter = n_iter,
        warmup = n_iter / 2, chains = chains, seed = 7)
  }
}

test_that("crr_ess and crr_rhat return one finite value per parameter", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit()

  ess <- crr_ess(fit)
  expect_length(ess, fit$p * fit$d)
  expect_true(all(ess > 0))

  rhat <- crr_rhat(fit)
  expect_length(rhat, fit$p * fit$d)
  expect_true(all(is.finite(rhat)))
  expect_true(all(rhat >= 1 - 1e-8))
})

test_that("crr_diagnostics builds a complete report", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit()
  diag <- crr_diagnostics(fit)

  expect_s3_class(diag, "crr_diagnostics")
  expect_equal(nrow(diag$table), fit$p * fit$d)
  expect_true(all(c("mean", "sd", "lower", "upper", "ess", "ess_per_sec",
                    "rhat") %in% names(diag$table)))
  expect_true(all(diag$table$lower <= diag$table$upper))
  expect_true(!is.null(diag$fit_stats))
  expect_true(diag$fit_stats$exact_match >= 0 &&
                diag$fit_stats$exact_match <= 1)
  expect_length(diag$fit_stats$coordinate_accuracy, fit$d)

  expect_output(print(diag), "Model")
  expect_output(print(diag), "Coefficients")
  expect_output(print(diag), "Regression fit")
})

test_that("fitted and residuals are consistent with predict", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit(chains = 1)

  Y_hat <- fitted(fit)
  expect_equal(Y_hat, predict(fit, type = "response"))
  res <- residuals(fit)
  expect_equal(dim(res), dim(fit$Y))
  expect_true(all(res %in% c(-1, 0, 1)))
})

test_that("predict and fitted work for the unconstrained baseline", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit(chains = 1, method = "unconstrained")

  Y_hat <- predict(fit, type = "response")
  expect_true(all(Y_hat %in% c(0, 1)))
  expect_equal(Y_hat, (predict(fit, type = "utility") > 0) * 1)
  expect_equal(residuals(fit), fit$Y - Y_hat)
})

test_that("all plot types run without error", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit()

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  })
  for (tp in c("trace", "acf", "violin", "ess", "ess_time", "residual")) {
    expect_invisible(plot(fit, type = tp))
  }
  expect_error(plot(fit, pars = "beta[9,9]"), "Unknown parameters")
})
