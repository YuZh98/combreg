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

test_that("crr_ppc calibrates fit statistics against the model", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit()
  ppc <- crr_ppc(fit, n_rep = 20, seed = 1)

  expect_s3_class(ppc, "crr_ppc")
  expect_equal(ppc$n_rep, 20)
  expect_equal(dim(ppc$replicated), c(20, 1 + 2 * fit$d))
  expect_named(ppc$observed,
               c("exact_match", paste0("accuracy[", 1:3, "]"),
                 paste0("marginal_freq[", 1:3, "]")))
  expect_true(all(ppc$p_value >= 0 & ppc$p_value <= 1))
  # well-specified simulation: headline statistic typical of the model
  expect_gt(ppc$p_value[["exact_match"]], 0.05)
  expect_output(print(ppc), "predictive replicates")
})

test_that("crr_ppc works for the unconstrained baseline", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit(chains = 1, method = "unconstrained")
  ppc <- crr_ppc(fit, n_rep = 10, seed = 2)
  expect_true(all(is.finite(ppc$observed)))
  expect_true(all(ppc$p_value >= 0 & ppc$p_value <= 1))
})

test_that("crr_diagnostics builds a complete report", {
  skip_if_not_installed("lpSolve")
  fit <- small_fit()
  diag <- crr_diagnostics(fit, n_rep = 10)

  expect_s3_class(diag, "crr_diagnostics")
  expect_s3_class(diag$ppc, "crr_ppc")
  expect_null(diag$truth)
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
  expect_output(print(diag), "Posterior predictive checks")
})

test_that("crr_diagnostics reports estimation error against known truth", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 60, p = 2, constraints = con, seed = 7)
  fit <- crr(sim$Y, sim$X, con, n_iter = 200, warmup = 100, seed = 7)

  diag <- crr_diagnostics(fit, beta = sim$beta, n_rep = 0)
  expect_null(diag$ppc)
  expect_equal(diag$truth$rmse,
               sqrt(mean((coef(fit) - sim$beta)^2)))
  expect_equal(diag$truth$error, coef(fit) - sim$beta)
  expect_true(diag$truth$coverage >= 0 && diag$truth$coverage <= 1)
  expect_output(print(diag), "Estimation vs known truth")
  expect_output(print(diag), "RMSE")

  expect_error(crr_diagnostics(fit, beta = matrix(0, 3, 3)))
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
  expect_invisible(plot(fit, type = "beta_diff",
                        beta = matrix(0, fit$p, fit$d)))
  expect_error(plot(fit, type = "beta_diff"), "beta argument")
  expect_error(plot(fit, type = "beta_diff", beta = matrix(0, 3, 3)),
               "matrix")
  expect_error(plot(fit, pars = "beta[9,9]"), "Unknown parameters")
})
