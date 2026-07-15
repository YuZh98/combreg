make_tight_con <- function(d, m, seed = 1) {
  set.seed(seed)
  A <- matrix(0, m, d)
  for (i in seq_len(m)) {
    idx <- sample.int(d, 2)
    A[i, idx[1]] <- 1; A[i, idx[2]] <- -1
  }
  crr_constraints(A, sample(0:1, m, replace = TRUE), check_tum = FALSE)
}

test_that("crr_control validates zeta_block and zeta_target_accept", {
  expect_silent(crr_control(zeta_block = "adaptive"))
  expect_silent(crr_control(zeta_block = 10))
  expect_error(crr_control(zeta_block = "foo"), "positive integer")
  expect_error(crr_control(zeta_block = 0), "positive integer")
  expect_error(crr_control(zeta_block = c(1, 2)), "positive integer")
  expect_error(crr_control(zeta_target_accept = 0))
  expect_error(crr_control(zeta_target_accept = 1))

  ctl <- crr_control()
  expect_identical(ctl$zeta_block, "adaptive")   # new default
  expect_equal(ctl$zeta_target_accept, 0.6)
  expect_identical(crr_control(zeta_block = 7)$zeta_block, 7L)
})

test_that("adaptive zeta_block tunes acceptance toward the target and is stored", {
  skip_if_not_installed("lpSolve")
  con <- make_tight_con(d = 15, m = 10, seed = 4)
  sim <- simulate_crr(n = 120, p = 3, constraints = con, seed = 4)

  fit <- crr(sim$Y, sim$X, con, n_iter = 700, warmup = 400, seed = 7)
  tuned <- fit$zeta_block_tuned
  expect_length(tuned, 1L)
  expect_true(tuned >= 1L && tuned <= 15L)
  # acceptance should land in a band around the 0.6 target (not collapsed).
  expect_gt(mean(fit$accept_rate), 0.4)
  expect_lt(mean(fit$accept_rate), 0.85)
})

test_that("adaptive block never accepts worse than the full-block rule", {
  skip_if_not_installed("lpSolve")
  con <- make_tight_con(d = 15, m = 10, seed = 4)
  sim <- simulate_crr(n = 120, p = 3, constraints = con, seed = 4)

  f_adapt <- crr(sim$Y, sim$X, con, n_iter = 700, warmup = 400, seed = 7)
  f_all <- crr(sim$Y, sim$X, con, n_iter = 700, warmup = 400, seed = 7,
               control = crr_control(zeta_block = 15))
  # Shrinking the block can only raise acceptance; adaptive should not be worse.
  expect_gte(mean(f_adapt$accept_rate), mean(f_all$accept_rate) - 0.05)
  expect_lte(f_adapt$zeta_block_tuned, 15L)
})

test_that("adaptive controller still tunes under a short warmup", {
  skip_if_not_installed("lpSolve")
  con <- make_tight_con(d = 15, m = 10, seed = 4)
  sim <- simulate_crr(n = 120, p = 3, constraints = con, seed = 4)

  # warmup (16) is shorter than the default 25-iteration window; the controller
  # must still adapt rather than silently staying at the untuned min(d, 100).
  f_adapt <- crr(sim$Y, sim$X, con, n_iter = 500, warmup = 16, seed = 7)
  f_all <- crr(sim$Y, sim$X, con, n_iter = 500, warmup = 16, seed = 7,
               control = crr_control(zeta_block = 15))
  # On a tight problem the full block over-rejects; a tuned block accepts more.
  expect_gt(mean(f_adapt$accept_rate), mean(f_all$accept_rate))
  expect_lt(f_adapt$zeta_block_tuned, 15L)
})

test_that("fixed integer zeta_block is respected and reported as tuned", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 80, p = 2, constraints = con, seed = 2)

  fit <- crr(sim$Y, sim$X, con, n_iter = 300, warmup = 150, seed = 2,
             control = crr_control(zeta_block = 2))
  expect_identical(fit$control$zeta_block, 2L)
  expect_identical(fit$zeta_block_tuned, 2L)   # min(d = 3, 2)
})

test_that("adaptive sampling is reproducible under a seed", {
  skip_if_not_installed("lpSolve")
  con <- make_tight_con(d = 10, m = 6, seed = 8)
  sim <- simulate_crr(n = 60, p = 2, constraints = con, seed = 8)

  f1 <- crr(sim$Y, sim$X, con, n_iter = 400, warmup = 200, seed = 9)
  f2 <- crr(sim$Y, sim$X, con, n_iter = 400, warmup = 200, seed = 9)
  expect_identical(f1$draws, f2$draws)
  expect_identical(f1$zeta_block_tuned, f2$zeta_block_tuned)
})

test_that("unconstrained baseline reports NA tuned block", {
  skip_if_not_installed("lpSolve")
  con <- crr_constraints(rbind(c(1, 1, 0), c(0, 1, 1)), b = c(1, 1))
  sim <- simulate_crr(n = 60, p = 2, constraints = con, seed = 2)
  suppressWarnings(
    fit <- crr(sim$Y, sim$X, con, method = "unconstrained",
               n_iter = 150, warmup = 75, seed = 2)
  )
  expect_true(is.na(fit$zeta_block_tuned))
})
