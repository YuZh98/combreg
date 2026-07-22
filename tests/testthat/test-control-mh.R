# The MH certificate draw has its own hit-and-run sweep count, separate from
# the dual refresh. The default must leave behaviour unchanged.

make_problem <- function(seed = 1) {
  set.seed(seed)
  con <- crr_constraints(rbind(c(1, 1, 0, 0), c(0, 0, 1, 1),
                               c(1, 0, 1, 0), c(0, 1, 0, 1)), b = rep(1, 4))
  skip_if_not_installed("lpSolve")
  simulate_crr(n = 40, p = 2, constraints = con, seed = seed)
}

test_that("n_iter_hit_and_run_mh defaults to n_iter_hit_and_run", {
  ctl <- crr_control(n_iter_hit_and_run = 17)
  expect_identical(ctl$n_iter_hit_and_run_mh, 17L)
  expect_identical(crr_control()$n_iter_hit_and_run_mh, 50L)
})

test_that("the default leaves crr() output bit-for-bit unchanged", {
  sim <- make_problem()
  a <- crr(sim$Y, sim$X, sim$constraints, n_iter = 60, warmup = 20,
           chains = 1, seed = 3,
           control = crr_control(n_iter_hit_and_run = 8))
  b <- crr(sim$Y, sim$X, sim$constraints, n_iter = 60, warmup = 20,
           chains = 1, seed = 3,
           control = crr_control(n_iter_hit_and_run = 8,
                                 n_iter_hit_and_run_mh = 8))
  expect_identical(a$draws, b$draws)
})

test_that("the new control reaches the MH step and changes the chain", {
  sim <- make_problem()
  base <- crr(sim$Y, sim$X, sim$constraints, n_iter = 60, warmup = 20,
              chains = 1, seed = 3,
              control = crr_control(n_iter_hit_and_run = 8))
  more <- crr(sim$Y, sim$X, sim$constraints, n_iter = 60, warmup = 20,
              chains = 1, seed = 3,
              control = crr_control(n_iter_hit_and_run = 8,
                                    n_iter_hit_and_run_mh = 80))
  expect_false(identical(base$draws, more$draws))
})

test_that("the dual refresh is unaffected by the MH sweep count", {
  sim <- make_problem(2)
  dual <- init_dual(sim$constraints, sim$Y)
  Mu <- sim$X %*% sim$beta
  zeta <- draw_utility(matrix(NA_real_, nrow(sim$Y), ncol(sim$Y)), Mu, sim$Y,
                       dual$U %*% sim$constraints$A)
  draw <- function(mh) {
    set.seed(11)
    sample_dual(sim$constraints, zeta, sim$Y, dual$U, dual$active, "exponential",
                crr_control(n_iter_hit_and_run = 6, n_iter_hit_and_run_mh = mh))
  }
  expect_identical(draw(6), draw(400))
})

test_that("sample_utility tolerates a control without the new field", {
  sim <- make_problem(3)
  dual <- init_dual(sim$constraints, sim$Y)
  Mu <- sim$X %*% sim$beta
  zeta <- draw_utility(matrix(NA_real_, nrow(sim$Y), ncol(sim$Y)), Mu, sim$Y,
                       dual$U %*% sim$constraints$A)
  legacy <- crr_control(n_iter_hit_and_run = 6)
  legacy$n_iter_hit_and_run_mh <- NULL
  expect_no_error(
    sample_utility(sim$constraints, zeta, sim$Y, dual$U, dual$active, Mu,
                   "exponential", legacy, block = 4)
  )
})

test_that("n_iter_hit_and_run_mh is validated", {
  expect_error(crr_control(n_iter_hit_and_run_mh = 1))
  expect_error(crr_control(n_iter_hit_and_run_mh = 0))
})
