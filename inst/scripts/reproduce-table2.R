# Reproduce Table 2 of "Statistical Modeling of Combinatorial Response Data"
# (Zheng, Ghosh & Duan, 2026+): coefficient RMSE across dimensionalities (d, m).
#
#
# NOTES
#   * Paper design per (d, m) cell: n = 1000, p = 5, exponential dual kernel,
#     100 inner hit-and-run steps, N(0, 1) coefficient prior. Constraints follow
#     random_constraints() (random TUM matrices with b = 1 for small d * m,
#     random network-incidence matrices otherwise).
#   * ZETA_BLOCK = "adaptive" (the package default) tunes the latent-utility
#     block to a healthy acceptance during warmup. The paper's fixed
#     min(d, 100) block mixes very slowly on tightly constrained (large-m) cells,
#     so the adaptive controller reaches the paper's RMSE within the budget
#     instead of lagging, while never doing worse on loose cells. Set
#     ZETA_BLOCK = 100L for literal fixed-block reproduction.

library(combreg)

## ======================= CONFIGURATION  ===============================

CELLS <- list(c(2, 1), c(5, 1), c(10, 5), c(20, 10))  # (d, m) cells to run
N_REP      <- 5           # number of repetitions per (d, m) cell
N_ITER     <- 10000       # MCMC iterations per fit 
WARMUP     <- 5000        # warmup iterations
ZETA_BLOCK <- "adaptive"  # "adaptive" (recommended) or an integer, e.g. 100L
N_THREADS  <- 4           # OpenMP threads for the dual updates

## =======================================================================

paper_rmse <- c("2,1" = 0.046, "5,1" = 0.066, "10,1" = 0.076, "20,1" = 0.065,
                "50,1" = 0.086, "100,1" = 0.079, "200,1" = 0.090,
                "500,1" = 0.094, "1000,1" = 0.091,
                "10,5" = 0.066, "20,5" = 0.084, "50,5" = 0.084,
                "100,5" = 0.081, "200,5" = 0.082, "500,5" = 0.088,
                "1000,5" = 0.087,
                "20,10" = 0.088, "50,10" = 0.099, "100,10" = 0.084,
                "200,10" = 0.095, "500,10" = 0.087, "1000,10" = 0.088,
                "50,20" = 0.116, "100,20" = 0.104, "200,20" = 0.090,
                "500,20" = 0.091, "1000,20" = 0.091,
                "100,50" = 0.209, "200,50" = 0.171, "500,50" = 0.125,
                "1000,50" = 0.104,
                "500,100" = 0.263, "1000,100" = 0.128)


## =======================================================================
# Optional: when called from a terminal, cells may be passed as arguments,
# e.g. `Rscript reproduce-table2.R 2,1 5,1 10,5`. Ignored under Source / Run All.
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) > 0) CELLS <- lapply(strsplit(.args, ","), as.integer)
## =======================================================================

control <- crr_control(n_iter_hit_and_run = 100, zeta_block = ZETA_BLOCK,
                       n_threads = N_THREADS)

results <- data.frame()
for (cell in CELLS) {
  d <- cell[1]; m <- cell[2]; key <- paste0(d, ",", m)
  rmses <- numeric(N_REP); accs <- numeric(N_REP)
  t0 <- proc.time()
  
  # Run N_REP independent realizations of the (d, m) cell
  for (r in seq_len(N_REP)) {
    set.seed(1000L * d + m + r)
    con <- random_constraints(d, m)
    sim <- simulate_crr(n = 1000, p = 5, constraints = con)
    fit <- crr(sim$Y, sim$X, con,
               kernel = "exponential", prior = crr_prior(sd = 1),
               n_iter = N_ITER, warmup = WARMUP, thin = 25,
               control = control, seed = r)
    rmses[r] <- sqrt(mean((coef(fit) - sim$beta)^2))
    accs[r]  <- mean(fit$accept_rate)
  }
  elapsed <- (proc.time() - t0)[["elapsed"]]
  results <- rbind(results, data.frame(
    d = d, m = m, rmse = round(mean(rmses), 3),
    paper = unname(paper_rmse[key]), rmse_sd = round(stats::sd(rmses), 3),
    accept = round(mean(accs), 2), minutes = round(elapsed / 60, 1)))
  cat(sprintf("[done] d=%d m=%d rmse=%.3f (sd %.3f, %d rep) (paper %.3f) accept=%.2f %.1f min\n",
              d, m, mean(rmses), stats::sd(rmses), N_REP, paper_rmse[key],
              mean(accs), elapsed / 60))
}

cat("\n== Table 2 reproduction (mean over ", N_REP, " realization(s), ",
    N_ITER, " iterations, zeta_block = ", as.character(ZETA_BLOCK), ") ==\n",
    sep = "")
print(results, row.names = FALSE)
