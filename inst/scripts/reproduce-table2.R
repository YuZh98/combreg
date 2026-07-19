# Reproduce Table 2 of "Statistical Modeling of Combinatorial Response Data"
# (RMSE for estimating beta with the MH-within-Gibbs sampler) using combreg.
#
# Paper settings per (d, m) cell: n = 1000, p = 5, 50000 iterations with 5000
# warmup and thinning 25, exponential dual kernel, 100 inner hit-and-run steps,
# N(0, 1) coefficient prior. Constraints follow the paper's design (random TUM
# matrices with b = 1 for small d * m, random network incidence matrices with
# Bernoulli b otherwise), as implemented by random_constraints().
#
# Two deliberate choices make this reproduction at least as accurate as the
# paper's fixed-block runs:
#
#   * Adaptive latent-utility block (zeta_block = "adaptive"). The paper uses a
#     fixed block of min(d, 100) coordinates, which for tightly constrained
#     cells (large m) is the entire response vector and mixes very slowly
#     (acceptance can fall below 0.1). The adaptive controller tunes the block
#     to a healthy acceptance, so those cells converge within the budget instead
#     of lagging. On loosely constrained cells it leaves the block at the full
#     size, so it never does worse than the fixed rule. Set ZETA_BLOCK <- 100L
#     for literal fixed-block reproduction.
#
#   * Averaging over COMBREG_NREP realizations. Each Table 2 entry is a single
#     simulated data set, so the published values carry realization noise
#     (especially for small cells). Averaging several independent realizations
#     reports a stable estimate rather than one noisy draw.
#
# RMSE at a reduced iteration budget reflects convergence, not bias: shortening
# the run inflates RMSE for both this code and the paper's. Use the full budget
# for a faithful comparison.
#
# Usage:
#   Rscript reproduce-table2.R [d,m ...]              # e.g. Rscript reproduce-table2.R 2,1 10,5
#   COMBREG_NREP=1 COMBREG_NITER=5000 Rscript reproduce-table2.R 2,1   # quick check
#
# Without arguments, runs the low-dimensional cells (2,1) (5,1) (10,5) (20,10).
# Large cells (d up to 1000) match the paper but take hours to days; pass them
# explicitly. Seeds are fixed per (cell, replication), so results are
# reproducible; they match the paper statistically, not bit-for-bit.

library(combreg)

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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("2,1", "5,1", "10,5", "20,10")
cells <- lapply(strsplit(args, ","), as.integer)

ZETA_BLOCK <- "adaptive"   # or 100L for literal fixed-block reproduction
N_REP  <- as.integer(Sys.getenv("COMBREG_NREP",  "5"))
N_ITER <- as.integer(Sys.getenv("COMBREG_NITER", "50000"))
WARMUP <- as.integer(Sys.getenv("COMBREG_WARMUP", as.character(N_ITER %/% 10L)))

control <- crr_control(n_iter_hit_and_run = 100, zeta_block = ZETA_BLOCK,
                       n_threads = 4)

results <- data.frame()
for (cell in cells) {
  d <- cell[1]; m <- cell[2]; key <- paste0(d, ",", m)
  rmses <- numeric(N_REP); accs <- numeric(N_REP)
  t0 <- proc.time()
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
