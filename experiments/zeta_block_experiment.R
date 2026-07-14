# =============================================================================
# Experiment: acceptance rate vs. number of zeta coordinates updated per MH step
# =============================================================================
#
# The MH-within-Gibbs sampler updates a random subset of the d response
# coordinates of zeta per iteration and accepts/rejects each observation's
# proposal as a unit. Acceptance collapses as the subset grows, so the block
# size trades acceptance against how many coordinates actually move. This
# experiment sweeps five block-size strategies across a (d, m) grid to find
# the strategy that maximizes mixing efficiency without starving acceptance.
#
# The sampler step reuses the package's exported Layer-1 primitives, so the
# chain here is identical to crr(method = "mhwg") except that the block-size
# rule is pluggable (verified in the smoke test against crr()).
#
# HOW TO RUN
#   Smoke (correctness, minutes):     SMOKE <- TRUE  (default)
#   Full study (hours, run yourself): set SMOKE <- FALSE below, then
#     Rscript experiments/zeta_block_experiment.R
#   Outputs: experiments/output/zeta_block_results.csv
#            experiments/output/zeta_block_report.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(combreg)
  library(coda)
})

# ------------------------------------------------------------------ CONFIG ---
SMOKE <- TRUE   # <- set FALSE for the full study

CONFIG <- if (SMOKE) {
  list(
    d_values = c(2, 10, 40),
    m_pool   = c(1, 5),
    n = 60, p = 3,
    n_iter = 400, warmup = 150, n_iter_har = 15,
    cores = 1, seed_base = 100L,
    out_dir = file.path("experiments", "output")
  )
} else {
  list(
    d_values = c(2, 5, 10, 20, 50, 100, 200),
    m_pool   = c(1, 5, 10, 20, 50),
    n = 1000, p = 5,
    n_iter = 5000, warmup = 1000, n_iter_har = 100,
    cores = max(1L, parallel::detectCores() - 1L), seed_base = 100L,
    out_dir = file.path("experiments", "output")
  )
}

# ------------------------------------------------------------- STRATEGIES ---
# Each rule maps d -> number of coordinates to update. Selection is uniform
# random for all five, so block SIZE is the sole treatment variable.
STRATEGIES <- list(
  "Default (min(d,100))" = function(d) min(d, 100),
  "All (d)"              = function(d) d,
  "Fixed 10"             = function(d) min(d, 10),
  "Proportional (0.25d)" = function(d) max(1, round(0.25 * d)),
  "Sqrt (2*sqrt(d))"     = function(d) max(1, round(2 * sqrt(d)))
)
STRAT_NAMES <- names(STRATEGIES)
STRAT_COL <- stats::setNames(
  c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"), STRAT_NAMES)
STRAT_PCH <- stats::setNames(c(16, 17, 15, 18, 8), STRAT_NAMES)

# ------------------------------------------------------------ CONSTRAINTS ---
# Directed-graph incidence design (one +1 and one -1 per constraint row):
# totally unimodular by construction and matches the paper's high-dimensional
# simulation design, avoiding the slow random-TUM search for moderate sizes.
make_constraints <- function(d, m) {
  A <- matrix(0, m, d)
  for (i in seq_len(m)) {
    idx <- sample.int(d, 2)
    A[i, idx[1]] <- 1
    A[i, idx[2]] <- -1
  }
  b <- if (m >= 2) sample(0:1, m, replace = TRUE) else 1
  crr_constraints(A, b, check_tum = FALSE)
}

m_values_for <- function(d) {
  ms <- CONFIG$m_pool[CONFIG$m_pool < d]
  if (length(ms) == 0) 1L else ms
}

# --------------------------------------------------------------- SAMPLER ---
# MH-within-Gibbs chain with a pluggable block-size rule. Mirrors
# combreg:::fit_mhwg / sample_utility exactly aside from block_fn.
run_chain <- function(Y, X, con, block_fn, kernel, cfg, beta_true = NULL) {
  n <- nrow(Y); d <- ncol(Y); p <- ncol(X)
  control <- crr_control(n_iter_hit_and_run = cfg$n_iter_har, n_threads = 1)

  precomp <- coef_precompute(X, crr_prior(sd = 1))
  dual <- init_dual(con, Y); U <- dual$U; active <- dual$active
  beta <- matrix(0, p, d); Mu <- X %*% beta
  UA <- U %*% con$A
  zeta <- draw_utility(matrix(NA_real_, n, d), Mu, Y, UA)

  draws <- matrix(NA_real_, cfg$n_iter, p * d)
  acc_trace <- numeric(cfg$n_iter)
  k_trace <- integer(cfg$n_iter)

  t0 <- proc.time()[["elapsed"]]
  for (iter in seq_len(cfg$n_iter)) {
    U <- sample_dual(con, zeta, Y, U, active, kernel, control)

    UA <- U %*% con$A
    k <- max(1L, min(d, as.integer(block_fn(d))))
    subset <- sample.int(d, k, replace = FALSE)
    zeta_new <- draw_utility(zeta, Mu, Y, UA, subset)
    zeta_tilde <- (Y > 0.5) * pmax(zeta, zeta_new) +
                  (Y < 0.5) * pmin(zeta, zeta_new)
    U_star <- sample_dual(con, zeta_tilde, Y, U, active, kernel, control)
    acc <- dual_feasible(con, zeta, Y, U_star, control)
    zeta <- zeta_new * acc + zeta * (1 - acc)

    beta <- update_coef(X, zeta, precomp); Mu <- X %*% beta
    draws[iter, ] <- as.vector(beta)
    acc_trace[iter] <- mean(acc)
    k_trace[iter] <- k
  }
  elapsed <- proc.time()[["elapsed"]] - t0

  keep <- seq.int(cfg$warmup + 1L, cfg$n_iter)
  post <- draws[keep, , drop = FALSE]
  ess <- tryCatch(coda::effectiveSize(coda::mcmc(post)),
                  error = function(e) rep(NA_real_, ncol(post)))
  beta_hat <- matrix(colMeans(post), p, d)

  list(
    acc = mean(acc_trace[keep]),
    k = mean(k_trace[keep]),
    ess_med = stats::median(ess, na.rm = TRUE),
    ess_min = min(ess, na.rm = TRUE),
    time = elapsed,
    rmse = if (is.null(beta_true)) NA_real_ else
      sqrt(mean((beta_hat - beta_true)^2))
  )
}

# ------------------------------------------------------------------- GRID ---
build_cells <- function() {
  cells <- list()
  for (d in CONFIG$d_values) {
    for (m in m_values_for(d)) {
      cells[[length(cells) + 1L]] <- list(d = d, m = m)
    }
  }
  cells
}

run_cell <- function(cell) {
  d <- cell$d; m <- cell$m
  cell_seed <- CONFIG$seed_base + d * 1000L + m
  set.seed(cell_seed)
  con <- make_constraints(d, m)
  sim <- simulate_crr(CONFIG$n, CONFIG$p, con)

  rows <- lapply(STRAT_NAMES, function(sname) {
    set.seed(cell_seed + 7L)  # identical chain start across strategies
    res <- run_chain(sim$Y, sim$X, con, STRATEGIES[[sname]],
                     kernel = "exponential", cfg = CONFIG,
                     beta_true = sim$beta)
    data.frame(
      strategy = sname, d = d, m = m,
      k = res$k, acc = res$acc,
      throughput = res$k * res$acc,        # coords moved per proposal per obs
      ess_med = res$ess_med, ess_min = res$ess_min,
      ess_per_sec = res$ess_med / res$time,
      time = res$time, rmse = res$rmse,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# --------------------------------------------------------------- SANITY -----
sanity_checks <- function(df) {
  msg <- function(ok, label) {
    cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
    if (!ok && SMOKE) stop("Sanity check failed: ", label, call. = FALSE)
  }
  cat("Sanity checks:\n")

  # Strategy size rules at a large d.
  msg(STRATEGIES[["Default (min(d,100))"]](200) == 100 &&
      STRATEGIES[["All (d)"]](200) == 200 &&
      STRATEGIES[["Fixed 10"]](200) == 10 &&
      STRATEGIES[["Proportional (0.25d)"]](200) == 50 &&
      STRATEGIES[["Sqrt (2*sqrt(d))"]](200) == 28,
      "block-size rules return expected k at d=200")

  msg(all(df$acc >= 0 & df$acc <= 1), "acceptance rates in [0,1]")
  msg(all(df$k >= 1 & df$k <= df$d), "block sizes in [1,d]")
  msg(all(df$ess_med > 0, na.rm = TRUE), "effective sample sizes positive")

  # Core monotonic relationship: acceptance falls as more coords are updated.
  sp <- suppressWarnings(stats::cor(df$k, df$acc, method = "spearman"))
  msg(!is.na(sp) && sp < 0,
      sprintf("acceptance decreases with block size (Spearman rho = %.2f)", sp))

  # Throughput should peak at an interior block size somewhere, i.e. neither
  # the smallest nor the largest strategy always wins.
  best <- do.call(rbind, by(df, list(df$d, df$m), function(g) {
    g[which.max(g$throughput), c("d", "m", "strategy")]
  }))
  interior <- any(!best$strategy %in% c("Fixed 10", "All (d)"))
  msg(interior, "an interior block size maximizes throughput in some cell")
}

verify_against_crr <- function() {
  # The pluggable Default chain must reproduce crr()'s acceptance (same
  # sampler, same block rule). Statistical agreement over a short chain.
  set.seed(4242)
  con <- make_constraints(10, 5)
  sim <- simulate_crr(80, 3, con)
  cfg <- modifyList(CONFIG, list(n = 80, p = 3, n_iter = 400, warmup = 100,
                                 n_iter_har = 20))

  set.seed(11)
  ours <- run_chain(sim$Y, sim$X, con, function(d) min(d, 100),
                    "exponential", cfg)
  set.seed(11)
  fit <- crr(sim$Y, sim$X, con, n_iter = 400, warmup = 100, seed = 11,
             control = crr_control(n_iter_hit_and_run = 20, n_threads = 1))
  diff <- abs(ours$acc - mean(fit$accept_rate))
  cat(sprintf("Cross-check vs crr(): harness acc = %.3f, crr acc = %.3f, |diff| = %.3f\n",
              ours$acc, mean(fit$accept_rate), diff))
  if (SMOKE && diff > 0.15)
    stop("Harness acceptance diverges from crr(); logic mismatch.",
         call. = FALSE)
  cat(sprintf("  [%s] harness reproduces crr() acceptance\n",
              if (diff <= 0.15) "PASS" else "FAIL"))
}

# ------------------------------------------------------------------ PLOTS ---
panel_by_m <- function(df, yvar, ylab, logy = FALSE) {
  ms <- sort(unique(df$m))
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  nr <- ceiling(length(ms) / 2)
  graphics::par(mfrow = c(nr, min(2, length(ms))), mar = c(4, 4, 2.5, 1),
                oma = c(0, 0, 2, 0))
  for (mm in ms) {
    sub <- df[df$m == mm, ]
    yl <- range(sub[[yvar]][is.finite(sub[[yvar]])], na.rm = TRUE)
    if (logy) yl[1] <- max(yl[1], min(sub[[yvar]][sub[[yvar]] > 0], na.rm = TRUE))
    plot(NA, xlim = range(df$d), ylim = yl, log = if (logy) "xy" else "x",
         xlab = "d (response dimension)", ylab = ylab,
         main = sprintf("m = %d", mm))
    for (sname in STRAT_NAMES) {
      g <- sub[sub$strategy == sname, ]
      g <- g[order(g$d), ]
      if (nrow(g)) {
        lines(g$d, g[[yvar]], col = STRAT_COL[sname], lwd = 2)
        points(g$d, g[[yvar]], col = STRAT_COL[sname], pch = STRAT_PCH[sname])
      }
    }
    if (mm == ms[1]) {
      legend("topright", legend = STRAT_NAMES, col = STRAT_COL,
             pch = STRAT_PCH, lwd = 2, cex = 0.7, bg = "white")
    }
  }
  graphics::mtext(ylab, outer = TRUE, cex = 1.1, font = 2)
}

plot_acc_vs_k <- function(df) {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mar = c(4, 4, 3, 1))
  ds <- sort(unique(df$d))
  pal <- grDevices::hcl.colors(length(ds), "viridis")
  col <- pal[match(df$d, ds)]
  plot(df$k, df$acc, col = col, pch = 19, log = "x",
       xlab = "number of coordinates updated (k, log scale)",
       ylab = "acceptance rate",
       main = "Acceptance vs. block size (color = d)")
  o <- order(df$k)
  lines(lowess(df$k[o], df$acc[o]), col = "grey30", lwd = 2, lty = 2)
  legend("bottomleft", legend = paste0("d=", ds), col = pal, pch = 19,
         cex = 0.7, bg = "white", title = "d")
}

text_page <- function(lines, title = NULL) {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mar = c(1, 1, 1, 1))
  plot.new()
  if (!is.null(title)) graphics::mtext(title, side = 3, line = -2, font = 2,
                                       cex = 1.3)
  y <- 0.9
  for (ln in lines) {
    graphics::text(0.02, y, ln, adj = 0, cex = 0.8, family = "mono")
    y <- y - 0.045
  }
}

table_page <- function(tbl, title) {
  gridExtra::grid.arrange(
    gridExtra::tableGrob(tbl, rows = NULL,
                         theme = gridExtra::ttheme_default(base_size = 8)),
    top = title)
}

make_report <- function(df, pdf_path) {
  # Headline pivots.
  acc_pivot <- tapply(df$acc, list(df$strategy, df$d),
                      function(v) round(mean(v), 3))
  acc_pivot <- as.data.frame.matrix(acc_pivot)
  acc_pivot <- cbind(strategy = rownames(acc_pivot), acc_pivot)

  best_eff <- do.call(rbind, by(df, list(df$d, df$m), function(g) {
    g[which.max(g$ess_per_sec), ]
  }))
  best_tbl <- data.frame(
    d = best_eff$d, m = best_eff$m,
    best_strategy = best_eff$strategy,
    k = round(best_eff$k, 1), acc = round(best_eff$acc, 3),
    ess_per_sec = round(best_eff$ess_per_sec, 2))
  best_tbl <- best_tbl[order(best_tbl$d, best_tbl$m), ]

  grDevices::pdf(pdf_path, width = 9, height = 6.5)
  on.exit(grDevices::dev.off())

  text_page(title = "Zeta block-size experiment", c(
    "",
    "Goal: relate MH acceptance to the number of zeta coordinates (k)",
    "updated per proposal, and pick a block-size rule that keeps",
    "acceptance healthy while maximizing mixing efficiency.",
    "",
    "Five strategies (uniform-random coordinate selection; k = f(d)):",
    "  Default (min(d,100))  - current package behavior",
    "  All (d)               - update every coordinate",
    "  Fixed 10              - small constant block",
    "  Proportional (0.25d)  - block scales with d",
    "  Sqrt (2*sqrt(d))      - sublinear scaling",
    "",
    "Metrics per (strategy, d, m):",
    "  acc          acceptance rate (per-observation, post-warmup)",
    "  k            coordinates updated per proposal",
    "  throughput   k * acc  (coords actually moved per proposal)",
    "  ess_per_sec  median beta ESS / wall-clock  (mixing efficiency)",
    "  rmse         ||beta_hat - beta_true||  (accuracy check)",
    "",
    sprintf("Grid: d in {%s}, m in {%s}, n=%d, p=%d, iters=%d (warmup %d).",
            paste(CONFIG$d_values, collapse = ","),
            paste(CONFIG$m_pool, collapse = ","),
            CONFIG$n, CONFIG$p, CONFIG$n_iter, CONFIG$warmup),
    if (SMOKE) "MODE: SMOKE (small; correctness only)." else "MODE: FULL."
  ))

  panel_by_m(df, "acc", "Acceptance rate")
  plot_acc_vs_k(df)
  panel_by_m(df, "throughput", "Throughput  (k * acc)")
  panel_by_m(df, "ess_per_sec", "Mixing efficiency  (median ESS / sec)", logy = TRUE)

  table_page(acc_pivot, "Acceptance rate by strategy x d  (averaged over m)")
  table_page(best_tbl,
             "Most efficient strategy per (d, m)  (max median ESS/sec)")
  invisible(NULL)
}

# ------------------------------------------------------------------- MAIN ---
main <- function() {
  dir.create(CONFIG$out_dir, recursive = TRUE, showWarnings = FALSE)
  cells <- build_cells()
  cat(sprintf("Running %d cells x %d strategies (%s mode, cores=%d)...\n",
              length(cells), length(STRATEGIES),
              if (SMOKE) "SMOKE" else "FULL", CONFIG$cores))

  results <- if (CONFIG$cores > 1L) {
    parallel::mclapply(cells, run_cell, mc.cores = CONFIG$cores,
                       mc.preschedule = FALSE)
  } else {
    lapply(cells, function(cl) {
      cat(sprintf("  cell d=%d m=%d\n", cl$d, cl$m)); run_cell(cl)
    })
  }
  df <- do.call(rbind, results)

  csv_path <- file.path(CONFIG$out_dir, "zeta_block_results.csv")
  pdf_path <- file.path(CONFIG$out_dir, "zeta_block_report.pdf")
  utils::write.csv(df, csv_path, row.names = FALSE)

  verify_against_crr()
  sanity_checks(df)
  make_report(df, pdf_path)

  cat(sprintf("\nDone.\n  results: %s\n  report:  %s\n", csv_path, pdf_path))
  invisible(df)
}

if (sys.nframe() == 0) invisible(main())
