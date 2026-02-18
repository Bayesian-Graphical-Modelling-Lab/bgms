# =============================================================================
# Benchmark: Hash-map cache vs Single-entry cache for NUTS memoization
# =============================================================================
# 
# Purpose: Compare two memoization strategies in src/mcmc/mcmc_memoization.h
#   1. Hash-map cache (uses std::unordered_map)
#   2. Single-entry cache (uses memcmp on last theta)
#
# Instructions:
#   1. Install the version you want to test (modify mcmc_memoization.h)
#   2. Run: R CMD INSTALL . 
#   3. Run this script: Rscript dev/benchmark_memoization.R
#
# =============================================================================

library(bgms)
library(psych)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
N_REPLICATES <- 5
SEED <- 123
WARMUP <- 1000
ITER <- 500

# -----------------------------------------------------------------------------
# Data setup
# -----------------------------------------------------------------------------
data(bfi)
bfi_data <- bfi[, 1:25]  # 25 personality items, NO na.omit

cat("=== Memoization Benchmark ===\n")
cat("Data: psych::bfi[, 1:25]\n")
cat("Dimensions:", nrow(bfi_data), "x", ncol(bfi_data), "\n")
cat("Warmup:", WARMUP, "\n")
cat("Iterations:", ITER, "\n")
cat("Replicates:", N_REPLICATES, "\n")
cat("Seed:", SEED, "\n")
cat("==============================\n\n")

# -----------------------------------------------------------------------------
# Run benchmark
# -----------------------------------------------------------------------------
set.seed(SEED)
times <- numeric(N_REPLICATES)

for (i in 1:N_REPLICATES) {
  cat("Replicate", i, "of", N_REPLICATES, "... ")
  t0 <- Sys.time()
  fit <- bgm(bfi_data, warmup = WARMUP, iter = ITER, verbose = FALSE)
  times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(round(times[i], 1), "s\n")
}

# -----------------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------------
cat("\n=== RESULTS ===\n")
cat("Mean:", round(mean(times), 2), "s\n")
cat("SD:", round(sd(times), 2), "s\n")
cat("Min:", round(min(times), 2), "s\n")
cat("Max:", round(max(times), 2), "s\n")
cat("Times:", paste(round(times, 1), collapse = ", "), "\n")
cat("===============\n")
