# =============================================================================
# R script to run constrained HMC for GGM precision matrix sampling via Python
# Uses subprocess approach to work with uv environments (no --enable-shared needed)
# =============================================================================

library(bgms)
library(mvtnorm)
library(jsonlite)
library(RcppCNPy)  # For reading numpy arrays
library(posterior)  # For ESS calculation

# -----------------------------------------------------------------------------
# Helper function to convert vectorized upper triangle to matrix
# -----------------------------------------------------------------------------
upper_tri_to_matrix <- function(vec, p) {
  mat <- matrix(0, p, p)
  idx <- 1
  for (j in 1:p) {
    for (i in 1:j) {
      mat[i, j] <- vec[idx]
      mat[j, i] <- vec[idx]
      idx <- idx + 1
    }
  }
  mat
}

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

# Path to the folder with the python file
python_dir <- file.path(getwd(), "dev/ggm-hmc")
# Path to the uv virtual environment
python_bin <- "dev/ggm-precision-constrained-hmc/.venv/bin/python"

# Temporary files for data exchange
tmp_dir <- tempdir()
input_file <- file.path(tmp_dir, "r_to_python_data.json")
output_file <- file.path(tmp_dir, "python_to_r_results.json")
samples_file <- file.path(tmp_dir, "python_P_samples.npy")

# -----------------------------------------------------------------------------
# 2. Generate data in R using mvtnorm
# -----------------------------------------------------------------------------

# Dimension and true precision
p <- 20

density <- 0.15

set.seed(123)
adj <- matrix(0, nrow = p, ncol = p)
adj[lower.tri(adj)] <- rbinom(p * (p - 1) / 2, size = 1, prob = density)
adj <- adj + t(adj)
# qgraph::qgraph(adj)
# sample true precision matrix from G-Wishart
Omega <- BDgraph::rgwish(1, adj = adj, b = p + sample(0:p, 1), D = diag(p))
Sigma <- solve(Omega)
zapsmall(Omega)

# Data
n <- 200
x <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)

# Identify zero elements in precision matrix (adj that don't exist)
zero_indices <- which(adj == 0 & row(adj) < col(Omega), arr.ind = TRUE)
zero_indices <- zero_indices - 1L  # Convert to 0-based indexing for Python
cat("\nNumber of zero index pairs:", nrow(zero_indices), "\n")


mh_timing <- system.time({
sampling_results <- bgms:::sample_ggm(
  inputFromR = list(X = x),
  prior_inclusion_prob = matrix(.5, p, p),
  initial_edge_indicators = adj,
  no_iter = 10000,
  no_warmup = 1000,
  no_chains = 4,
  edge_selection = FALSE,
  no_threads = 1,
  seed = 123,
  progress_type = 1
)})
cat("\nTime taken for bgms sampling:", mh_timing[3], "seconds\n")


# -----------------------------------------------------------------------------
# 3. Save data to JSON for Python
# -----------------------------------------------------------------------------

# Compute sufficient statistics for GGM
# S = X'X (scatter matrix)
S <- crossprod(x)
stopifnot(all(dim(S) == c(p, p)))

input_data <- list(
  n_variable = p,
  n_obs = as.integer(n),
  scatter_matrix = as.matrix(S),
  zero_indices = as.matrix(zero_indices),
  n_warm_up_iter = 100L,
  n_main_iter = 500L,
  n_chain = 1L,
  seed = 1234L,
  output_file = output_file,
  samples_file = samples_file
)

write_json(input_data, input_file, auto_unbox = TRUE, matrix = "rowmajor")
cat("\nData saved to:", input_file, "\n")

# -----------------------------------------------------------------------------
# 4. Python script location
# -----------------------------------------------------------------------------

python_script <- file.path(python_dir, "run_constrained_hmc.py")

# -----------------------------------------------------------------------------
# 5. Run Python script
# -----------------------------------------------------------------------------

cat("\nRunning constrained HMC sampling...\n")
cat("This may take a few minutes.\n\n")

# Run Python with the uv environment
# Use system() instead of system2() for real-time progress bar output
exit_code <- system(
  paste(shQuote(python_bin), shQuote(python_script), shQuote(input_file)),
  ignore.stdout = FALSE,
  ignore.stderr = FALSE
)

# -----------------------------------------------------------------------------
# 6. Load results back into R
# -----------------------------------------------------------------------------

if (exit_code == 0 && file.exists(output_file)) {
  results <- fromJSON(output_file)
  min_ess <- results$min_ess
  max_rhat <- results$max_rhat
  sampling_duration <- results$sampling_duration_seconds
  P_mean <- matrix(unlist(results$P_mean), nrow = p, ncol = p, byrow = TRUE)

  # Load full posterior samples from numpy file
  # Saved as 2D: (n_chain * n_iter, p * p) due to RcppCNPy limitations
  # Reshape back to 4D: (n_chain, n_iter, p, p)
  P_samples_flat <- npyLoad(samples_file)
  n_chain_loaded <- results$n_chain
  n_iter_loaded <- results$n_iter
  P_samples <- array(
    t(P_samples_flat),  # Transpose because R is column-major
    dim = c(p, p, n_iter_loaded, n_chain_loaded)
  )
  # Reorder dimensions to match expected (n_chain, n_iter, p, p)
  P_samples <- aperm(P_samples, c(4, 3, 1, 2))
  cat("\nLoaded posterior samples with shape:", dim(P_samples), "\n")

  cat("\n=== Results ===\n")
  cat(sprintf("Sampling duration: %.2f seconds\n", sampling_duration))

  cat("\nTrue precision matrix (first 5x5):\n")
  print(round(Omega[1:5, 1:5], 3))

  # -----------------------------------------------------------------------------
  # 7. Compare estimated vs true precision
  # -----------------------------------------------------------------------------

  # Check if zeros are recovered
  estimated_zeros <- which(abs(P_mean) < 0.1 & row(P_mean) != col(P_mean), arr.ind = TRUE)
  true_zeros <- which(adj == 0 & row(adj) != col(adj), arr.ind = TRUE)

  cat("\nNumber of near-zero elements in estimate:", nrow(estimated_zeros), "\n")
  cat("Number of true zero elements:", nrow(true_zeros), "\n")

  # Correlation between true and estimated
  upper_tri <- upper.tri(Omega)
  cor_precision <- cor(Omega[upper_tri], P_mean[upper_tri])
  cat(sprintf("\nCorrelation between true and estimated precision (upper triangle): %.3f\n", cor_precision))

  # -----------------------------------------------------------------------------
  # 8. Compare with bgms MH sampler results
  # -----------------------------------------------------------------------------

  cat("\n=== Comparison: HMC vs MH (bgms) ===\n")

  # Extract bgms samples - stored as vectorized upper triangle
  # samples matrix has dimensions: (p*(p+1)/2, n_iter) per chain
  bgms_samples_list <- lapply(sampling_results, function(chain) {
    if (!is.null(chain$samples)) {
      t(chain$samples)  # Transpose to (n_iter, p*(p+1)/2)
    } else {
      NULL
    }
  })
  bgms_samples_list <- Filter(Negate(is.null), bgms_samples_list)

  # Combine chains: (n_chains * n_iter, p*(p+1)/2)
  bgms_samples_combined <- do.call(rbind, bgms_samples_list)
  n_bgms_samples <- nrow(bgms_samples_combined)
  n_bgms_chains <- length(bgms_samples_list)
  n_bgms_iter <- n_bgms_samples / n_bgms_chains

  cat(sprintf("bgms: %d chains x %d iterations = %d samples\n",
              n_bgms_chains, n_bgms_iter, n_bgms_samples))
  cat(sprintf("HMC:  %d chains x %d iterations = %d samples\n",
              n_chain_loaded, n_iter_loaded, n_chain_loaded * n_iter_loaded))

  # Compute posterior means for bgms
  bgms_mean_vec <- colMeans(bgms_samples_combined)
  bgms_P_mean <- upper_tri_to_matrix(bgms_mean_vec, p)

  cat("\nbgms posterior mean precision (first 5x5):\n")
  print(round(bgms_P_mean[1:5, 1:5], 3))

  cat("\nHMC posterior mean precision (first 5x5):\n")
  print(round(P_mean[1:5, 1:5], 3))

  # -----------------------------------------------------------------------------
  # 9. Scatter plot: Compare edge estimates between methods
  # -----------------------------------------------------------------------------

  # Get upper triangle indices (off-diagonal only for edges)
  edge_idx <- which(upper.tri(Omega), arr.ind = TRUE)
  n_edges <- nrow(edge_idx)

  # Extract edge values from both methods
  hmc_edges <- sapply(1:n_edges, function(k) P_mean[edge_idx[k, 1], edge_idx[k, 2]])
  bgms_edges <- sapply(1:n_edges, function(k) bgms_P_mean[edge_idx[k, 1], edge_idx[k, 2]])
  true_edges <- sapply(1:n_edges, function(k) Omega[edge_idx[k, 1], edge_idx[k, 2]])
  true_edges2 <- sapply(1:n_edges, function(k) adj[edge_idx[k, 1], edge_idx[k, 2]])
  true_zero <- true_edges2 == 0

  # Create scatter plot comparing posterior means
  par(mfrow = c(1, 3))

  # Plot 1: HMC vs bgms posterior means
  plot(bgms_edges, hmc_edges,
       xlab = "bgms (MH) posterior mean",
       ylab = "HMC posterior mean",
       main = "Edge Estimates: HMC vs MH",
       pch = ifelse(true_zero, 1, 19),
       col = ifelse(true_zero, "red", "blue"))
  abline(0, 1, lty = 2, col = "gray")
  legend("topleft",
         legend = c("True zero", "True non-zero"),
         pch = c(1, 19), col = c("red", "blue"), bty = "n")

  # Plot 2: Both methods vs truth
  plot(true_edges, hmc_edges,
       xlab = "True precision",
       ylab = "Estimated precision",
       main = "Estimates vs Truth",
       pch = 19, col = "blue")
  points(true_edges, bgms_edges, pch = 17, col = "red")
  abline(0, 1, lty = 2, col = "gray")
  legend("topleft",
         legend = c("HMC", "bgms (MH)"),
         pch = c(19, 17), col = c("blue", "red"), bty = "n")

  bgms_incl_prob_mat <- upper_tri_to_matrix(colMeans(zapsmall(bgms_samples_combined) != 0), p)
  bgms_incl_prob <- bgms_incl_prob_mat[lower.tri(bgms_incl_prob_mat)]

  # Compute HMC inclusion probabilities from P_samples
  # P_samples has shape (n_chain, n_iter, p, p)
  hmc_incl_prob <- numeric(n_edges)
  for (k in 1:n_edges) {
    i <- edge_idx[k, 1]
    j <- edge_idx[k, 2]
    edge_samples <- as.vector(P_samples[, , i, j])  # All chains, all iters
    hmc_incl_prob[k] <- mean(zapsmall(edge_samples) != 0)
  }

  # Plot 3: Posterior mean vs inclusion probability for both methods
  xlim_range <- range(c(hmc_edges, bgms_edges))
  ylim_range <- c(0, 1)

  plot(hmc_edges, hmc_incl_prob,
       xlab = "Posterior mean (edge weight)",
       ylab = "Inclusion probability",
       main = "Posterior Mean vs Inclusion Probability",
       pch = 19, col = "blue",
       xlim = xlim_range, ylim = ylim_range)
  points(bgms_edges, bgms_incl_prob, pch = 17, col = "red")
  abline(v = 0, lty = 2, col = "gray")
  legend("topright",
         legend = c("HMC", "bgms (MH)"),
         pch = c(19, 17), col = c("blue", "red"), bty = "n")

  par(mfrow = c(1, 1))

  # Correlation between methods
  cat(sprintf("\nCorrelation between HMC and bgms edge estimates: %.4f\n",
              cor(hmc_edges, bgms_edges)))

  # -----------------------------------------------------------------------------
  # 10. ESS per second comparison
  # -----------------------------------------------------------------------------

  cat("\n=== ESS/second Comparison ===\n")

  # HMC ESS (from Python, computed by arviz)
  hmc_ess <- min_ess
  hmc_time <- sampling_duration
  hmc_ess_per_sec <- hmc_ess / hmc_time

  # bgms ESS - compute for each parameter using posterior package
  # Reshape bgms samples to draws_array format: (n_iter, n_chains, n_params)
  bgms_samples_array <- array(
    NA,
    dim = c(n_bgms_iter, n_bgms_chains, ncol(bgms_samples_combined))
  )
  for (c in 1:n_bgms_chains) {
    start_idx <- (c - 1) * n_bgms_iter + 1
    end_idx <- c * n_bgms_iter
    bgms_samples_array[, c, ] <- bgms_samples_combined[start_idx:end_idx, ]
  }

  # Compute bulk ESS for each parameter
  bgms_ess_bulk <- apply(bgms_samples_array, 3, function(x) {
    posterior::ess_bulk(x)
  })
  bgms_min_ess <- min(bgms_ess_bulk, na.rm = TRUE)
  bgms_time <- mh_timing[3]
  bgms_ess_per_sec <- bgms_min_ess / bgms_time

  print_efficiency <- function(hmc_time, hmc_ess, hmc_ess_per_sec, bgms_time, bgms_min_ess, bgms_ess_per_sec) {
    cat(sprintf("\nHMC:\n"))
    cat(sprintf("  Time: %.2f seconds\n", hmc_time))
    cat(sprintf("  Min ESS: %.1f\n", hmc_ess))
    cat(sprintf("  ESS/sec: %.2f\n", hmc_ess_per_sec))

    cat(sprintf("\nbgms (MH):\n"))
    cat(sprintf("  Time: %.2f seconds\n", bgms_time))
    cat(sprintf("  Min ESS: %.1f\n", bgms_min_ess))
    cat(sprintf("  ESS/sec: %.2f\n", bgms_ess_per_sec))

    cat(sprintf("\nEfficiency ratio (HMC / MH): %.5fx\n", hmc_ess_per_sec / bgms_ess_per_sec))
  }
  print_efficiency(hmc_time, hmc_ess, hmc_ess_per_sec, bgms_time, bgms_min_ess, bgms_ess_per_sec)

} else {
  cat("ERROR: Python script failed (exit code:", exit_code, ")\n")
  cat("Check the Python output above for errors.\n")
}
