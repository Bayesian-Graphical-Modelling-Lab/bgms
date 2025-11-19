summarize_nuts_diagnostics <- function(out, nuts_max_depth = 10, verbose = TRUE) {
  # Filter chains that include NUTS diagnostics
  nuts_chains <- Filter(function(chain) {
    all(c("treedepth__", "divergent__", "energy__") %in% names(chain))
  }, out)

  if(length(nuts_chains) == 0) {
    stop("No NUTS diagnostics found in output.")
  }

  # Combine per-field into matrices (chains × iterations)
  combine_diag <- function(field) {
    do.call(rbind, lapply(nuts_chains, function(chain) as.numeric(chain[[field]])))
  }

  treedepth_mat <- combine_diag("treedepth__")
  divergent_mat <- combine_diag("divergent__")
  energy_mat <- combine_diag("energy__")

  # Per-chain E-BFMI calculation
  compute_ebfmi <- function(energy) {
    diffs <- diff(energy)
    mean(diffs^2) / stats::var(energy)
  }
  ebfmi_per_chain <- apply(energy_mat, 1, compute_ebfmi)

  # Summaries
  total_divergences <- sum(divergent_mat)
  max_tree_depth_hits <- sum(treedepth_mat == nuts_max_depth)
  min_ebfmi <- min(ebfmi_per_chain)

  if(verbose) {
    cat("NUTS Diagnostics Summary:\n")
    cat("  Total divergences:        ", total_divergences, "\n")
    cat("  Max tree depth hits:      ", max_tree_depth_hits, "\n")
    cat("  Min E-BFMI across chains: ", round(min_ebfmi, 3), "\n")

    divergence_rate <- total_divergences / (nrow(divergent_mat) * ncol(divergent_mat))
    if(divergence_rate > 0.001) {
      warning(sprintf(
        "About %.3f%% of transitions ended with a divergence (%d out of %d).\n",
        100 * divergence_rate,
        total_divergences,
        nrow(divergent_mat) * ncol(divergent_mat)
      ), "Consider increasing the target acceptance rate.")
    } else if(divergence_rate > 0) {
      message(
        sprintf(
          "Note: %.3f%% of transitions ended with a divergence (%d of %d).\n",
          100 * divergence_rate,
          total_divergences,
          nrow(divergent_mat) * ncol(divergent_mat)
        ),
        "Check R-hat and effective sample size (ESS) to ensure the chains are\n",
        "mixing well."
      )
    }

    depth_hit_rate <- max_tree_depth_hits / (nrow(treedepth_mat) * ncol(treedepth_mat))
    if(depth_hit_rate > 0.01) {
      warning(paste0(
        sprintf(
          "About %.2f%% of transitions hit the maximum tree depth (%d out of %d).\n",
          100 * depth_hit_rate,
          max_tree_depth_hits,
          nrow(treedepth_mat) * ncol(treedepth_mat)
        ),
        "Consider increasing max_depth."
      ))
    } else if(depth_hit_rate > 0) {
      message(paste0(
        sprintf(
          "Note: %.2f%% of transitions hit the maximum tree depth (%d of %d).\n",
          100 * depth_hit_rate,
          max_tree_depth_hits,
          nrow(treedepth_mat) * ncol(treedepth_mat)
        ),
        "Check efficiency metrics such as effective sample size (ESS) to ensure\n",
        "sufficient exploration of the posterior."
      ))
    }


    low_ebfmi_chains <- which(ebfmi_per_chain < 0.3)
    min_ebfmi <- min(ebfmi_per_chain)

    if(length(low_ebfmi_chains) > 0) {
      warning(
        sprintf(
          "E-BFMI below 0.3 detected in %d chain(s): %s.\n",
          length(low_ebfmi_chains),
          paste(low_ebfmi_chains, collapse = ", ")
        ),
        "This suggests inefficient momentum resampling in those chains.\n",
        "Sampling efficiency may be reduced. Consider longer chains or checking convergence diagnostics."
      )
    }
  }

  # Return structured summary
  list(
    treedepth = treedepth_mat,
    divergent = divergent_mat,
    energy = energy_mat,
    ebfmi = ebfmi_per_chain,
    summary = list(
      total_divergences = total_divergences,
      max_tree_depth_hits = max_tree_depth_hits,
      min_ebfmi = min_ebfmi
    )
  )
}
