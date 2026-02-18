# Assess warmup adequacy using energy stationarity
check_warmup_complete <- function(energy_mat) {
  nchains <- nrow(energy_mat)
  n <- ncol(energy_mat)

  if (n < 20) {
    return(list(
      warmup_incomplete = rep(FALSE, nchains),
      energy_slope = rep(NA_real_, nchains),
      slope_significant = rep(FALSE, nchains),
      ebfmi_first_half = rep(NA_real_, nchains),
      ebfmi_second_half = rep(NA_real_, nchains),
      var_ratio = rep(NA_real_, nchains)
    ))
  }

  mid <- floor(n / 2)

  results <- lapply(seq_len(nchains), function(chain) {
    energy <- energy_mat[chain, ]
    first_half <- energy[1:mid]
    second_half <- energy[(mid + 1):n]

    # Linear trend in energy
    time_idx <- seq_len(n)
    trend_lm <- stats::lm(energy ~ time_idx)
    slope <- stats::coef(trend_lm)[2]
    slope_se <- summary(trend_lm)$coefficients[2, 2]
    slope_significant <- abs(slope / slope_se) > 2.58

    # E-BFMI per half
    ebfmi_first <- mean(diff(first_half)^2) / stats::var(first_half)
    ebfmi_second <- mean(diff(second_half)^2) / stats::var(second_half)

    # Variance ratio
    var_ratio <- stats::var(first_half) / stats::var(second_half)

    # Flag if any criterion triggered
    warmup_incomplete <- slope_significant || ebfmi_first < 0.3 || var_ratio > 2.0

    list(
      warmup_incomplete = warmup_incomplete,
      energy_slope = slope,
      slope_significant = slope_significant,
      ebfmi_first_half = ebfmi_first,
      ebfmi_second_half = ebfmi_second,
      var_ratio = var_ratio
    )
  })

  list(
    warmup_incomplete = sapply(results, `[[`, "warmup_incomplete"),
    energy_slope = sapply(results, `[[`, "energy_slope"),
    slope_significant = sapply(results, `[[`, "slope_significant"),
    ebfmi_first_half = sapply(results, `[[`, "ebfmi_first_half"),
    ebfmi_second_half = sapply(results, `[[`, "ebfmi_second_half"),
    var_ratio = sapply(results, `[[`, "var_ratio")
  )
}

# Summarize NUTS diagnostics across chains
summarize_nuts_diagnostics <- function(out, nuts_max_depth = 10, verbose = TRUE) {
  nuts_chains <- Filter(function(chain) {
    all(c("treedepth__", "divergent__", "energy__") %in% names(chain))
  }, out)

  if(length(nuts_chains) == 0) {
    stop("No NUTS diagnostics found in output.")
  }

  # Combine fields into matrices (chains x iterations)
  combine_diag <- function(field) {
    do.call(rbind, lapply(nuts_chains, function(chain) as.numeric(chain[[field]])))
  }

  treedepth_mat <- combine_diag("treedepth__")
  divergent_mat <- combine_diag("divergent__")
  energy_mat <- combine_diag("energy__")

  # E-BFMI per chain
  compute_ebfmi <- function(energy) {
    mean(diff(energy)^2) / stats::var(energy)
  }
  ebfmi_per_chain <- apply(energy_mat, 1, compute_ebfmi)

  warmup_check <- check_warmup_complete(energy_mat)

  # Summaries
  n_total <- nrow(divergent_mat) * ncol(divergent_mat)
  total_divergences <- sum(divergent_mat)
  max_tree_depth_hits <- sum(treedepth_mat == nuts_max_depth)
  min_ebfmi <- min(ebfmi_per_chain)
  low_ebfmi_chains <- which(ebfmi_per_chain < 0.2)

  divergence_rate <- total_divergences / n_total
  depth_hit_rate <- max_tree_depth_hits / n_total

  if(verbose) {
    issues <- character(0)

    if (total_divergences > 0) {
      if (divergence_rate > 0.001) {
        issues <- c(issues, sprintf(
          "Divergences: %d (%.2f%%) - increase target acceptance or use adaptive-metropolis",
          total_divergences, 100 * divergence_rate))
      } else {
        issues <- c(issues, sprintf(
          "Divergences: %d (%.3f%%) - check R-hat and ESS",
          total_divergences, 100 * divergence_rate))
      }
    }

    if (max_tree_depth_hits > 0) {
      if (depth_hit_rate > 0.01) {
        issues <- c(issues, sprintf(
          "Tree depth: %d hits (%.1f%%) - consider max_depth > %d",
          max_tree_depth_hits, 100 * depth_hit_rate, nuts_max_depth))
      } else {
        issues <- c(issues, sprintf(
          "Tree depth: %d hits (%.2f%%) - check ESS",
          max_tree_depth_hits, 100 * depth_hit_rate))
      }
    }

    if (length(low_ebfmi_chains) > 0) {
      issues <- c(issues, sprintf(
        "E-BFMI: %.3f in chain%s %s - see vignette('diagnostics') for guidance",
        min_ebfmi,
        if(length(low_ebfmi_chains) > 1) "s" else "",
        paste(low_ebfmi_chains, collapse = ", ")))
    }

    if (length(issues) > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
      cat("NUTS issues:\n")
      for (issue in issues) {
        cat("  -", issue, "\n")
      }
    }
  }

  invisible(list(
    treedepth = treedepth_mat,
    divergent = divergent_mat,
    energy = energy_mat,
    ebfmi = ebfmi_per_chain,
    warmup_check = warmup_check,
    summary = list(
      total_divergences = total_divergences,
      max_tree_depth_hits = max_tree_depth_hits,
      min_ebfmi = min_ebfmi,
      warmup_incomplete = any(warmup_check$warmup_incomplete)
    )
  ))
}
