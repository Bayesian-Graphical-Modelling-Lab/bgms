# ==============================================================================
# Output Utilities
#
# Contains:
#   - transform_ggm_backend_output()   : Reshape GGM C++ output
#   - transform_new_backend_output()   : Reshape OMRF C++ output
#   - generate_param_names_bgmCompare(): Build parameter names for compare models
#
# Historical note: prepare_output_bgm(), prepare_output_ggm(), and
# prepare_output_bgmCompare() lived here until Phase C of the bgm_spec
# refactor, when they were replaced by build_output.R.
# ==============================================================================


# --- REMOVED: prepare_output_bgm (replaced by build_output_bgm in build_output.R) ---
# --- REMOVED: prepare_output_ggm (replaced by build_output_bgm in build_output.R) ---
# --- REMOVED: prepare_output_bgmCompare (replaced by build_output_compare in build_output.R) ---


# Transform sample_ggm output to match the old backend format.
#
# The GGM backend returns a `samples` matrix (params x iters) where params
# are the upper triangle of the precision matrix stored row-by-row:
# (0,0), (0,1), (0,2), ..., (1,1), (1,2), ..., (p-1,p-1)
# We split these into diagonal elements ("main") and off-diagonal ("pairwise").
transform_ggm_backend_output = function(out, p) {
  # Build index maps for upper triangle (row-major, matching C++ ordering)
  diag_idx = integer(p)
  offdiag_idx = integer(p * (p - 1) / 2)
  pos = 0L
  d = 0L
  od = 0L
  for (i in seq_len(p)) {
    for (j in i:p) {
      pos = pos + 1L
      if (i == j) {
        d = d + 1L
        diag_idx[d] = pos
      } else {
        od = od + 1L
        offdiag_idx[od] = pos
      }
    }
  }

  lapply(out, function(chain) {
    samples_t = t(chain$samples)  # (params x iters) -> (iters x params)

    res = list(
      main_samples = samples_t[, diag_idx, drop = FALSE],
      pairwise_samples = samples_t[, offdiag_idx, drop = FALSE],
      userInterrupt = isTRUE(chain$userInterrupt),
      chain_id = chain$chain_id
    )

    if (!is.null(chain$indicator_samples)) {
      indic_t = t(chain$indicator_samples)  # (params x iters) -> (iters x params)
      res$indicator_samples = indic_t[, offdiag_idx, drop = FALSE]
    }

    if (!is.null(chain$allocation_samples)) {
      res$allocations = t(chain$allocation_samples)  # (variables x iters) -> (iters x variables)
    }

    res
  })
}


# Transform sample_omrf output to match the old backend format.
#
# The new backend returns a flat `samples` matrix (params x iters) containing
# all main + pairwise parameters concatenated. The old backend stores separate
# `main_samples` and `pairwise_samples` matrices (iters x params). NUTS fields
# also differ in naming. This function bridges the gap so that
# `prepare_output_bgm()` can process both backends identically.
transform_new_backend_output = function(out, num_thresholds) {
  lapply(out, function(chain) {
    samples_t = t(chain$samples)  # (params x iters) -> (iters x params)
    n_params = ncol(samples_t)

    res = list(
      main_samples = samples_t[, seq_len(num_thresholds), drop = FALSE],
      pairwise_samples = samples_t[, seq(num_thresholds + 1, n_params), drop = FALSE],
      userInterrupt = isTRUE(chain$userInterrupt),
      chain_id = chain$chain_id
    )

    if (!is.null(chain$indicator_samples)) {
      res$indicator_samples = t(chain$indicator_samples)
    }

    if (!is.null(chain$allocation_samples)) {
      res$allocations = t(chain$allocation_samples)  # (variables x iters) -> (iters x variables)
    }

    # Rename NUTS diagnostics to match old backend convention (trailing __)
    if (!is.null(chain$treedepth)) res[["treedepth__"]] = chain$treedepth
    if (!is.null(chain$divergent))  res[["divergent__"]]  = chain$divergent
    if (!is.null(chain$energy))     res[["energy__"]]     = chain$energy

    res
  })
}


# Generate names for bgmCompare parameters
generate_param_names_bgmCompare = function(
  data_columnnames,
  num_categories,
  is_ordinal_variable,
  num_variables,
  num_groups
) {
  # --- main baselines
  names_main_baseline = character()
  for(v in seq_len(num_variables)) {
    if(is_ordinal_variable[v]) {
      cats = seq_len(num_categories[v])
      names_main_baseline = c(
        names_main_baseline,
        paste0(data_columnnames[v], " (", cats, ")")
      )
    } else {
      names_main_baseline = c(
        names_main_baseline,
        paste0(data_columnnames[v], " (linear)"),
        paste0(data_columnnames[v], " (quadratic)")
      )
    }
  }

  # --- main differences
  names_main_diff = character()
  for(g in 2:num_groups) {
    for(v in seq_len(num_variables)) {
      if(is_ordinal_variable[v]) {
        cats = seq_len(num_categories[v])
        names_main_diff = c(
          names_main_diff,
          paste0(data_columnnames[v], " (diff", g - 1, "; ", cats, ")")
        )
      } else {
        names_main_diff = c(
          names_main_diff,
          paste0(data_columnnames[v], " (diff", g - 1, "; linear)"),
          paste0(data_columnnames[v], " (diff", g - 1, "; quadratic)")
        )
      }
    }
  }

  # --- pairwise baselines
  names_pairwise_baseline = character()
  for(i in 1:(num_variables - 1)) {
    for(j in (i + 1):num_variables) {
      names_pairwise_baseline = c(
        names_pairwise_baseline,
        paste0(data_columnnames[i], "-", data_columnnames[j])
      )
    }
  }

  # --- pairwise differences
  names_pairwise_diff = character()
  for(g in 2:num_groups) {
    for(i in 1:(num_variables - 1)) {
      for(j in (i + 1):num_variables) {
        names_pairwise_diff = c(
          names_pairwise_diff,
          paste0(data_columnnames[i], "-", data_columnnames[j], " (diff", g - 1, ")")
        )
      }
    }
  }

  # --- indicators
  generate_indicator_names <- function(data_columnnames) {
    V <- length(data_columnnames)
    out <- character()
    for(i in seq_len(V)) {
      # main (diagonal)
      out <- c(out, paste0(data_columnnames[i], " (main)"))
      # then all pairs with i as the first index
      if(i < V) {
        for(j in seq.int(i + 1L, V)) {
          out <- c(out, paste0(data_columnnames[i], "-", data_columnnames[j], " (pairwise)"))
        }
      }
    }
    # optional sanity check: length must be V*(V+1)/2
    stopifnot(length(out) == V * (V + 1L) / 2L)
    out
  }
  names_indicators <- generate_indicator_names(data_columnnames)

  list(
    main_baseline = names_main_baseline,
    main_diff = names_main_diff,
    pairwise_baseline = names_pairwise_baseline,
    pairwise_diff = names_pairwise_diff,
    indicators = names_indicators
  )
}
