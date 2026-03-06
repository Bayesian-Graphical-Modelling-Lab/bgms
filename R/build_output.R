# ==============================================================================
# build_output: assemble fit objects from bgm_spec + raw sampler output
# ==============================================================================
#
# build_output()         — thin dispatcher
# build_output_bgm()     — unified GGM + OMRF builder (normalizes raw C++
#                           output, computes MCMC summaries, assembles fit)
# build_output_compare() — Compare-specific builder
#
# Both use build_arguments() from bgm_spec.R for the $arguments list.
# ==============================================================================


# ==============================================================================
# build_output()  — dispatcher
# ==============================================================================
build_output = function(spec, raw) {
  stopifnot(inherits(spec, "bgm_spec"))

  switch(spec$model_type,
    ggm       = build_output_bgm(spec, raw),
    omrf      = build_output_bgm(spec, raw),
    mixed_mrf = build_output_mixed_mrf(spec, raw),
    compare   = build_output_compare(spec, raw),
    stop("Unknown model_type: ", spec$model_type)
  )
}


# ==============================================================================
# build_output_bgm()  — unified GGM + OMRF
# ==============================================================================
#
# The two paths share ~80% of logic. Differences:
#   1. Parameter naming: GGM uses "Var (precision)", OMRF uses "Var (k)"
#   2. Main posterior mean shape: GGM = p×1, OMRF = p×max_categories
# ==============================================================================
build_output_bgm = function(spec, raw) {
  d = spec$data
  v = spec$variables
  p = spec$prior
  s = spec$sampler

  is_continuous = v$is_continuous
  num_variables = d$num_variables
  data_columnnames = d$data_columnnames
  edge_selection = p$edge_selection
  edge_prior = p$edge_prior

  # --- Normalize raw C++ output -----------------------------------------------
  # The C++ GGM/OMRF backends return a flat `samples` matrix (params x iters)
  # via convert_results_to_list(). Split into main and pairwise components and
  # transpose to (iters x params) — the layout that MCMC summary functions
  # expect.
  if(is_continuous) {
    # GGM: samples contain the upper triangle of the precision matrix
    # (row-major). Diagonal entries are "main"; off-diagonal are "pairwise".
    diag_idx = integer(num_variables)
    offdiag_idx = integer(num_variables * (num_variables - 1L) / 2L)
    pos = 0L
    di = 0L
    oi = 0L
    for(i in seq_len(num_variables)) {
      for(j in i:num_variables) {
        pos = pos + 1L
        if(i == j) {
          di = di + 1L
          diag_idx[di] = pos
        } else {
          oi = oi + 1L
          offdiag_idx[oi] = pos
        }
      }
    }
    raw = lapply(raw, function(chain) {
      samples_t = t(chain$samples)
      res = list(
        main_samples     = samples_t[, diag_idx, drop = FALSE],
        pairwise_samples = samples_t[, offdiag_idx, drop = FALSE],
        userInterrupt    = isTRUE(chain$userInterrupt),
        chain_id         = chain$chain_id
      )
      if(!is.null(chain$indicator_samples)) {
        res$indicator_samples = t(chain$indicator_samples)[, offdiag_idx, drop = FALSE]
      }
      if(!is.null(chain$allocation_samples)) {
        res$allocations = t(chain$allocation_samples)
      }
      res
    })
  } else {
    # OMRF: the first num_thresholds params are main effects, the rest are
    # pairwise. NUTS diagnostics use bare names from C++; rename to the
    # trailing-__ convention expected by summarize_nuts_diagnostics().
    num_thresholds = spec$precomputed$num_thresholds
    raw = lapply(raw, function(chain) {
      samples_t = t(chain$samples)
      n_params = ncol(samples_t)
      res = list(
        main_samples     = samples_t[, seq_len(num_thresholds), drop = FALSE],
        pairwise_samples = samples_t[, seq(num_thresholds + 1L, n_params), drop = FALSE],
        userInterrupt    = isTRUE(chain$userInterrupt),
        chain_id         = chain$chain_id
      )
      if(!is.null(chain$indicator_samples)) {
        res$indicator_samples = t(chain$indicator_samples)
      }
      if(!is.null(chain$allocation_samples)) {
        res$allocations = t(chain$allocation_samples)
      }
      if(!is.null(chain$treedepth)) res[["treedepth__"]] = chain$treedepth
      if(!is.null(chain$divergent)) res[["divergent__"]] = chain$divergent
      if(!is.null(chain$energy)) res[["energy__"]] = chain$energy
      res
    })
  }

  # --- Parameter names --------------------------------------------------------
  if(is_continuous) {
    # GGM: one "precision" per variable
    names_main = paste0(data_columnnames, " (precision)")
    is_ordinal_variable = NULL
    num_categories = NULL
  } else {
    # OMRF: per-category thresholds or BC linear/quadratic
    is_ordinal_variable = v$is_ordinal
    num_categories = d$num_categories
    names_main = character()
    for(vi in seq_len(num_variables)) {
      if(is_ordinal_variable[vi]) {
        cats = seq_len(num_categories[vi])
        names_main = c(
          names_main,
          paste0(data_columnnames[vi], " (", cats, ")")
        )
      } else {
        names_main = c(
          names_main,
          paste0(data_columnnames[vi], " (linear)"),
          paste0(data_columnnames[vi], " (quadratic)")
        )
      }
    }
  }

  edge_names = character()
  for(i in seq_len(num_variables - 1)) {
    for(j in seq(i + 1, num_variables)) {
      edge_names = c(
        edge_names,
        paste0(data_columnnames[i], "-", data_columnnames[j])
      )
    }
  }

  # --- MCMC summaries ---------------------------------------------------------
  summary_list = summarize_fit(raw, edge_selection = edge_selection)
  main_summary = summary_list$main[, -1]
  pairwise_summary = summary_list$pairwise[, -1]

  rownames(main_summary) = names_main
  rownames(pairwise_summary) = edge_names

  results = list()
  results$posterior_summary_main = main_summary
  results$posterior_summary_pairwise = pairwise_summary

  # --- Edge selection summaries -----------------------------------------------
  has_sbm = FALSE
  if(edge_selection) {
    indicator_summary = summarize_indicator(raw, param_names = edge_names)[, -1]
    rownames(indicator_summary) = edge_names
    results$posterior_summary_indicator = indicator_summary

    has_sbm = identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])

    if(has_sbm) {
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(raw, `[[`, "allocations"),
        node_names  = data_columnnames
      )
      results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
    }
  }

  # --- Posterior mean: main ---------------------------------------------------
  if(is_continuous) {
    # GGM: p × 1 matrix
    results$posterior_mean_main = matrix(
      main_summary$mean,
      nrow     = num_variables,
      ncol     = 1,
      dimnames = list(data_columnnames, "precision_diag")
    )
  } else {
    # OMRF: p × max_categories matrix
    num_params = ifelse(is_ordinal_variable, num_categories, 2L)
    max_num_categories = max(num_params)

    pmm = matrix(NA, nrow = num_variables, ncol = max_num_categories)
    start = 0L
    stop = 0L
    for(vi in seq_len(num_variables)) {
      if(is_ordinal_variable[vi]) {
        start = stop + 1L
        stop = start + num_categories[vi] - 1L
        pmm[vi, seq_len(num_categories[vi])] = main_summary$mean[start:stop]
      } else {
        start = stop + 1L
        stop = start + 1L
        pmm[vi, 1:2] = main_summary$mean[start:stop]
      }
    }
    results$posterior_mean_main = pmm
    rownames(results$posterior_mean_main) = data_columnnames
    colnames(results$posterior_mean_main) = paste0("cat (", seq_len(ncol(pmm)), ")")
  }

  # --- Posterior mean: pairwise -----------------------------------------------
  results$posterior_mean_pairwise = matrix(0,
    nrow = num_variables, ncol = num_variables,
    dimnames = list(data_columnnames, data_columnnames)
  )
  results$posterior_mean_pairwise[lower.tri(results$posterior_mean_pairwise)] =
    pairwise_summary$mean
  results$posterior_mean_pairwise = results$posterior_mean_pairwise +
    t(results$posterior_mean_pairwise)

  # --- Posterior mean: indicator + SBM ----------------------------------------
  if(edge_selection) {
    indicator_means = indicator_summary$mean
    results$posterior_mean_indicator = matrix(0,
      nrow = num_variables, ncol = num_variables,
      dimnames = list(data_columnnames, data_columnnames)
    )
    results$posterior_mean_indicator[lower.tri(results$posterior_mean_indicator)] =
      indicator_means
    results$posterior_mean_indicator = results$posterior_mean_indicator +
      t(results$posterior_mean_indicator)

    if(has_sbm) {
      sbm_convergence2 = summarize_alloc_pairs(
        allocations = lapply(raw, `[[`, "allocations"),
        node_names  = data_columnnames
      )
      results$posterior_mean_coclustering_matrix = sbm_convergence2$co_occur_matrix

      arguments = build_arguments(spec)
      sbm_summary = posterior_summary_SBM(
        allocations = lapply(raw, `[[`, "allocations"),
        arguments   = arguments
      )
      results$posterior_mean_allocations = sbm_summary$allocations_mean
      results$posterior_mode_allocations = sbm_summary$allocations_mode
      results$posterior_num_blocks = sbm_summary$blocks
    }
  }

  # --- arguments + class ------------------------------------------------------
  results$arguments = build_arguments(spec)
  class(results) = "bgms"

  # --- raw_samples ------------------------------------------------------------
  results$raw_samples = list(
    main = lapply(raw, function(chain) chain$main_samples),
    pairwise = lapply(raw, function(chain) chain$pairwise_samples),
    indicator = if(edge_selection) {
      lapply(raw, function(chain) chain$indicator_samples)
    } else {
      NULL
    },
    allocations = if(edge_selection &&
      identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])) {
      lapply(raw, `[[`, "allocations")
    } else {
      NULL
    },
    nchains = length(raw),
    niter = nrow(raw[[1]]$main_samples),
    parameter_names = list(
      main = names_main,
      pairwise = edge_names,
      indicator = if(edge_selection) edge_names else NULL,
      allocations = if(identical(edge_prior, "Stochastic-Block")) {
        if(is_continuous) data_columnnames else edge_names
      } else {
        NULL
      }
    )
  )

  # --- easybgm compat shim (OMRF only) ---------------------------------------
  if(!is_continuous && "easybgm" %in% loadedNamespaces()) {
    ebgm_version = utils::packageVersion("easybgm")
    if(ebgm_version <= "0.2.1") {
      warning(
        "bgms is running in compatibility mode for easybgm (<= 0.2.1). ",
        "This will be removed once easybgm >= 0.2.2 is on CRAN."
      )
      results$arguments$save = TRUE
      if(edge_selection) {
        results$indicator = extract_indicators(results)
      }
      results$interactions = extract_pairwise_interactions(results)
      results$thresholds = extract_category_thresholds(results)
    }
  }

  # --- NUTS diagnostics -------------------------------------------------------
  if(s$update_method == "nuts") {
    results$nuts_diag = summarize_nuts_diagnostics(
      raw,
      nuts_max_depth = s$nuts_max_depth
    )
  }

  results
}


# ==============================================================================
# build_output_mixed_mrf()  — Mixed MRF builder
# ==============================================================================
#
# Handles the mixed discrete + continuous parameter layout:
#   C++ flat vector: [mux | Kxx_ut | muy | Kyy_ut | Kxy]
#   C++ indicators:  [Gxx_ut | Gyy_ut | Gxy]
#
# Splits into main (mux, muy, Kyy_diag) and pairwise (Kxx, Kyy_offdiag, Kxy),
# builds (p+q)×(p+q) interaction and indicator matrices, and maps internal
# variable order (discrete first, continuous second) back to original column
# order from the user's data.
# ==============================================================================
build_output_mixed_mrf = function(spec, raw) {
  d = spec$data
  v = spec$variables
  pr = spec$prior
  s = spec$sampler

  p = d$num_discrete
  q = d$num_continuous
  num_variables = d$num_variables
  data_columnnames = d$data_columnnames
  disc_names = d$data_columnnames_discrete
  cont_names = d$data_columnnames_continuous
  disc_idx = d$discrete_indices
  cont_idx = d$continuous_indices
  is_ordinal = v$is_ordinal
  num_categories = d$num_categories
  edge_selection = pr$edge_selection

  # --- Compute index layout in flat parameter vector --------------------------
  nt = spec$precomputed$num_thresholds
  nxx = as.integer(p * (p - 1) / 2)
  nyy_total = as.integer(q * (q + 1) / 2)
  nyy_offdiag = as.integer(q * (q - 1) / 2)
  nxy = as.integer(p * q)

  # Offsets in the flat vector (1-based)
  mux_start = 1L
  mux_end = nt
  kxx_start = nt + 1L
  kxx_end = nt + nxx
  muy_start = nt + nxx + 1L
  muy_end = nt + nxx + q
  kyy_start = nt + nxx + q + 1L
  kyy_end = nt + nxx + q + nyy_total
  kxy_start = nt + nxx + q + nyy_total + 1L
  kxy_end = nt + nxx + q + nyy_total + nxy

  # Kyy diagonal indices within the Kyy block
  kyy_diag_within = integer(q)
  kyy_offdiag_within = integer(nyy_offdiag)
  k_diag = 0L
  k_off = 0L
  pos = 0L
  for(i in seq_len(q)) {
    for(j in i:q) {
      pos = pos + 1L
      if(i == j) {
        k_diag = k_diag + 1L
        kyy_diag_within[k_diag] = pos
      } else {
        k_off = k_off + 1L
        kyy_offdiag_within[k_off] = pos
      }
    }
  }
  kyy_diag_abs = kyy_start - 1L + kyy_diag_within
  kyy_offdiag_abs = kyy_start - 1L + kyy_offdiag_within

  # Main indices: mux + muy + Kyy diagonal
  main_idx = c(seq(mux_start, mux_end),
               seq(muy_start, muy_end),
               kyy_diag_abs)

  # Pairwise indices: Kxx + Kyy off-diag + Kxy
  pairwise_idx = c(
    if(nxx > 0) seq(kxx_start, kxx_end) else integer(0),
    kyy_offdiag_abs,
    if(nxy > 0) seq(kxy_start, kxy_end) else integer(0)
  )

  # --- Indicator index layout -------------------------------------------------
  # C++ indicator vector: [Gxx_ut | Gyy_ut | Gxy]
  # All are pairwise, so indicator_samples maps directly to pairwise order:
  # Kxx edges, Kyy edges, Kxy edges — same order as pairwise_idx above.

  # --- Normalize raw output per chain -----------------------------------------
  raw = lapply(raw, function(chain) {
    samples_t = t(chain$samples)
    res = list(
      main_samples     = samples_t[, main_idx, drop = FALSE],
      pairwise_samples = samples_t[, pairwise_idx, drop = FALSE],
      userInterrupt    = isTRUE(chain$userInterrupt),
      chain_id         = chain$chain_id
    )
    if(!is.null(chain$indicator_samples)) {
      res$indicator_samples = t(chain$indicator_samples)
    }
    if(!is.null(chain$allocation_samples)) {
      res$allocations = t(chain$allocation_samples)
    }
    res
  })

  # --- Parameter names --------------------------------------------------------
  # Main effect names (in internal order: discrete first, continuous second)
  names_main = character()
  for(si in seq_len(p)) {
    if(is_ordinal[si]) {
      cats = seq_len(num_categories[si])
      names_main = c(names_main, paste0(disc_names[si], " (", cats, ")"))
    } else {
      names_main = c(names_main,
        paste0(disc_names[si], " (linear)"),
        paste0(disc_names[si], " (quadratic)"))
    }
  }
  for(ji in seq_len(q)) {
    names_main = c(names_main, paste0(cont_names[ji], " (mean)"))
  }
  for(ji in seq_len(q)) {
    names_main = c(names_main, paste0(cont_names[ji], " (precision)"))
  }

  # Pairwise edge names — internal order, mapped to original column names
  # We need a mapping from internal index to original variable name
  # Internal variables: [disc_1, ..., disc_p, cont_1, ..., cont_q]
  # Their original names: c(disc_names, cont_names)
  all_internal_names = c(disc_names, cont_names)

  edge_names = character()
  # Kxx edges (discrete-discrete)
  if(p > 1) {
    for(i in seq_len(p - 1)) {
      for(j in seq(i + 1, p)) {
        edge_names = c(edge_names,
          paste0(disc_names[i], "-", disc_names[j]))
      }
    }
  }
  # Kyy edges (continuous-continuous, off-diagonal)
  if(q > 1) {
    for(i in seq_len(q - 1)) {
      for(j in seq(i + 1, q)) {
        edge_names = c(edge_names,
          paste0(cont_names[i], "-", cont_names[j]))
      }
    }
  }
  # Kxy edges (discrete-continuous)
  if(p > 0 && q > 0) {
    for(i in seq_len(p)) {
      for(j in seq_len(q)) {
        edge_names = c(edge_names,
          paste0(disc_names[i], "-", cont_names[j]))
      }
    }
  }

  # --- MCMC summaries ---------------------------------------------------------
  summary_list = summarize_fit(raw, edge_selection = edge_selection)
  main_summary = summary_list$main[, -1]
  pairwise_summary = summary_list$pairwise[, -1]

  rownames(main_summary) = names_main
  rownames(pairwise_summary) = edge_names

  results = list()
  results$posterior_summary_main = main_summary
  results$posterior_summary_pairwise = pairwise_summary

  # --- Edge selection summaries -----------------------------------------------
  edge_prior = pr$edge_prior
  has_sbm = FALSE
  if(edge_selection) {
    indicator_summary = summarize_indicator(raw, param_names = edge_names)[, -1]
    rownames(indicator_summary) = edge_names
    results$posterior_summary_indicator = indicator_summary

    has_sbm = identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])

    if(has_sbm) {
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(raw, `[[`, "allocations"),
        node_names  = all_internal_names
      )
      results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
    }
  }

  # --- Posterior mean: main ---------------------------------------------------
  # Discrete main effects: p × max_cats matrix (like OMRF)
  num_params_disc = ifelse(is_ordinal, num_categories, 2L)
  max_num_cats = max(num_params_disc)
  pmm_disc = matrix(NA, nrow = p, ncol = max_num_cats)
  start = 0L
  stop = 0L
  for(si in seq_len(p)) {
    if(is_ordinal[si]) {
      start = stop + 1L
      stop = start + num_categories[si] - 1L
      pmm_disc[si, seq_len(num_categories[si])] = main_summary$mean[start:stop]
    } else {
      start = stop + 1L
      stop = start + 1L
      pmm_disc[si, 1:2] = main_summary$mean[start:stop]
    }
  }
  rownames(pmm_disc) = disc_names
  colnames(pmm_disc) = paste0("cat (", seq_len(max_num_cats), ")")

  # Continuous main effects: q × 2 matrix (mean, precision)
  pmm_cont = matrix(NA, nrow = q, ncol = 2)
  muy_means = main_summary$mean[nt + seq_len(q)]
  kyy_diag_means = main_summary$mean[nt + q + seq_len(q)]
  pmm_cont[, 1] = muy_means
  pmm_cont[, 2] = kyy_diag_means
  rownames(pmm_cont) = cont_names
  colnames(pmm_cont) = c("mean", "precision")

  results$posterior_mean_main = list(
    discrete = pmm_disc,
    continuous = pmm_cont
  )

  # --- Posterior mean: pairwise as (p+q) × (p+q) matrix -----------------------
  # Map from internal block indices to original column positions
  pmat = matrix(0, nrow = num_variables, ncol = num_variables,
    dimnames = list(data_columnnames, data_columnnames))

  pw_means = pairwise_summary$mean
  idx = 0L

  # Kxx block
  if(p > 1) {
    for(i in seq_len(p - 1)) {
      for(j in seq(i + 1, p)) {
        idx = idx + 1L
        oi = disc_idx[i]
        oj = disc_idx[j]
        pmat[oi, oj] = pw_means[idx]
        pmat[oj, oi] = pw_means[idx]
      }
    }
  }

  # Kyy off-diagonal block
  if(q > 1) {
    for(i in seq_len(q - 1)) {
      for(j in seq(i + 1, q)) {
        idx = idx + 1L
        oi = cont_idx[i]
        oj = cont_idx[j]
        pmat[oi, oj] = pw_means[idx]
        pmat[oj, oi] = pw_means[idx]
      }
    }
  }

  # Kxy block
  if(p > 0 && q > 0) {
    for(i in seq_len(p)) {
      for(j in seq_len(q)) {
        idx = idx + 1L
        oi = disc_idx[i]
        oj = cont_idx[j]
        pmat[oi, oj] = pw_means[idx]
        pmat[oj, oi] = pw_means[idx]
      }
    }
  }

  results$posterior_mean_pairwise = pmat

  # --- Posterior mean: indicator -----------------------------------------------
  if(edge_selection) {
    ind_means = indicator_summary$mean
    imat = matrix(0, nrow = num_variables, ncol = num_variables,
      dimnames = list(data_columnnames, data_columnnames))

    idx = 0L
    if(p > 1) {
      for(i in seq_len(p - 1)) {
        for(j in seq(i + 1, p)) {
          idx = idx + 1L
          oi = disc_idx[i]
          oj = disc_idx[j]
          imat[oi, oj] = ind_means[idx]
          imat[oj, oi] = ind_means[idx]
        }
      }
    }
    if(q > 1) {
      for(i in seq_len(q - 1)) {
        for(j in seq(i + 1, q)) {
          idx = idx + 1L
          oi = cont_idx[i]
          oj = cont_idx[j]
          imat[oi, oj] = ind_means[idx]
          imat[oj, oi] = ind_means[idx]
        }
      }
    }
    if(p > 0 && q > 0) {
      for(i in seq_len(p)) {
        for(j in seq_len(q)) {
          idx = idx + 1L
          oi = disc_idx[i]
          oj = cont_idx[j]
          imat[oi, oj] = ind_means[idx]
          imat[oj, oi] = ind_means[idx]
        }
      }
    }
    results$posterior_mean_indicator = imat

    if(has_sbm) {
      sbm_convergence2 = summarize_alloc_pairs(
        allocations = lapply(raw, `[[`, "allocations"),
        node_names  = all_internal_names
      )
      results$posterior_mean_coclustering_matrix = sbm_convergence2$co_occur_matrix

      arguments = build_arguments(spec)
      sbm_summary = posterior_summary_SBM(
        allocations = lapply(raw, `[[`, "allocations"),
        arguments   = arguments
      )
      results$posterior_mean_allocations = sbm_summary$allocations_mean
      results$posterior_mode_allocations = sbm_summary$allocations_mode
      results$posterior_num_blocks = sbm_summary$blocks
    }
  }

  # --- arguments + class ------------------------------------------------------
  results$arguments = build_arguments(spec)
  class(results) = "bgms"

  # --- raw_samples ------------------------------------------------------------
  results$raw_samples = list(
    main = lapply(raw, function(chain) chain$main_samples),
    pairwise = lapply(raw, function(chain) chain$pairwise_samples),
    indicator = if(edge_selection) {
      lapply(raw, function(chain) chain$indicator_samples)
    } else {
      NULL
    },
    allocations = if(edge_selection &&
      identical(edge_prior, "Stochastic-Block") &&
      "allocations" %in% names(raw[[1]])) {
      lapply(raw, `[[`, "allocations")
    } else {
      NULL
    },
    nchains = length(raw),
    niter = nrow(raw[[1]]$main_samples),
    parameter_names = list(
      main = names_main,
      pairwise = edge_names,
      indicator = if(edge_selection) edge_names else NULL
    )
  )

  results
}


# ==============================================================================
# build_output_compare()
# ==============================================================================
build_output_compare = function(spec, raw) {
  d = spec$data
  v = spec$variables
  p = spec$prior
  s = spec$sampler
  pc = spec$precomputed

  num_variables = d$num_variables
  num_groups = d$num_groups
  data_columnnames = d$data_columnnames
  num_categories = d$num_categories
  is_ordinal_variable = v$is_ordinal
  difference_selection = p$difference_selection

  # --- Parameter names --------------------------------------------------------
  names_all = generate_param_names_bgmCompare(
    data_columnnames    = data_columnnames,
    num_categories      = num_categories,
    is_ordinal_variable = is_ordinal_variable,
    num_variables       = num_variables,
    num_groups          = num_groups
  )

  # --- MCMC summaries ---------------------------------------------------------
  summary_list = summarize_fit_compare(
    fit = raw,
    main_effect_indices = pc$main_effect_indices,
    pairwise_effect_indices = pc$pairwise_effect_indices,
    num_variables = num_variables,
    num_groups = num_groups,
    difference_selection = difference_selection,
    param_names_main = names_all$main_baseline,
    param_names_pairwise = names_all$pairwise_baseline,
    param_names_main_diff = names_all$main_diff,
    param_names_pairwise_diff = names_all$pairwise_diff,
    param_names_indicators = names_all$indicators
  )

  results = list(
    posterior_summary_main_baseline        = summary_list$main_baseline,
    posterior_summary_pairwise_baseline    = summary_list$pairwise_baseline,
    posterior_summary_main_differences     = summary_list$main_differences,
    posterior_summary_pairwise_differences = summary_list$pairwise_differences
  )

  if(difference_selection) {
    results$posterior_summary_indicator = summary_list$indicators
  }

  # --- Posterior mean: main baseline ------------------------------------------
  num_params = ifelse(is_ordinal_variable, num_categories, 2L)
  max_num_categories = max(num_params, na.rm = TRUE)

  pmm = matrix(NA, nrow = num_variables, ncol = max_num_categories)
  start = 0L
  stop = 0L
  for(vi in seq_len(num_variables)) {
    if(is_ordinal_variable[vi]) {
      start = stop + 1L
      stop = start + num_categories[vi] - 1L
      pmm[vi, seq_len(num_categories[vi])] =
        summary_list$main_baseline$mean[start:stop]
    } else {
      start = stop + 1L
      stop = start + 1L
      pmm[vi, 1:2] = summary_list$main_baseline$mean[start:stop]
    }
  }
  results$posterior_mean_main_baseline = pmm
  rownames(results$posterior_mean_main_baseline) = data_columnnames
  colnames(results$posterior_mean_main_baseline) =
    paste0("cat (", seq_len(ncol(pmm)), ")")

  # --- Posterior mean: pairwise baseline --------------------------------------
  results$posterior_mean_pairwise_baseline = matrix(0,
    nrow = num_variables, ncol = num_variables,
    dimnames = list(data_columnnames, data_columnnames)
  )
  results$posterior_mean_pairwise_baseline[
    lower.tri(results$posterior_mean_pairwise_baseline)
  ] = summary_list$pairwise_baseline$mean
  results$posterior_mean_pairwise_baseline =
    results$posterior_mean_pairwise_baseline +
    t(results$posterior_mean_pairwise_baseline)

  # --- raw_samples ------------------------------------------------------------
  results$raw_samples = list(
    main = lapply(raw, function(chain) chain$main_samples),
    pairwise = lapply(raw, function(chain) chain$pairwise_samples),
    indicator = if(difference_selection) {
      lapply(raw, function(chain) chain$indicator_samples)
    } else {
      NULL
    },
    nchains = length(raw),
    niter = nrow(raw[[1]]$main_samples),
    parameter_names = names_all
  )

  # --- arguments + class ------------------------------------------------------
  results$arguments = build_arguments(spec)
  class(results) = "bgmCompare"

  # --- NUTS diagnostics -------------------------------------------------------
  if(s$update_method == "nuts") {
    results$nuts_diag = summarize_nuts_diagnostics(
      raw,
      nuts_max_depth = s$nuts_max_depth
    )
  }

  results
}


# ==============================================================================
# generate_param_names_bgmCompare()
# ==============================================================================
#
# Build parameter names for bgmCompare models. Used by build_output_compare().
#
# @param data_columnnames  Character vector: variable names.
# @param num_categories  Integer vector: max category per variable.
# @param is_ordinal_variable  Logical vector: TRUE = ordinal, FALSE = BC.
# @param num_variables  Integer: number of variables.
# @param num_groups  Integer: number of groups.
#
# Returns: named list with main_baseline, main_diff, pairwise_baseline,
#   pairwise_diff, and indicators character vectors.
# ==============================================================================
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
  generate_indicator_names = function(data_columnnames) {
    V = length(data_columnnames)
    out = character()
    for(i in seq_len(V)) {
      # main (diagonal)
      out = c(out, paste0(data_columnnames[i], " (main)"))
      # then all pairs with i as the first index
      if(i < V) {
        for(j in seq.int(i + 1L, V)) {
          out = c(out, paste0(data_columnnames[i], "-", data_columnnames[j], " (pairwise)"))
        }
      }
    }
    # optional sanity check: length must be V*(V+1)/2
    stopifnot(length(out) == V * (V + 1L) / 2L)
    out
  }
  names_indicators = generate_indicator_names(data_columnnames)

  list(
    main_baseline = names_main_baseline,
    main_diff = names_main_diff,
    pairwise_baseline = names_pairwise_baseline,
    pairwise_diff = names_pairwise_diff,
    indicators = names_indicators
  )
}
