prepare_output_bgm = function(
  out, x, num_categories, iter, data_columnnames, is_ordinal_variable,
  warmup, pairwise_scale, main_alpha, main_beta,
  na_action, na_impute, edge_selection, edge_prior, inclusion_probability,
  beta_bernoulli_alpha, beta_bernoulli_beta, beta_bernoulli_alpha_between,
  beta_bernoulli_beta_between, dirichlet_alpha, lambda,
  variable_type, update_method, target_accept, hmc_num_leapfrogs,
  nuts_max_depth, learn_mass_matrix, num_chains
) {
  arguments = list(
    prepared_data = x,
    num_variables = ncol(x),
    num_cases = nrow(x),
    na_impute = na_impute,
    variable_type = variable_type,
    iter = iter,
    warmup = warmup,
    pairwise_scale = pairwise_scale,
    main_alpha = main_alpha,
    main_beta = main_beta,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda,
    na_action = na_action,
    version = packageVersion("bgms"),
    update_method = update_method,
    target_accept = target_accept,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    num_chains = num_chains,
    num_categories = num_categories,
    data_columnnames = data_columnnames,
    no_variables = ncol(x) # backwards compatibility easybgm
  )

  num_variables = ncol(x)
  results = list()

  # ======= Parameter name generation =======
  names_variable_categories = character()
  for(v in seq_len(num_variables)) {
    if(is_ordinal_variable[v]) {
      cats = seq_len(num_categories[v])
      names_variable_categories = c(
        names_variable_categories,
        paste0(data_columnnames[v], " (", cats, ")")
      )
    } else {
      names_variable_categories = c(
        names_variable_categories,
        paste0(data_columnnames[v], " (linear)"),
        paste0(data_columnnames[v], " (quadratic)")
      )
    }
  }

  edge_names = character()
  for(i in 1:(num_variables - 1)) {
    for(j in (i + 1):num_variables) {
      edge_names = c(edge_names, paste0(data_columnnames[i], "-", data_columnnames[j]))
    }
  }

  # ======= Summarize MCMC chains =======
  summary_list = summarize_fit(out, edge_selection = edge_selection)
  main_summary = summary_list$main[, -1]
  pairwise_summary = summary_list$pairwise[, -1]

  rownames(main_summary) = names_variable_categories
  rownames(pairwise_summary) = edge_names

  results$posterior_summary_main = main_summary
  results$posterior_summary_pairwise = pairwise_summary

  if(edge_selection) {
    indicator_summary = summarize_indicator(out, param_names = edge_names)[, -1]
    rownames(indicator_summary) = edge_names
    results$posterior_summary_indicator = indicator_summary
    if(identical(edge_prior, "Stochastic-Block") && "allocations" %in% names(out[[1]])) {
      # convergence diagnostics of the co-apperance of the nodes
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(out, `[[`, "allocations"),
        node_names = data_columnnames
      )
      # posterior summary of each pairwise cluster coclustering
      results$posterior_summary_pairwise_allocations = sbm_convergence$sbm_summary
    }
  }

  # ======= Posterior mean matrices =======
  if(any(is_ordinal_variable)) {
    max_num_categories = max(2, num_categories[is_ordinal_variable])
    pmm = matrix(NA, nrow = num_variables, ncol = max_num_categories)
  } else {
    pmm = matrix(NA, nrow = num_variables, ncol = 2)
  }
  start = stop = 0
  for(v in 1:num_variables) {
    if(is_ordinal_variable[v]) {
      start = stop + 1
      stop = start + num_categories[v] - 1
      pmm[v, 1:num_categories[v]] = main_summary$mean[start:stop]
    } else {
      start = stop + 1
      stop = start + 1
      pmm[v, 1:2] = main_summary$mean[start:stop]
    }
  }

  results$posterior_mean_main = pmm
  rownames(results$posterior_mean_main) = data_columnnames
  colnames(results$posterior_mean_main) = paste0("cat (", 1:ncol(pmm), ")")

  results$posterior_mean_pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
  results$posterior_mean_pairwise[lower.tri(results$posterior_mean_pairwise)] = pairwise_summary$mean
  results$posterior_mean_pairwise = results$posterior_mean_pairwise + t(results$posterior_mean_pairwise)
  rownames(results$posterior_mean_pairwise) = data_columnnames
  colnames(results$posterior_mean_pairwise) = data_columnnames

  if(edge_selection) {
    indicator_means = indicator_summary$mean
    results$posterior_mean_indicator = matrix(0, nrow = num_variables, ncol = num_variables)
    results$posterior_mean_indicator[upper.tri(results$posterior_mean_indicator)] = indicator_means
    results$posterior_mean_indicator[lower.tri(results$posterior_mean_indicator)] =
      t(results$posterior_mean_indicator)[lower.tri(results$posterior_mean_indicator)]
    rownames(results$posterior_mean_indicator) = data_columnnames
    colnames(results$posterior_mean_indicator) = data_columnnames

    if(identical(edge_prior, "Stochastic-Block") && "allocations" %in% names(out[[1]])) {
      # convergence diagnostics of the co-apperance of the nodes
      sbm_convergence = summarize_alloc_pairs(
        allocations = lapply(out, `[[`, "allocations"),
        node_names = data_columnnames
      )
      results$posterior_mean_coclustering_matrix = sbm_convergence$co_occur_matrix
      # calculate the estimated clustering and block probabilities
      sbm_summary = posterior_summary_SBM(
        allocations = lapply(out, `[[`, "allocations"),
        arguments = arguments
      ) # check if only arguments would work
      # extract the posterior mean and median
      results$posterior_mean_allocations = sbm_summary$allocations_mean
      results$posterior_mode_allocations = sbm_summary$allocations_mode

      # extract the number of blocks and their estimated posterior probabilties
      results$posterior_num_blocks = sbm_summary$blocks
    }
  }


  results$arguments = arguments
  class(results) = "bgms"

  results$raw_samples = list(
    main = lapply(out, function(chain) chain$main_samples),
    pairwise = lapply(out, function(chain) chain$pairwise_samples),
    indicator = if(edge_selection) lapply(out, function(chain) chain$indicator_samples) else NULL,
    allocations = if(edge_selection && identical(edge_prior, "Stochastic-Block") && "allocations" %in% names(out[[1]])) lapply(out, `[[`, "allocations") else NULL,
    nchains = length(out),
    niter = nrow(out[[1]]$main_samples),
    parameter_names = list(
      main = names_variable_categories,
      pairwise = edge_names,
      indicator = if(edge_selection) edge_names else NULL,
      allocations = if(identical(edge_prior, "Stochastic-Block")) edge_names else NULL
    )
  )

  return(results)
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


prepare_output_bgmCompare = function(
  out, observations, num_categories, is_ordinal_variable,
  num_groups, group, iter, warmup,
  main_effect_indices, pairwise_effect_indices,
  data_columnnames, difference_selection,
  difference_prior, difference_selection_alpha, difference_selection_beta,
  inclusion_probability,
  pairwise_scale, difference_scale,
  update_method, target_accept, nuts_max_depth, hmc_num_leapfrogs,
  learn_mass_matrix, num_chains, projection
) {
  num_variables = ncol(observations)

  arguments = list(
    prepared_data = observations,
    num_variables = num_variables,
    num_cases = nrow(observations),
    iter = iter,
    warmup = warmup,
    pairwise_scale = pairwise_scale,
    difference_scale = difference_scale,
    difference_selection = difference_selection,
    difference_prior = difference_prior,
    difference_selection_alpha = difference_selection_alpha,
    difference_selection_beta = difference_selection_beta,
    inclusion_probability = inclusion_probability,
    version = packageVersion("bgms"),
    update_method = update_method,
    target_accept = target_accept,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    num_chains = num_chains,
    num_groups = num_groups,
    data_columnnames = data_columnnames,
    projection = projection,
    num_categories = num_categories,
    is_ordinal_variable = is_ordinal_variable,
    group = group
  )

  # --- parameter names
  names_all = generate_param_names_bgmCompare(
    data_columnnames, num_categories, is_ordinal_variable,
    num_variables, num_groups
  )

  # --- summaries
  summary_list = summarize_fit_compare(
    fit = out,
    main_effect_indices = main_effect_indices,
    pairwise_effect_indices = pairwise_effect_indices,
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
    posterior_summary_main_baseline = summary_list$main_baseline,
    posterior_summary_pairwise_baseline = summary_list$pairwise_baseline,
    posterior_summary_main_differences = summary_list$main_differences,
    posterior_summary_pairwise_differences = summary_list$pairwise_differences
  )

  if(difference_selection) {
    results$posterior_summary_indicator = summary_list$indicators
  }

  # --- posterior means (legacy-style matrices)
  # baselines
  if(any(is_ordinal_variable)) {
    max_num_categories = max(2, num_categories[is_ordinal_variable])
    pmm = matrix(NA, nrow = num_variables, ncol = max_num_categories)
  } else {
    pmm = matrix(NA, nrow = num_variables, ncol = 2)
  }
  start = stop = 0
  for(v in 1:num_variables) {
    if(is_ordinal_variable[v]) {
      start = stop + 1
      stop = start + num_categories[v] - 1
      pmm[v, 1:num_categories[v]] = summary_list$main_baseline$mean[start:stop]
    } else {
      start = stop + 1
      stop = start + 1
      pmm[v, 1:2] = summary_list$main_baseline$mean[start:stop]
    }
  }

  results$posterior_mean_main_baseline = pmm
  rownames(results$posterior_mean_main_baseline) = data_columnnames
  colnames(results$posterior_mean_main_baseline) = paste0("cat (", 1:ncol(pmm), ")")

  results$posterior_mean_pairwise_baseline = matrix(0, num_variables, num_variables)
  results$posterior_mean_pairwise_baseline[lower.tri(results$posterior_mean_pairwise_baseline)] =
    summary_list$pairwise_baseline$mean
  results$posterior_mean_pairwise_baseline = results$posterior_mean_pairwise_baseline +
    t(results$posterior_mean_pairwise_baseline)
  rownames(results$posterior_mean_pairwise_baseline) = data_columnnames
  colnames(results$posterior_mean_pairwise_baseline) = data_columnnames

  # --- raw samples (like in prepare_output_bgm)
  results$raw_samples = list(
    main = lapply(out, function(chain) chain$main_samples),
    pairwise = lapply(out, function(chain) chain$pairwise_samples),
    indicator = if(difference_selection) lapply(out, function(chain) chain$indicator_samples) else NULL,
    nchains = length(out),
    niter = nrow(out[[1]]$main_samples),
    parameter_names = names_all
  )

  results$arguments = arguments
  class(results) = c("bgmCompare")
  results
}
