# ==============================================================================
# run_sampler: dispatch from bgm_spec to C++ backends
# ==============================================================================
#
# Thin dispatcher that reads a validated bgm_spec and calls the appropriate
# C++ sampling function. Returns the raw per-chain output lists from C++.
# ==============================================================================


# ==============================================================================
# run_sampler()  — main dispatcher
# ==============================================================================
run_sampler = function(spec) {
  stopifnot(inherits(spec, "bgm_spec"))

  raw = switch(spec$model_type,
    ggm       = run_sampler_ggm(spec),
    omrf      = run_sampler_omrf(spec),
    mixed_mrf = run_sampler_mixed_mrf(spec),
    compare   = run_sampler_compare(spec),
    stop("Unknown model_type: ", spec$model_type)
  )

  # Check for user interrupt across all chains
  userInterrupt = any(vapply(raw, `[[`, logical(1L), "userInterrupt"))
  attr(raw, "userInterrupt") = userInterrupt
  if(userInterrupt) {
    warning("Stopped sampling after user interrupt, results are likely uninterpretable.")
  }

  raw
}


# ==============================================================================
# run_sampler_ggm()
# ==============================================================================
run_sampler_ggm = function(spec) {
  d = spec$data
  p = spec$prior
  s = spec$sampler
  m = spec$missing

  # C++ expects -1 for "no between-cluster prior"
  bb_alpha_between = if(is.null(p$beta_bernoulli_alpha_between)) {
    -1.0
  } else {
    p$beta_bernoulli_alpha_between
  }
  bb_beta_between = if(is.null(p$beta_bernoulli_beta_between)) {
    -1.0
  } else {
    p$beta_bernoulli_beta_between
  }

  out_raw = sample_ggm(
    inputFromR = list(X = d$x),
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L,
      nrow = d$num_variables,
      ncol = d$num_variables
    ),
    no_iter = s$iter,
    no_warmup = s$warmup,
    no_chains = s$chains,
    edge_selection = p$edge_selection,
    seed = s$seed,
    no_threads = s$cores,
    progress_type = s$progress_type,
    edge_prior = p$edge_prior,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda,
    na_impute = m$na_impute,
    missing_index_nullable = m$missing_index
  )

  out_raw
}


# ==============================================================================
# run_sampler_omrf()
# ==============================================================================
run_sampler_omrf = function(spec) {
  d = spec$data
  v = spec$variables
  m = spec$missing
  p = spec$prior
  s = spec$sampler

  # C++ expects -1 for "no between-cluster prior"
  bb_alpha_between = if(is.null(p$beta_bernoulli_alpha_between)) {
    -1.0
  } else {
    p$beta_bernoulli_alpha_between
  }
  bb_beta_between = if(is.null(p$beta_bernoulli_beta_between)) {
    -1.0
  } else {
    p$beta_bernoulli_beta_between
  }

  input_list = list(
    observations        = d$x,
    num_categories      = d$num_categories,
    is_ordinal_variable = v$is_ordinal,
    baseline_category   = v$baseline_category,
    main_alpha          = p$main_alpha,
    main_beta           = p$main_beta,
    pairwise_scale      = p$pairwise_scale
  )

  out_raw = sample_omrf(
    inputFromR = input_list,
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L,
      nrow = d$num_variables,
      ncol = d$num_variables
    ),
    no_iter = s$iter,
    no_warmup = s$warmup,
    no_chains = s$chains,
    no_threads = s$cores,
    progress_type = s$progress_type,
    edge_selection = p$edge_selection,
    sampler_type = s$update_method,
    seed = s$seed,
    edge_prior = p$edge_prior,
    na_impute = m$na_impute,
    missing_index_nullable = m$missing_index,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda,
    target_acceptance = s$target_accept,
    max_tree_depth = s$nuts_max_depth,
    num_leapfrogs = s$hmc_num_leapfrogs,
    pairwise_scaling_factors_nullable = p$pairwise_scaling_factors
  )

  out_raw
}


# ==============================================================================
# run_sampler_mixed_mrf()
# ==============================================================================
run_sampler_mixed_mrf = function(spec) {
  d = spec$data
  v = spec$variables
  p = spec$prior
  s = spec$sampler

  # C++ expects -1 for "no between-cluster prior"
  bb_alpha_between = if(is.null(p$beta_bernoulli_alpha_between)) {
    -1.0
  } else {
    p$beta_bernoulli_alpha_between
  }
  bb_beta_between = if(is.null(p$beta_bernoulli_beta_between)) {
    -1.0
  } else {
    p$beta_bernoulli_beta_between
  }

  input_list = list(
    discrete_observations   = d$x_discrete,
    continuous_observations = d$x_continuous,
    num_categories          = d$num_categories,
    is_ordinal_variable     = as.integer(v$is_ordinal),
    baseline_category       = v$baseline_category,
    main_alpha              = p$main_alpha,
    main_beta               = p$main_beta,
    pairwise_scale          = p$pairwise_scale,
    pseudolikelihood        = p$pseudolikelihood
  )

  out_raw = sample_mixed_mrf(
    inputFromR = input_list,
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L,
      nrow = d$num_variables,
      ncol = d$num_variables
    ),
    no_iter = s$iter,
    no_warmup = s$warmup,
    no_chains = s$chains,
    edge_selection = p$edge_selection,
    seed = s$seed,
    no_threads = s$cores,
    progress_type = s$progress_type,
    edge_prior = p$edge_prior,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda
  )

  out_raw
}


# ==============================================================================
# run_sampler_compare()
# ==============================================================================
run_sampler_compare = function(spec) {
  d = spec$data
  v = spec$variables
  m = spec$missing
  p = spec$prior
  s = spec$sampler
  pc = spec$precomputed

  run_bgmCompare_parallel(
    observations = d$x,
    num_groups = d$num_groups,
    counts_per_category = pc$counts_per_category,
    blume_capel_stats = pc$blume_capel_stats,
    pairwise_stats = pc$pairwise_stats,
    num_categories = d$num_categories,
    main_alpha = p$main_alpha,
    main_beta = p$main_beta,
    pairwise_scale = p$pairwise_scale,
    pairwise_scaling_factors = p$pairwise_scaling_factors,
    difference_scale = p$difference_scale,
    difference_selection_alpha = p$beta_bernoulli_alpha,
    difference_selection_beta = p$beta_bernoulli_beta,
    difference_prior = p$difference_prior,
    iter = s$iter,
    warmup = s$warmup,
    na_impute = m$na_impute,
    missing_data_indices = m$missing_index,
    is_ordinal_variable = v$is_ordinal,
    baseline_category = v$baseline_category,
    difference_selection = p$difference_selection,
    main_difference_selection = p$main_difference_selection,
    main_effect_indices = pc$main_effect_indices,
    pairwise_effect_indices = pc$pairwise_effect_indices,
    target_accept = s$target_accept,
    nuts_max_depth = s$nuts_max_depth,
    learn_mass_matrix = s$learn_mass_matrix,
    projection = d$projection,
    group_membership = sort(d$group) - 1L,
    group_indices = d$group_indices,
    interaction_index_matrix = pc$interaction_index_matrix,
    inclusion_probability = p$inclusion_probability_difference,
    num_chains = s$chains,
    nThreads = s$cores,
    seed = s$seed,
    update_method = s$update_method,
    hmc_num_leapfrogs = s$hmc_num_leapfrogs,
    progress_type = s$progress_type
  )
}
