# ==============================================================================
# bgm_spec: validated intermediate specification object
# ==============================================================================
#
# Central construction point for all bgm/bgmCompare models. Three layers:
#   bgm_spec()          --- user-facing: validates inputs, assembles sub-lists
#   new_bgm_spec()      --- low-level: type/presence assertions per field
#   validate_bgm_spec() --- cross-field invariant checks
#
# The result is an S3 list of class "bgm_spec" consumed by run_sampler()
# and build_output().
# ==============================================================================


# ==============================================================================
# new_bgm_spec()  --- low-level constructor
# ==============================================================================
#
# Asserts presence and type of every field. Does NOT validate values
# (that's done upstream by the individual validators) or cross-field
# invariants (that's validate_bgm_spec).
# ==============================================================================
new_bgm_spec = function(model_type, data, variables, missing, prior,
                        sampler, precomputed = list()) {
  # --- top-level structure ---
  stopifnot(
    is.character(model_type), length(model_type) == 1L,
    model_type %in% c("ggm", "omrf", "compare", "mixed_mrf")
  )

  # --- data sub-list ---
  stopifnot(is.list(data))
  if(model_type == "mixed_mrf") {
    stopifnot(is.matrix(data$x_discrete))
    stopifnot(is.matrix(data$x_continuous))
  } else {
    stopifnot(is.matrix(data$x))
  }
  stopifnot(is.character(data$data_columnnames))
  stopifnot(is.integer(data$num_variables), length(data$num_variables) == 1L)
  stopifnot(is.integer(data$num_cases), length(data$num_cases) == 1L)

  if(model_type == "omrf" || model_type == "compare") {
    stopifnot(
      is.integer(data$num_categories),
      length(data$num_categories) == data$num_variables
    )
  }
  if(model_type == "mixed_mrf") {
    stopifnot(
      is.integer(data$num_categories),
      length(data$num_categories) == data$num_discrete
    )
  }


  if(model_type == "compare") {
    stopifnot(is.integer(data$group), length(data$group) == data$num_cases)
    stopifnot(is.integer(data$num_groups), length(data$num_groups) == 1L)
    stopifnot(is.matrix(data$group_indices))
    stopifnot(is.matrix(data$projection))
  }

  # --- variables sub-list ---
  stopifnot(is.list(variables))
  stopifnot(is.character(variables$variable_type))
  stopifnot(is.logical(variables$is_ordinal))
  stopifnot(is.logical(variables$is_continuous), length(variables$is_continuous) == 1L)
  stopifnot(is.integer(variables$baseline_category))

  # --- missing sub-list ---
  stopifnot(is.list(missing))
  stopifnot(
    is.character(missing$na_action), length(missing$na_action) == 1L,
    missing$na_action %in% c("listwise", "impute")
  )
  stopifnot(is.logical(missing$na_impute), length(missing$na_impute) == 1L)
  # missing_index can be NULL (no missing) or a matrix
  if(!is.null(missing$missing_index)) {
    stopifnot(is.matrix(missing$missing_index))
  }
  # mixed MRF uses separate indices for discrete and continuous
  if(!is.null(missing$missing_index_discrete)) {
    stopifnot(is.matrix(missing$missing_index_discrete))
  }
  if(!is.null(missing$missing_index_continuous)) {
    stopifnot(is.matrix(missing$missing_index_continuous))
  }

  # --- prior sub-list ---
  stopifnot(is.list(prior))
  # All model types now carry interaction_prior_type
  stopifnot(
    is.character(prior$interaction_prior_type),
    length(prior$interaction_prior_type) == 1L
  )
  stopifnot(is.numeric(prior$pairwise_scale), length(prior$pairwise_scale) == 1L)
  if(model_type %in% c("omrf", "compare", "mixed_mrf")) {
    stopifnot(
      is.character(prior$threshold_prior_type),
      length(prior$threshold_prior_type) == 1L
    )
    stopifnot(is.logical(prior$standardize), length(prior$standardize) == 1L)
  }
  if(model_type %in% c("omrf", "compare")) {
    stopifnot(is.matrix(prior$pairwise_scaling_factors))
  }
  if(model_type == "mixed_mrf") {
    stopifnot(is.logical(prior$standardize), length(prior$standardize) == 1L)
  }
  if(model_type %in% c("ggm", "omrf", "mixed_mrf")) {
    stopifnot(is.logical(prior$edge_selection), length(prior$edge_selection) == 1L)
    stopifnot(is.character(prior$edge_prior), length(prior$edge_prior) == 1L)
    stopifnot(is.matrix(prior$inclusion_probability))
  }
  if(model_type == "compare") {
    stopifnot(
      is.logical(prior$difference_selection),
      length(prior$difference_selection) == 1L
    )
    stopifnot(
      is.logical(prior$main_difference_selection),
      length(prior$main_difference_selection) == 1L
    )
    stopifnot(
      is.character(prior$difference_prior),
      length(prior$difference_prior) == 1L
    )
    stopifnot(
      is.numeric(prior$difference_scale),
      length(prior$difference_scale) == 1L
    )
    stopifnot(is.matrix(prior$inclusion_probability_difference))
  }

  # --- sampler sub-list ---
  stopifnot(is.list(sampler))
  stopifnot(is.character(sampler$update_method), length(sampler$update_method) == 1L)
  stopifnot(is.numeric(sampler$target_accept), length(sampler$target_accept) == 1L)
  stopifnot(is.integer(sampler$iter), length(sampler$iter) == 1L)
  stopifnot(is.integer(sampler$warmup), length(sampler$warmup) == 1L)
  stopifnot(is.integer(sampler$chains), length(sampler$chains) == 1L)
  stopifnot(is.integer(sampler$cores), length(sampler$cores) == 1L)
  stopifnot(is.integer(sampler$nuts_max_depth), length(sampler$nuts_max_depth) == 1L)
  stopifnot(is.logical(sampler$learn_mass_matrix), length(sampler$learn_mass_matrix) == 1L)
  stopifnot(is.integer(sampler$seed), length(sampler$seed) == 1L)
  stopifnot(is.integer(sampler$progress_type), length(sampler$progress_type) == 1L)

  # --- precomputed sub-list ---
  stopifnot(is.list(precomputed))

  structure(
    list(
      model_type  = model_type,
      data        = data,
      variables   = variables,
      missing     = missing,
      prior       = prior,
      sampler     = sampler,
      precomputed = precomputed
    ),
    class = "bgm_spec"
  )
}


# ==============================================================================
# validate_bgm_spec()  --- cross-field invariant checks
# ==============================================================================
validate_bgm_spec = function(spec) {
  mt = spec$model_type

  # GGM invariants
  if(mt == "ggm") {
    if(!isTRUE(spec$variables$is_continuous)) {
      stop("bgm_spec: model_type = 'ggm' requires is_continuous = TRUE.")
    }
  }

  # Compare invariants
  if(mt == "compare") {
    if(is.null(spec$data$group)) {
      stop("bgm_spec: model_type = 'compare' requires data$group.")
    }
    if(spec$data$num_groups < 2L) {
      stop("bgm_spec: model_type = 'compare' requires num_groups >= 2.")
    }
  }

  # Edge selection consistency
  if(mt %in% c("ggm", "omrf", "mixed_mrf")) {
    if(spec$prior$edge_selection && spec$prior$edge_prior == "Not Applicable") {
      stop("bgm_spec: edge_selection = TRUE but edge_prior = 'Not Applicable'.")
    }
  }

  # Scaling factors dimensions
  if(mt %in% c("omrf", "compare")) {
    nv = spec$data$num_variables
    sf = spec$prior$pairwise_scaling_factors
    if(nrow(sf) != nv || ncol(sf) != nv) {
      stop(
        "bgm_spec: pairwise_scaling_factors dimensions (",
        nrow(sf), "x", ncol(sf), ") don't match num_variables (", nv, ")."
      )
    }
  }

  # num_categories length (OMRF / compare)
  if(mt == "omrf" || mt == "compare") {
    if(length(spec$data$num_categories) != spec$data$num_variables) {
      stop("bgm_spec: num_categories length doesn't match num_variables.")
    }
  }
  if(mt == "mixed_mrf") {
    if(length(spec$data$num_categories) != spec$data$num_discrete) {
      stop("bgm_spec: num_categories length doesn't match num_discrete.")
    }
    allowed = c("adaptive-metropolis", "nuts")
    if(!(spec$sampler$update_method %in% allowed)) {
      stop(
        "bgm_spec: model_type = 'mixed_mrf' requires update_method in ",
        paste(sQuote(allowed), collapse = " or "), ". Got '",
        spec$sampler$update_method, "'."
      )
    }
  }

  invisible(spec)
}


# ==============================================================================
# bgm_spec()  --- user-facing constructor
# ==============================================================================
#
# Validates all user inputs via dedicated validators, assembles sub-lists,
# and passes through new_bgm_spec() and validate_bgm_spec().
#
# Parameters mirror the union of bgm() and bgmCompare() arguments.
# ==============================================================================
bgm_spec = function(x,
                    model_type = c("omrf", "ggm", "compare", "mixed_mrf"),
                    # Variable specification
                    variable_type = "ordinal",
                    baseline_category = 0L,
                    # Data (compare-specific)
                    y = NULL,
                    group_indicator = NULL,
                    # Missing data
                    na_action = c("listwise", "impute"),
                    # Priors (new: prior objects unpacked by bgm())
                    interaction_prior_type = "cauchy",
                    pairwise_scale = 1,
                    interaction_alpha = NA_real_,
                    interaction_beta = NA_real_,
                    threshold_prior_type = "beta-prime",
                    main_alpha = 0.5,
                    main_beta = 0.5,
                    threshold_scale = NA_real_,
                    means_prior_type = "normal",
                    means_scale = 1,
                    means_alpha = NA_real_,
                    means_beta = NA_real_,
                    scale_prior_type = "gamma",
                    scale_shape = 1,
                    scale_rate = 1,
                    delta = NULL,
                    prior_factorization = NULL,
                    standardize = FALSE,
                    edge_selection = TRUE,
                    edge_prior = bernoulli_prior(0.5),
                    # Legacy edge prior params (accepted for backward compat)
                    inclusion_probability = 0.5,
                    beta_bernoulli_alpha_between = 1,
                    beta_bernoulli_beta_between = 1,
                    dirichlet_alpha = 1,
                    lambda = 1,
                    # Priors (compare-specific)
                    difference_selection = TRUE,
                    main_difference_selection = FALSE,
                    difference_prior = c(
                      "Bernoulli", "Beta-Bernoulli", "Stochastic-Block"
                    ),
                    difference_scale = 2.5,
                    difference_probability = 0.5,
                    # Compare difference prior hyperparameters
                    beta_bernoulli_alpha = 1,
                    beta_bernoulli_beta = 1,
                    difference_beta_bernoulli_alpha_between = 1,
                    difference_beta_bernoulli_beta_between = 1,
                    difference_dirichlet_alpha = 1,
                    difference_lambda = 1,
                    # Sampler
                    update_method = c(
                      "nuts",
                      "adaptive-metropolis"
                    ),
                    target_accept = NULL,
                    iter = 10000L,
                    warmup = 1000L,
                    nuts_max_depth = 10L,
                    learn_mass_matrix = TRUE,
                    chains = 4L,
                    cores = parallel::detectCores(),
                    seed = NULL,
                    display_progress = c("per-chain", "total", "none"),
                    verbose = TRUE,
                    progress_callback = NULL) {
  model_type = match.arg(model_type)
  na_action = tryCatch(match.arg(na_action), error = function(e) {
    stop(paste0(
      "The na_action argument should be one of \"listwise\" or \"impute\", not \"",
      na_action, "\"."
    ), call. = FALSE)
  })

  # --- Data validation --------------------------------------------------------
  x = data_check(x, "x")
  data_columnnames = if(is.null(colnames(x))) {
    paste0("Variable ", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
  num_variables = ncol(x)

  # --- Variable types ---------------------------------------------------------
  allow_continuous = (model_type != "compare")
  vt = validate_variable_types(
    variable_type    = variable_type,
    num_variables    = num_variables,
    allow_continuous = allow_continuous,
    allow_mixed      = (model_type != "compare"),
    caller           = if(model_type == "compare") "bgmCompare" else "bgm"
  )
  variable_type = vt$variable_type
  is_ordinal = vt$variable_bool
  is_continuous = vt$is_continuous
  is_mixed = vt$is_mixed

  # Resolve model_type if "omrf" default was kept but data is continuous
  if(model_type == "omrf" && is_continuous) {
    model_type = "ggm"
  }
  if(model_type == "omrf" && is_mixed) {
    model_type = "mixed_mrf"
  }

  # Auto-resolve delta = NULL to the dimension-adaptive default 0.5 * log(p),
  # where p is the dimension of the continuous precision matrix. For models
  # without a continuous block (omrf, compare) the tilt has no target, so
  # NULL resolves to 0.
  if(is.null(delta)) {
    delta = if(model_type == "ggm") {
      0.5 * log(max(num_variables, 1))
    } else if(model_type == "mixed_mrf") {
      0.5 * log(max(sum(!is_ordinal), 1))
    } else {
      0
    }
  }

  # Validate determinant-tilt exponent and reject for pure-ordinal models
  if(!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
    !is.finite(delta) || delta < 0) {
    stop("'delta' must be a single finite non-negative numeric, or NULL.")
  }
  # Resolve prior_factorization default. Pick "hierarchical" when it is
  # actually applicable: model has a continuous precision block AND the
  # encompassing slab/diag prior pair matches what the L-space SD step
  # supports. The SD path supports a Normal slab directly and a Cauchy
  # slab via the scale-mixture-of-normals representation
  # (K_ij | omega ~ N(0, sigma^2 omega), omega ~ IG(1/2, 1/2)); both
  # require a Gamma diagonal. Otherwise fall back to "joint" so the
  # historical Roverato path keeps working. omrf/compare have no
  # Z(Gamma) factor, so "joint" is the only valid path there anyway.
  if(is.null(prior_factorization)) {
    prior_factorization = if(model_type %in% c("ggm", "mixed_mrf") &&
      interaction_prior_type %in% c("normal", "cauchy") &&
      scale_prior_type == "gamma") {
      "hierarchical"
    } else {
      "joint"
    }
  }
  if(!(length(prior_factorization) == 1L &&
    is.character(prior_factorization) &&
    prior_factorization %in% c("joint", "hierarchical"))) {
    stop("'prior_factorization' must be \"joint\" or \"hierarchical\".")
  }
  if(prior_factorization == "hierarchical" &&
    !(model_type %in% c("ggm", "mixed_mrf"))) {
    stop(
      "'prior_factorization = \"hierarchical\"' requires continuous data; ",
      "the current model_type is '", model_type, "', which has no ",
      "continuous precision block. Use \"joint\" or supply continuous data."
    )
  }
  if(prior_factorization == "hierarchical" &&
    !(interaction_prior_type %in% c("normal", "cauchy"))) {
    stop(
      "'prior_factorization = \"hierarchical\"' requires a Normal or Cauchy ",
      "slab prior (interaction_prior_type in {\"normal\", \"cauchy\"}). ",
      "Re-fit with interaction_prior = normal_prior(scale = ...) or ",
      "cauchy_prior(scale = ...)."
    )
  }
  if(prior_factorization == "hierarchical" &&
    scale_prior_type != "gamma") {
    stop(
      "'prior_factorization = \"hierarchical\"' requires a Gamma diagonal ",
      "prior (scale_prior_type = \"gamma\")."
    )
  }
  if(delta > 0 && model_type %in% c("omrf", "compare")) {
    stop(
      "'delta' (determinant tilt) requires continuous variables; the ",
      "current model_type is '", model_type, "', which has no precision ",
      "matrix to tilt. Pass delta = 0 or use continuous data."
    )
  }

  # --- Sampler (needs is_continuous and edge_selection early) ------------------
  sampler = validate_sampler(
    update_method = update_method,
    target_accept = target_accept,
    iter = iter,
    warmup = warmup,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains = chains,
    cores = cores,
    seed = seed,
    display_progress = display_progress,
    is_continuous = is_continuous,
    edge_selection = if(model_type == "compare") FALSE else edge_selection,
    verbose = verbose,
    progress_callback = progress_callback
  )

  # --- Resolve edge prior object -----------------------------------------------
  if(inherits(edge_prior, "bgms_indicator_prior")) {
    ep_flat = unpack_indicator_prior(edge_prior, num_variables)
  } else if(is.character(edge_prior)) {
    # Legacy string path (tests and bgmCompare may call bgm_spec directly)
    edge_prior_str = match.arg(edge_prior,
      choices = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
    )
    ep_flat = validate_edge_prior(
      edge_selection = edge_selection, edge_prior = edge_prior_str,
      inclusion_probability = inclusion_probability,
      num_variables = num_variables,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha, lambda = lambda
    )
    ep_flat$beta_bernoulli_alpha = beta_bernoulli_alpha
    ep_flat$beta_bernoulli_beta = beta_bernoulli_beta
    ep_flat$beta_bernoulli_alpha_between = beta_bernoulli_alpha_between
    ep_flat$beta_bernoulli_beta_between = beta_bernoulli_beta_between
    ep_flat$dirichlet_alpha = dirichlet_alpha
    ep_flat$lambda = lambda
  } else {
    stop(
      "'edge_prior' must be a bgms_indicator_prior object.",
      " Use bernoulli_prior(), beta_bernoulli_prior(), or sbm_prior()."
    )
  }
  # Override edge_selection if explicitly FALSE
  if(!edge_selection) {
    ep_flat$edge_selection = FALSE
    ep_flat$edge_prior = "Not Applicable"
    ep_flat$inclusion_probability = matrix(0.5, nrow = 1, ncol = 1)
  }

  # Joint factorization with a hyperprior on inclusion probabilities is
  # not well-defined: the marginal on Gamma under "joint" is
  # pi(Gamma) * Z(Gamma), and the intractable G-Wishart normalising
  # constant Z(Gamma) propagates into the conditional update for the
  # inclusion hyperparameters (Beta-Bernoulli's p, the SBM block-level
  # probabilities). Hierarchical factors out Z(Gamma) cleanly. Only
  # matters for models with a continuous precision block.
  if(edge_selection &&
    prior_factorization == "joint" &&
    model_type %in% c("ggm", "mixed_mrf") &&
    ep_flat$edge_prior %in% c("Beta-Bernoulli", "Stochastic-Block")) {
    warning(
      "'prior_factorization = \"joint\"' combined with a ", ep_flat$edge_prior,
      " inclusion prior is not properly defined: under the joint ",
      "factorization the marginal on Gamma is pi(Gamma) * Z(Gamma), ",
      "so the intractable G-Wishart normalising constant Z(Gamma) ",
      "propagates into the likelihood for the inclusion hyperparameters. ",
      "Use 'prior_factorization = \"hierarchical\"' for a clean ",
      "factorization."
    )
  }

  # --- Build by model type ----------------------------------------------------
  if(model_type == "ggm") {
    spec = build_spec_ggm(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = as.integer(rep(0L, num_variables)),
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      scale_prior_type = scale_prior_type,
      scale_shape = scale_shape,
      scale_rate = scale_rate,
      delta = delta,
      prior_factorization = prior_factorization,
      edge_prior_flat = ep_flat
    )
  } else if(model_type == "mixed_mrf") {
    spec = build_spec_mixed_mrf(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha, main_beta = main_beta,
      threshold_scale = threshold_scale,
      means_prior_type = means_prior_type,
      means_scale = means_scale,
      means_alpha = means_alpha,
      means_beta = means_beta,
      scale_prior_type = scale_prior_type,
      scale_shape = scale_shape,
      scale_rate = scale_rate,
      delta = delta,
      standardize = standardize,
      edge_prior_flat = ep_flat
    )
  } else if(model_type == "omrf") {
    spec = build_spec_omrf(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha, main_beta = main_beta,
      threshold_scale = threshold_scale,
      standardize = standardize,
      edge_prior_flat = ep_flat
    )
  } else {
    spec = build_spec_compare(
      x = x, y = y, group_indicator = group_indicator,
      data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha, main_beta = main_beta,
      threshold_scale = threshold_scale,
      standardize = standardize,
      difference_selection = difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = difference_prior,
      difference_scale = difference_scale,
      difference_probability = difference_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = difference_beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = difference_beta_bernoulli_beta_between,
      dirichlet_alpha = difference_dirichlet_alpha,
      lambda = difference_lambda
    )
  }

  validate_bgm_spec(spec)
}


# ==============================================================================
# Internal builders (one per model type)
# ==============================================================================

build_spec_ggm = function(x, data_columnnames, num_variables,
                          variable_type, is_ordinal, is_continuous,
                          baseline_category,
                          na_action, sampler,
                          interaction_prior_type, pairwise_scale,
                          interaction_alpha, interaction_beta,
                          scale_prior_type, scale_shape, scale_rate,
                          delta = 0,
                          prior_factorization = "hierarchical",
                          edge_prior_flat) {
  # Missing data
  md = validate_missing_data(
    x = x, na_action = na_action,
    is_continuous = TRUE
  )
  x = md$x

  # Center continuous data (GGM likelihood assumes zero mean)
  x = center_continuous_data(x)

  ep = edge_prior_flat

  new_bgm_spec(
    model_type = "ggm",
    data = list(
      x                = x,
      data_columnnames = data_columnnames,
      num_variables    = as.integer(ncol(x)),
      num_cases        = as.integer(nrow(x))
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = is_ordinal,
      is_continuous     = TRUE,
      baseline_category = baseline_category
    ),
    missing = list(
      na_action     = na_action,
      na_impute     = md$na_impute,
      missing_index = md$missing_index
    ),
    prior = list(
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      scale_prior_type = scale_prior_type,
      scale_shape = scale_shape,
      scale_rate = scale_rate,
      delta = delta,
      prior_factorization = prior_factorization,
      edge_selection = ep$edge_selection,
      edge_prior = ep$edge_prior,
      inclusion_probability = ep$inclusion_probability,
      beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
      beta_bernoulli_beta = ep$beta_bernoulli_beta,
      beta_bernoulli_alpha_between = ep$beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = ep$beta_bernoulli_beta_between,
      dirichlet_alpha = ep$dirichlet_alpha,
      lambda = ep$lambda
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list()
  )
}


build_spec_omrf = function(x, data_columnnames, num_variables,
                           variable_type, is_ordinal, is_continuous,
                           baseline_category,
                           na_action, sampler,
                           interaction_prior_type, pairwise_scale,
                           interaction_alpha, interaction_beta,
                           threshold_prior_type, main_alpha, main_beta,
                           threshold_scale,
                           standardize,
                           edge_prior_flat) {
  # Baseline category
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !identical(baseline_category, 0L),
    x = x,
    variable_bool = is_ordinal
  )

  # Missing data + ordinal recoding
  md = validate_missing_data(
    x = x, na_action = na_action,
    is_continuous = FALSE
  )
  x_clean = md$x
  ord = reformat_ordinal_data(
    x = x_clean, is_ordinal = is_ordinal,
    baseline_category = bc
  )
  x_recoded = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  missing_index = md$missing_index

  ep = edge_prior_flat

  # Scaling factors
  varnames = if(is.null(colnames(x))) {
    paste0("Variable ", seq_len(num_variables))
  } else {
    colnames(x)
  }
  psf = compute_scaling_factors(
    num_variables     = num_variables,
    is_ordinal        = is_ordinal,
    num_categories    = num_categories,
    baseline_category = bc_final,
    standardize       = standardize,
    varnames          = varnames
  )

  num_thresholds = sum(ifelse(is_ordinal, num_categories, 2L))

  new_bgm_spec(
    model_type = "omrf",
    data = list(
      x                = x_recoded,
      data_columnnames = data_columnnames,
      num_variables    = as.integer(num_variables),
      num_cases        = as.integer(nrow(x_recoded)),
      num_categories   = as.integer(num_categories)
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = is_ordinal,
      is_continuous     = FALSE,
      baseline_category = as.integer(bc_final)
    ),
    missing = list(
      na_action     = na_action,
      na_impute     = md$na_impute,
      missing_index = missing_index
    ),
    prior = list(
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha,
      main_beta = main_beta,
      threshold_scale = threshold_scale,
      standardize = standardize,
      pairwise_scaling_factors = psf,
      edge_selection = ep$edge_selection,
      edge_prior = ep$edge_prior,
      inclusion_probability = ep$inclusion_probability,
      beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
      beta_bernoulli_beta = ep$beta_bernoulli_beta,
      beta_bernoulli_alpha_between = ep$beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = ep$beta_bernoulli_beta_between,
      dirichlet_alpha = ep$dirichlet_alpha,
      lambda = ep$lambda
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list(
      num_thresholds = as.integer(num_thresholds)
    )
  )
}


# ------------------------------------------------------------------
# build_spec_mixed_mrf
# ------------------------------------------------------------------
# Builds a bgm_spec for the mixed MRF model (discrete + continuous).
# Splits the input data matrix into discrete and continuous parts,
# validates and recodes discrete variables (ordinal/BC), and assembles
# the spec with metadata needed by sample_mixed_mrf() and
# build_output_mixed_mrf().
# ------------------------------------------------------------------
build_spec_mixed_mrf = function(x, data_columnnames, num_variables,
                                variable_type, is_ordinal,
                                baseline_category,
                                na_action, sampler,
                                interaction_prior_type, pairwise_scale,
                                interaction_alpha, interaction_beta,
                                threshold_prior_type, main_alpha, main_beta,
                                threshold_scale,
                                means_prior_type, means_scale,
                                means_alpha, means_beta,
                                scale_prior_type, scale_shape, scale_rate,
                                delta = 0,
                                standardize,
                                edge_prior_flat) {
  # Identify discrete vs continuous columns
  cont_idx = which(variable_type == "continuous")
  disc_idx = which(variable_type != "continuous")
  p = length(disc_idx)
  q = length(cont_idx)

  # Split data
  x_disc = x[, disc_idx, drop = FALSE]
  x_cont = x[, cont_idx, drop = FALSE]

  # Ensure integer matrix for discrete data
  storage.mode(x_disc) = "integer"
  # Ensure numeric matrix for continuous data
  storage.mode(x_cont) = "double"

  # Discrete variable properties (subset to discrete columns)
  is_ordinal_disc = is_ordinal[disc_idx]
  vtype_disc = variable_type[disc_idx]

  # Subset baseline_category to discrete columns when the user supplies a

  # full-length vector (one entry per variable, including continuous ones).
  if(length(baseline_category) == num_variables && num_variables != p) {
    baseline_category = baseline_category[disc_idx]
  }

  # Baseline category for discrete variables
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !identical(baseline_category, 0L),
    x = x_disc,
    variable_bool = is_ordinal_disc
  )

  # Missing data handling
  na_impute = FALSE
  missing_index_discrete = NULL
  missing_index_continuous = NULL

  if(na_action == "listwise") {
    missing_rows = apply(x_disc, 1, anyNA) | apply(x_cont, 1, anyNA)
    if(all(missing_rows)) {
      stop(paste0(
        "All rows in x contain at least one missing response.\n",
        "You could try option na_action = \"impute\"."
      ))
    }
    n_removed = sum(missing_rows)
    if(n_removed > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
      n_remaining = nrow(x_disc) - n_removed
      message(
        n_removed, " row", if(n_removed > 1) "s" else "",
        " with missing values excluded (n = ", n_remaining, " remaining).\n",
        "To impute missing values instead, use na_action = \"impute\"."
      )
    }
    x_disc = x_disc[!missing_rows, , drop = FALSE]
    x_cont = x_cont[!missing_rows, , drop = FALSE]
    if(nrow(x_disc) < 2) {
      stop(paste0(
        "After removing missing observations from the input matrix x,\n",
        "there were less than two rows left in x."
      ))
    }
  } else {
    # Impute path: handle discrete and continuous sub-matrices separately
    md_disc = handle_impute(x_disc)
    md_cont = handle_impute(x_cont)
    x_disc = md_disc$x
    x_cont = md_cont$x
    na_impute = md_disc$na_impute || md_cont$na_impute
    if(md_disc$na_impute) missing_index_discrete = md_disc$missing_index
    if(md_cont$na_impute) missing_index_continuous = md_cont$missing_index
  }

  # Ordinal recoding (reformat discrete data)
  ord = reformat_ordinal_data(
    x = x_disc, is_ordinal = is_ordinal_disc,
    baseline_category = bc
  )
  x_disc_recoded = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  ep = edge_prior_flat

  num_thresholds = sum(ifelse(is_ordinal_disc, num_categories, 2L))

  new_bgm_spec(
    model_type = "mixed_mrf",
    data = list(
      x_discrete = x_disc_recoded,
      x_continuous = x_cont,
      data_columnnames = data_columnnames,
      data_columnnames_discrete = data_columnnames[disc_idx],
      data_columnnames_continuous = data_columnnames[cont_idx],
      num_variables = as.integer(num_variables),
      num_discrete = as.integer(p),
      num_continuous = as.integer(q),
      num_cases = as.integer(nrow(x_disc_recoded)),
      num_categories = as.integer(num_categories),
      discrete_indices = disc_idx,
      continuous_indices = cont_idx
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = is_ordinal_disc,
      is_continuous     = FALSE,
      is_mixed          = TRUE,
      baseline_category = as.integer(bc_final)
    ),
    missing = list(
      na_action = na_action,
      na_impute = na_impute,
      missing_index_discrete = missing_index_discrete,
      missing_index_continuous = missing_index_continuous
    ),
    prior = list(
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha,
      main_beta = main_beta,
      threshold_scale = threshold_scale,
      means_prior_type = means_prior_type,
      means_scale = means_scale,
      means_alpha = means_alpha,
      means_beta = means_beta,
      scale_prior_type = scale_prior_type,
      scale_shape = scale_shape,
      scale_rate = scale_rate,
      delta = delta,
      standardize = standardize,
      edge_selection = ep$edge_selection,
      edge_prior = ep$edge_prior,
      inclusion_probability = ep$inclusion_probability,
      beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
      beta_bernoulli_beta = ep$beta_bernoulli_beta,
      beta_bernoulli_alpha_between = ep$beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = ep$beta_bernoulli_beta_between,
      dirichlet_alpha = ep$dirichlet_alpha,
      lambda = ep$lambda
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list(
      num_thresholds = as.integer(num_thresholds)
    )
  )
}


build_spec_compare = function(x, y, group_indicator,
                              data_columnnames, num_variables,
                              variable_type, is_ordinal, is_continuous,
                              baseline_category,
                              na_action, sampler,
                              interaction_prior_type, pairwise_scale,
                              interaction_alpha, interaction_beta,
                              threshold_prior_type, main_alpha, main_beta,
                              threshold_scale,
                              standardize,
                              difference_selection, main_difference_selection,
                              difference_prior,
                              difference_scale, difference_probability,
                              beta_bernoulli_alpha, beta_bernoulli_beta,
                              beta_bernoulli_alpha_between = 1,
                              beta_bernoulli_beta_between = 1,
                              dirichlet_alpha = 1,
                              lambda = 1) {
  # --- Combine x/y and create group vector ------------------------------------
  if(!is.null(y)) {
    y = data_check(y, "y")
    if(ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }
  if(is.null(y) && is.null(group_indicator)) {
    stop(paste0(
      "For multi-group designs, the bgmCompare function requires input for\n",
      "either y (group 2 data) or group_indicator (group indicator)."
    ))
  }

  if(!is.null(group_indicator)) {
    group_indicator = as.vector(group_indicator)
    if(anyNA(group_indicator)) {
      stop("group_indicator cannot contain missing values.")
    }
    if(length(group_indicator) != nrow(x)) {
      stop("Length of group_indicator must match number of rows in x.")
    }

    unique_g = unique(group_indicator)
    if(length(unique_g) == 0L) {
      stop("The bgmCompare function expects at least two groups, but the input group_indicator contains no group value.")
    }
    if(length(unique_g) == 1L) {
      stop("The bgmCompare function expects at least two groups, but the input group_indicator contains only one group value.")
    }
    if(length(unique_g) == length(group_indicator)) {
      stop("The input group_indicator contains only unique group values.")
    }

    group = group_indicator
    for(u in unique_g) {
      group[group_indicator == u] = which(unique_g == u)
    }
    tab = tabulate(group)
    if(any(tab < 2L)) {
      stop("One or more groups only had one member in the input group_indicator.")
    }
  } else {
    group = c(rep.int(1L, nrow(x)), rep.int(2L, nrow(y)))
    x = rbind(x, y)
  }

  num_variables = ncol(x)

  # --- Baseline category (needs combined x) -----------------------------------
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !identical(baseline_category, 0L),
    x = x,
    variable_bool = is_ordinal
  )

  # --- Difference prior -------------------------------------------------------
  dp = validate_difference_prior(
    difference_selection = difference_selection,
    difference_prior = difference_prior,
    difference_probability = difference_probability,
    num_variables = num_variables,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda
  )

  # --- Missing data (compare path) --------------------------------------------
  md = validate_missing_data(
    x             = x,
    na_action     = na_action,
    is_continuous = FALSE,
    group         = group
  )
  x = md$x
  na_impute = md$na_impute
  missing_index = md$missing_index
  group = md$group

  # Post-listwise group validation (bgmCompare-specific) -----------------------
  if(na_action == "listwise" && md$n_removed > 0) {
    unique_g = unique(group)
    if(length(unique_g) == length(group)) {
      stop(paste0(
        "After rows with missing observations were excluded, there were no groups, as \n",
        "there were only unique values in the input g left."
      ))
    }
    if(length(unique_g) == 1) {
      stop(paste0(
        "After rows with missing observations were excluded, there were no groups, as \n",
        "there was only one value in the input g left."
      ))
    }
    g = group
    for(u in unique_g) {
      group[g == u] = which(unique_g == u)
    }
    tab = tabulate(group)
    if(any(tab < 2)) {
      stop(paste0(
        "After rows with missing observations were excluded, one or more groups, only \n",
        "had one member in the input g."
      ))
    }
  }

  # --- Ordinal recoding (compare path) ----------------------------------------
  ord = reformat_ordinal_data(
    x                 = x,
    is_ordinal        = is_ordinal,
    baseline_category = bc
  )
  x = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  # --- Collapse categories across groups (compare-specific) -------------------
  col = collapse_categories_across_groups(
    x                 = x,
    group             = group,
    is_ordinal        = is_ordinal,
    num_categories    = num_categories,
    baseline_category = bc_final
  )
  x_recoded = col$x
  num_categories = col$num_categories
  bc_final = col$baseline_category
  ordinal_variable = is_ordinal

  num_variables = ncol(x_recoded)
  num_groups = length(unique(group))

  # Compute precomputed structures
  counts_per_category = compute_counts_per_category(
    x_recoded, num_categories, group
  )
  blume_capel_stats = compute_blume_capel_stats(
    x_recoded, bc_final, ordinal_variable, group
  )

  # Center BC variables for pairwise stats
  x_centered = x_recoded
  for(i in which(!ordinal_variable)) {
    x_centered[, i] = x_centered[, i] - bc_final[i]
  }
  pairwise_stats = compute_pairwise_stats(x_centered, group)

  # Index structures
  num_interactions = as.integer(num_variables * (num_variables - 1) / 2)

  main_effect_indices = matrix(NA_integer_, nrow = num_variables, ncol = 2)
  for(variable in seq_len(num_variables)) {
    if(variable > 1) {
      main_effect_indices[variable, 1] = 1L + main_effect_indices[variable - 1, 2]
    } else {
      main_effect_indices[variable, 1] = 0L
    }
    if(ordinal_variable[variable]) {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] +
        num_categories[variable] - 1L
    } else {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + 1L
    }
  }

  pairwise_effect_indices = matrix(NA_integer_,
    nrow = num_variables, ncol = num_variables
  )
  tel = 0L
  for(v1 in seq_len(num_variables - 1)) {
    for(v2 in seq(v1 + 1, num_variables)) {
      pairwise_effect_indices[v1, v2] = tel
      pairwise_effect_indices[v2, v1] = tel
      tel = tel + 1L
    }
  }

  # Interaction index matrix (used by C++ to iterate edges in random order)
  interaction_index_matrix = matrix(0L, nrow = num_interactions, ncol = 3)
  counter = 0L
  for(v1 in seq_len(num_variables - 1)) {
    for(v2 in seq(v1 + 1, num_variables)) {
      counter = counter + 1L
      interaction_index_matrix[counter, ] = c(counter, v1 - 1L, v2 - 1L)
    }
  }

  # Scaling factors
  varnames = if(is.null(colnames(x_recoded))) {
    paste0("Variable ", seq_len(num_variables))
  } else {
    colnames(x_recoded)
  }
  psf = compute_scaling_factors(
    num_variables     = num_variables,
    is_ordinal        = ordinal_variable,
    num_categories    = num_categories,
    baseline_category = bc_final,
    standardize       = standardize,
    varnames          = varnames
  )

  # Group indices and projection
  group_indices = matrix(NA_integer_, nrow = num_groups, ncol = 2)
  observations = x_centered
  sorted_group = sort(group)
  for(g in unique(group)) {
    observations[which(sorted_group == g), ] = x_centered[which(group == g), ]
    group_indices[g, 1] = as.integer(min(which(sorted_group == g)) - 1)
    group_indices[g, 2] = as.integer(max(which(sorted_group == g)) - 1)
  }

  one = matrix(1, nrow = num_groups, ncol = num_groups)
  V = diag(num_groups) - one / num_groups
  projection = eigen(V)$vectors[, -num_groups, drop = FALSE]
  if(num_groups == 2) {
    projection = projection / sqrt(2)
  }

  new_bgm_spec(
    model_type = "compare",
    data = list(
      x                = observations,
      data_columnnames = data_columnnames,
      num_variables    = as.integer(num_variables),
      num_cases        = as.integer(nrow(observations)),
      num_categories   = as.integer(num_categories),
      group            = as.integer(group),
      num_groups       = as.integer(num_groups),
      group_indices    = group_indices,
      projection       = projection
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = ordinal_variable,
      is_continuous     = FALSE,
      baseline_category = as.integer(bc_final)
    ),
    missing = list(
      na_action     = na_action,
      na_impute     = na_impute,
      missing_index = missing_index
    ),
    prior = list(
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha,
      main_beta = main_beta,
      threshold_scale = threshold_scale,
      standardize = standardize,
      pairwise_scaling_factors = psf,
      difference_selection = dp$difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = dp$difference_prior,
      difference_scale = difference_scale,
      inclusion_probability_difference = dp$inclusion_probability_difference,
      beta_bernoulli_alpha = dp$beta_bernoulli_alpha,
      beta_bernoulli_beta = dp$beta_bernoulli_beta,
      beta_bernoulli_alpha_between = dp$beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = dp$beta_bernoulli_beta_between,
      dirichlet_alpha = dp$dirichlet_alpha,
      lambda = dp$lambda
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list(
      counts_per_category      = counts_per_category,
      blume_capel_stats        = blume_capel_stats,
      pairwise_stats           = pairwise_stats,
      main_effect_indices      = main_effect_indices,
      pairwise_effect_indices  = pairwise_effect_indices,
      interaction_index_matrix = interaction_index_matrix
    )
  )
}


# ==============================================================================
# sampler_sublist()  --- extract validated sampler list for new_bgm_spec()
# ==============================================================================
sampler_sublist = function(s) {
  list(
    update_method     = s$update_method,
    target_accept     = s$target_accept,
    iter              = as.integer(s$iter),
    warmup            = as.integer(s$warmup),
    chains            = as.integer(s$chains),
    cores             = as.integer(s$cores),
    nuts_max_depth    = as.integer(s$nuts_max_depth),
    learn_mass_matrix = s$learn_mass_matrix,
    seed              = as.integer(s$seed),
    progress_type     = as.integer(s$progress_type),
    progress_callback = s$progress_callback
  )
}


# ==============================================================================
# build_arguments()  --- convert spec -> arguments list for fit object
# ==============================================================================
#
# Produces the $arguments list stored in every bgms/bgmCompare fit object.
# Downstream code (extractor functions, simulate, predict, print, summary)
# reads this list to determine model properties.
# ==============================================================================
build_arguments = function(spec) {
  stopifnot(inherits(spec, "bgm_spec"))
  mt = spec$model_type

  if(mt == "ggm") {
    build_arguments_ggm(spec)
  } else if(mt == "omrf") {
    build_arguments_omrf(spec)
  } else if(mt == "mixed_mrf") {
    build_arguments_mixed_mrf(spec)
  } else {
    build_arguments_compare(spec)
  }
}


build_arguments_ggm = function(spec) {
  list(
    num_variables                = spec$data$num_variables,
    num_cases                    = spec$data$num_cases,
    na_impute                    = spec$missing$na_impute,
    variable_type                = spec$variables$variable_type,
    iter                         = spec$sampler$iter,
    warmup                       = spec$sampler$warmup,
    edge_selection               = spec$prior$edge_selection,
    edge_prior                   = spec$prior$edge_prior,
    inclusion_probability        = spec$prior$inclusion_probability,
    beta_bernoulli_alpha         = spec$prior$beta_bernoulli_alpha,
    beta_bernoulli_beta          = spec$prior$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = spec$prior$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between  = spec$prior$beta_bernoulli_beta_between,
    dirichlet_alpha              = spec$prior$dirichlet_alpha,
    lambda                       = spec$prior$lambda,
    na_action                    = spec$missing$na_action,
    version                      = packageVersion("bgms"),
    update_method                = spec$sampler$update_method,
    target_accept                = spec$sampler$target_accept,
    num_chains                   = spec$sampler$chains,
    data_columnnames             = spec$data$data_columnnames,
    no_variables                 = spec$data$num_variables,
    is_continuous                = TRUE,
    model_type                   = "ggm"
  )
}


build_arguments_omrf = function(spec) {
  # Legacy stores user-facing scalar (e.g. "ordinal") when all the same.
  vt = spec$variables$variable_type
  if(length(unique(vt)) == 1L) vt = unique(vt)

  list(
    num_variables                = spec$data$num_variables,
    num_cases                    = spec$data$num_cases,
    na_impute                    = spec$missing$na_impute,
    variable_type                = vt,
    iter                         = spec$sampler$iter,
    warmup                       = spec$sampler$warmup,
    pairwise_scale               = spec$prior$pairwise_scale,
    standardize                  = spec$prior$standardize,
    main_alpha                   = spec$prior$main_alpha,
    main_beta                    = spec$prior$main_beta,
    edge_selection               = spec$prior$edge_selection,
    edge_prior                   = spec$prior$edge_prior,
    inclusion_probability        = spec$prior$inclusion_probability,
    beta_bernoulli_alpha         = spec$prior$beta_bernoulli_alpha,
    beta_bernoulli_beta          = spec$prior$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = spec$prior$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between  = spec$prior$beta_bernoulli_beta_between,
    dirichlet_alpha              = spec$prior$dirichlet_alpha,
    lambda                       = spec$prior$lambda,
    na_action                    = spec$missing$na_action,
    version                      = packageVersion("bgms"),
    update_method                = spec$sampler$update_method,
    target_accept                = spec$sampler$target_accept,
    nuts_max_depth               = spec$sampler$nuts_max_depth,
    learn_mass_matrix            = spec$sampler$learn_mass_matrix,
    num_chains                   = spec$sampler$chains,
    num_categories               = spec$data$num_categories,
    data_columnnames             = spec$data$data_columnnames,
    baseline_category            = spec$variables$baseline_category,
    pairwise_scaling_factors     = spec$prior$pairwise_scaling_factors,
    no_variables                 = spec$data$num_variables,
    model_type                   = "omrf"
  )
}


build_arguments_mixed_mrf = function(spec) {
  list(
    num_variables                = spec$data$num_variables,
    num_discrete                 = spec$data$num_discrete,
    num_continuous               = spec$data$num_continuous,
    num_cases                    = spec$data$num_cases,
    variable_type                = spec$variables$variable_type,
    iter                         = spec$sampler$iter,
    warmup                       = spec$sampler$warmup,
    pairwise_scale               = spec$prior$pairwise_scale,
    standardize                  = spec$prior$standardize,
    main_alpha                   = spec$prior$main_alpha,
    main_beta                    = spec$prior$main_beta,
    edge_selection               = spec$prior$edge_selection,
    edge_prior                   = spec$prior$edge_prior,
    inclusion_probability        = spec$prior$inclusion_probability,
    beta_bernoulli_alpha         = spec$prior$beta_bernoulli_alpha,
    beta_bernoulli_beta          = spec$prior$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = spec$prior$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between  = spec$prior$beta_bernoulli_beta_between,
    dirichlet_alpha              = spec$prior$dirichlet_alpha,
    lambda                       = spec$prior$lambda,
    na_action                    = spec$missing$na_action,
    version                      = packageVersion("bgms"),
    update_method                = spec$sampler$update_method,
    target_accept                = spec$sampler$target_accept,
    nuts_max_depth               = spec$sampler$nuts_max_depth,
    num_chains                   = spec$sampler$chains,
    num_categories               = spec$data$num_categories,
    data_columnnames             = spec$data$data_columnnames,
    data_columnnames_discrete    = spec$data$data_columnnames_discrete,
    data_columnnames_continuous  = spec$data$data_columnnames_continuous,
    discrete_indices             = spec$data$discrete_indices,
    continuous_indices           = spec$data$continuous_indices,
    baseline_category            = spec$variables$baseline_category,
    is_ordinal                   = spec$variables$is_ordinal,
    no_variables                 = spec$data$num_variables,
    is_mixed                     = TRUE,
    model_type                   = "mixed_mrf"
  )
}


build_arguments_compare = function(spec) {
  list(
    num_variables                      = spec$data$num_variables,
    num_cases                          = spec$data$num_cases,
    iter                               = spec$sampler$iter,
    warmup                             = spec$sampler$warmup,
    pairwise_scale                     = spec$prior$pairwise_scale,
    difference_scale                   = spec$prior$difference_scale,
    standardize                        = spec$prior$standardize,
    difference_selection               = spec$prior$difference_selection,
    main_difference_selection          = spec$prior$main_difference_selection,
    difference_prior                   = spec$prior$difference_prior,
    difference_selection_alpha         = spec$prior$beta_bernoulli_alpha,
    difference_selection_beta          = spec$prior$beta_bernoulli_beta,
    difference_selection_alpha_between = spec$prior$beta_bernoulli_alpha_between,
    difference_selection_beta_between  = spec$prior$beta_bernoulli_beta_between,
    difference_dirichlet_alpha         = spec$prior$dirichlet_alpha,
    difference_lambda                  = spec$prior$lambda,
    inclusion_probability              = spec$prior$inclusion_probability_difference,
    version                            = packageVersion("bgms"),
    update_method                      = spec$sampler$update_method,
    target_accept                      = spec$sampler$target_accept,
    nuts_max_depth                     = spec$sampler$nuts_max_depth,
    learn_mass_matrix                  = spec$sampler$learn_mass_matrix,
    num_chains                         = spec$sampler$chains,
    num_groups                         = spec$data$num_groups,
    data_columnnames                   = spec$data$data_columnnames,
    projection                         = spec$data$projection,
    num_categories                     = spec$data$num_categories,
    is_ordinal_variable                = spec$variables$is_ordinal,
    group                              = sort(spec$data$group),
    pairwise_scaling_factors           = spec$prior$pairwise_scaling_factors,
    model_type                         = "compare"
  )
}


# ==============================================================================
# print.bgm_spec()  --- debugging summary
# ==============================================================================
#' @export
print.bgm_spec = function(x, ...) {
  s = x
  cat("bgm_spec object\n")
  cat("  model_type:", s$model_type, "\n")
  cat("  variables: ", s$data$num_variables, " (", s$data$num_cases, " cases)\n",
    sep = ""
  )
  cat(
    "  variable_type:",
    if(s$variables$is_continuous) {
      "continuous"
    } else {
      paste0(
        sum(s$variables$is_ordinal), " ordinal, ",
        sum(!s$variables$is_ordinal), " blume-capel"
      )
    },
    "\n"
  )
  cat("  sampler:", s$sampler$update_method,
    "(iter=", s$sampler$iter, ", warmup=", s$sampler$warmup,
    ", chains=", s$sampler$chains, ")\n",
    sep = ""
  )
  if(s$model_type %in% c("ggm", "omrf")) {
    cat(
      "  edge_selection:", s$prior$edge_selection,
      if(s$prior$edge_selection) paste0(" (", s$prior$edge_prior, ")"),
      "\n"
    )
  }
  if(s$model_type == "compare") {
    cat("  groups:", s$data$num_groups, "\n")
    cat(
      "  difference_selection:", s$prior$difference_selection,
      if(s$prior$difference_selection) paste0(" (", s$prior$difference_prior, ")"),
      "\n"
    )
  }
  cat(
    "  na_action:", s$missing$na_action,
    if(s$missing$na_impute) "(imputing)" else "(complete cases)", "\n"
  )
  invisible(s)
}
