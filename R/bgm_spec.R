# ==============================================================================
# bgm_spec: validated intermediate specification object
# ==============================================================================
#
# Central construction point for all bgm/bgmCompare models. Three layers:
#   bgm_spec()          — user-facing: validates inputs, assembles sub-lists
#   new_bgm_spec()      — low-level: type/presence assertions per field
#   validate_bgm_spec() — cross-field invariant checks
#
# The result is an S3 list of class "bgm_spec" consumed by run_sampler()
# and build_output().
# ==============================================================================


# ==============================================================================
# new_bgm_spec()  — low-level constructor
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

  # --- prior sub-list ---
  stopifnot(is.list(prior))
  if(model_type %in% c("omrf", "compare")) {
    stopifnot(is.numeric(prior$pairwise_scale), length(prior$pairwise_scale) == 1L)
    stopifnot(is.numeric(prior$main_alpha), length(prior$main_alpha) == 1L)
    stopifnot(is.numeric(prior$main_beta), length(prior$main_beta) == 1L)
    stopifnot(is.logical(prior$standardize), length(prior$standardize) == 1L)
    stopifnot(is.matrix(prior$pairwise_scaling_factors))
  }
  if(model_type == "mixed_mrf") {
    stopifnot(is.numeric(prior$pairwise_scale), length(prior$pairwise_scale) == 1L)
    stopifnot(is.numeric(prior$main_alpha), length(prior$main_alpha) == 1L)
    stopifnot(is.numeric(prior$main_beta), length(prior$main_beta) == 1L)
    stopifnot(is.logical(prior$standardize), length(prior$standardize) == 1L)
    stopifnot(is.character(prior$pseudolikelihood), length(prior$pseudolikelihood) == 1L)
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
  stopifnot(is.integer(sampler$hmc_num_leapfrogs), length(sampler$hmc_num_leapfrogs) == 1L)
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
# validate_bgm_spec()  — cross-field invariant checks
# ==============================================================================
validate_bgm_spec = function(spec) {
  mt = spec$model_type

  # GGM invariants
  if(mt == "ggm") {
    if(!isTRUE(spec$variables$is_continuous)) {
      stop("bgm_spec: model_type = 'ggm' requires is_continuous = TRUE.")
    }
    if(spec$sampler$update_method != "adaptive-metropolis") {
      stop("bgm_spec: model_type = 'ggm' requires update_method = 'adaptive-metropolis'.")
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
    if(spec$sampler$update_method != "adaptive-metropolis") {
      stop("bgm_spec: model_type = 'mixed_mrf' requires update_method = 'adaptive-metropolis'.")
    }
  }

  invisible(spec)
}


# ==============================================================================
# bgm_spec()  — user-facing constructor
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
                    # Priors (bgm / shared)
                    pairwise_scale = 2.5,
                    main_alpha = 0.5,
                    main_beta = 0.5,
                    standardize = FALSE,
                    edge_selection = TRUE,
                    edge_prior = c(
                      "Bernoulli", "Beta-Bernoulli",
                      "Stochastic-Block"
                    ),
                    inclusion_probability = 0.5,
                    beta_bernoulli_alpha = 1,
                    beta_bernoulli_beta = 1,
                    beta_bernoulli_alpha_between = 1,
                    beta_bernoulli_beta_between = 1,
                    dirichlet_alpha = 1,
                    lambda = 1,
                    # Priors (compare-specific)
                    difference_selection = TRUE,
                    main_difference_selection = FALSE,
                    difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                    difference_scale = 2.5,
                    difference_probability = 0.5,
                    # Sampler
                    update_method = c(
                      "nuts", "adaptive-metropolis",
                      "hamiltonian-mc"
                    ),
                    target_accept = NULL,
                    iter = 10000L,
                    warmup = 1000L,
                    hmc_num_leapfrogs = 100L,
                    nuts_max_depth = 10L,
                    learn_mass_matrix = TRUE,
                    chains = 4L,
                    cores = parallel::detectCores(),
                    seed = NULL,
                    display_progress = c("per-chain", "total", "none"),
                    verbose = TRUE,
                    pseudolikelihood = c("conditional", "marginal")) {
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

  # --- Sampler (needs is_continuous and edge_selection early) ------------------
  # Mixed MRF is MH-only (like GGM): force adaptive-metropolis via is_continuous
  sampler_is_continuous = is_continuous || is_mixed
  sampler = validate_sampler(
    update_method     = update_method,
    target_accept     = target_accept,
    iter              = iter,
    warmup            = warmup,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    nuts_max_depth    = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains            = chains,
    cores             = cores,
    seed              = seed,
    display_progress  = display_progress,
    is_continuous     = sampler_is_continuous,
    edge_selection    = if(model_type == "compare") FALSE else edge_selection,
    verbose           = verbose
  )

  # --- Build by model type ----------------------------------------------------
  if(model_type == "ggm") {
    spec = build_spec_ggm(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = as.integer(rep(0L, num_variables)),
      na_action = na_action, sampler = sampler,
      edge_selection = edge_selection,
      edge_prior = edge_prior,
      inclusion_probability = inclusion_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
    )
  } else if(model_type == "mixed_mrf") {
    pseudolikelihood = match.arg(pseudolikelihood)
    spec = build_spec_mixed_mrf(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      pairwise_scale = pairwise_scale, main_alpha = main_alpha,
      main_beta = main_beta, standardize = standardize,
      pseudolikelihood = pseudolikelihood,
      edge_selection = edge_selection,
      edge_prior = edge_prior,
      inclusion_probability = inclusion_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
    )
  } else if(model_type == "omrf") {
    spec = build_spec_omrf(
      x = x, data_columnnames = data_columnnames,
      num_variables = num_variables,
      variable_type = variable_type, is_ordinal = is_ordinal,
      is_continuous = is_continuous,
      baseline_category = baseline_category,
      na_action = na_action, sampler = sampler,
      pairwise_scale = pairwise_scale, main_alpha = main_alpha,
      main_beta = main_beta, standardize = standardize,
      edge_selection = edge_selection,
      edge_prior = edge_prior,
      inclusion_probability = inclusion_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
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
      pairwise_scale = pairwise_scale, main_alpha = main_alpha,
      main_beta = main_beta, standardize = standardize,
      difference_selection = difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = difference_prior,
      difference_scale = difference_scale,
      difference_probability = difference_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta
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
                          edge_selection, edge_prior,
                          inclusion_probability,
                          beta_bernoulli_alpha, beta_bernoulli_beta,
                          beta_bernoulli_alpha_between,
                          beta_bernoulli_beta_between,
                          dirichlet_alpha, lambda) {
  # Missing data
  md = validate_missing_data(
    x = x, na_action = na_action,
    is_continuous = TRUE
  )
  x = md$x

  # Center continuous data (GGM likelihood assumes zero mean)
  x = center_continuous_data(x)

  # Edge prior
  ep = validate_edge_prior(
    edge_selection = edge_selection, edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    num_variables = num_variables,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda
  )

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
      edge_selection = ep$edge_selection,
      edge_prior = ep$edge_prior,
      inclusion_probability = ep$inclusion_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list()
  )
}


build_spec_omrf = function(x, data_columnnames, num_variables,
                           variable_type, is_ordinal, is_continuous,
                           baseline_category,
                           na_action, sampler,
                           pairwise_scale, main_alpha, main_beta,
                           standardize,
                           edge_selection, edge_prior,
                           inclusion_probability,
                           beta_bernoulli_alpha, beta_bernoulli_beta,
                           beta_bernoulli_alpha_between,
                           beta_bernoulli_beta_between,
                           dirichlet_alpha, lambda) {
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

  # Edge prior
  ep = validate_edge_prior(
    edge_selection = edge_selection, edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    num_variables = num_variables,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda
  )

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
      pairwise_scale = pairwise_scale,
      main_alpha = main_alpha,
      main_beta = main_beta,
      standardize = standardize,
      pairwise_scaling_factors = psf,
      edge_selection = ep$edge_selection,
      edge_prior = ep$edge_prior,
      inclusion_probability = ep$inclusion_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
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
                                pairwise_scale, main_alpha, main_beta,
                                standardize, pseudolikelihood,
                                edge_selection, edge_prior,
                                inclusion_probability,
                                beta_bernoulli_alpha, beta_bernoulli_beta,
                                beta_bernoulli_alpha_between,
                                beta_bernoulli_beta_between,
                                dirichlet_alpha, lambda) {
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

  # Baseline category for discrete variables
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !identical(baseline_category, 0L),
    x = x_disc,
    variable_bool = is_ordinal_disc
  )

  # Missing data — not supported for mixed MRF (Phase H)
  if(na_action == "impute") {
    stop("Missing data imputation is not yet supported for mixed models.",
         call. = FALSE)
  }
  if(anyNA(x)) {
    stop("Missing data detected. Mixed models do not yet support missing data handling.",
         call. = FALSE)
  }

  # Ordinal recoding (reformat discrete data)
  ord = reformat_ordinal_data(
    x = x_disc, is_ordinal = is_ordinal_disc,
    baseline_category = bc
  )
  x_disc_recoded = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  # Edge prior (total variables = p + q)
  ep = validate_edge_prior(
    edge_selection = edge_selection, edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    num_variables = num_variables,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda
  )

  num_thresholds = sum(ifelse(is_ordinal_disc, num_categories, 2L))

  new_bgm_spec(
    model_type = "mixed_mrf",
    data = list(
      x_discrete       = x_disc_recoded,
      x_continuous     = x_cont,
      data_columnnames = data_columnnames,
      data_columnnames_discrete   = data_columnnames[disc_idx],
      data_columnnames_continuous = data_columnnames[cont_idx],
      num_variables    = as.integer(num_variables),
      num_discrete     = as.integer(p),
      num_continuous   = as.integer(q),
      num_cases        = as.integer(nrow(x)),
      num_categories   = as.integer(num_categories),
      discrete_indices   = disc_idx,
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
      na_action     = na_action,
      na_impute     = FALSE,
      missing_index = NULL
    ),
    prior = list(
      pairwise_scale = pairwise_scale,
      main_alpha = main_alpha,
      main_beta = main_beta,
      standardize = standardize,
      pseudolikelihood = pseudolikelihood,
      edge_selection = ep$edge_selection,
      edge_prior = ep$edge_prior,
      inclusion_probability = ep$inclusion_probability,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta,
      beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = beta_bernoulli_beta_between,
      dirichlet_alpha = dirichlet_alpha,
      lambda = lambda
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
                              pairwise_scale, main_alpha, main_beta,
                              standardize,
                              difference_selection, main_difference_selection,
                              difference_prior,
                              difference_scale, difference_probability,
                              beta_bernoulli_alpha, beta_bernoulli_beta) {
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
    beta_bernoulli_beta = beta_bernoulli_beta
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
      pairwise_scale = pairwise_scale,
      main_alpha = main_alpha,
      main_beta = main_beta,
      standardize = standardize,
      pairwise_scaling_factors = psf,
      difference_selection = dp$difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = dp$difference_prior,
      difference_scale = difference_scale,
      inclusion_probability_difference = dp$inclusion_probability_difference,
      beta_bernoulli_alpha = beta_bernoulli_alpha,
      beta_bernoulli_beta = beta_bernoulli_beta
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
# sampler_sublist()  — extract validated sampler list for new_bgm_spec()
# ==============================================================================
sampler_sublist = function(s) {
  list(
    update_method     = s$update_method,
    target_accept     = s$target_accept,
    iter              = as.integer(s$iter),
    warmup            = as.integer(s$warmup),
    chains            = as.integer(s$chains),
    cores             = as.integer(s$cores),
    hmc_num_leapfrogs = as.integer(s$hmc_num_leapfrogs),
    nuts_max_depth    = as.integer(s$nuts_max_depth),
    learn_mass_matrix = s$learn_mass_matrix,
    seed              = as.integer(s$seed),
    progress_type     = as.integer(s$progress_type)
  )
}


# ==============================================================================
# build_arguments()  — convert spec → arguments list for fit object
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
    is_continuous                = TRUE
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
    hmc_num_leapfrogs            = spec$sampler$hmc_num_leapfrogs,
    nuts_max_depth               = spec$sampler$nuts_max_depth,
    learn_mass_matrix            = spec$sampler$learn_mass_matrix,
    num_chains                   = spec$sampler$chains,
    num_categories               = spec$data$num_categories,
    data_columnnames             = spec$data$data_columnnames,
    baseline_category            = spec$variables$baseline_category,
    pairwise_scaling_factors     = spec$prior$pairwise_scaling_factors,
    no_variables                 = spec$data$num_variables
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
    pseudolikelihood             = spec$prior$pseudolikelihood,
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
    is_mixed                     = TRUE
  )
}


build_arguments_compare = function(spec) {
  list(
    num_variables                = spec$data$num_variables,
    num_cases                    = spec$data$num_cases,
    iter                         = spec$sampler$iter,
    warmup                       = spec$sampler$warmup,
    pairwise_scale               = spec$prior$pairwise_scale,
    difference_scale             = spec$prior$difference_scale,
    standardize                  = spec$prior$standardize,
    difference_selection         = spec$prior$difference_selection,
    main_difference_selection    = spec$prior$main_difference_selection,
    difference_prior             = spec$prior$difference_prior,
    difference_selection_alpha   = spec$prior$beta_bernoulli_alpha,
    difference_selection_beta    = spec$prior$beta_bernoulli_beta,
    inclusion_probability        = spec$prior$inclusion_probability_difference,
    version                      = packageVersion("bgms"),
    update_method                = spec$sampler$update_method,
    target_accept                = spec$sampler$target_accept,
    hmc_num_leapfrogs            = spec$sampler$hmc_num_leapfrogs,
    nuts_max_depth               = spec$sampler$nuts_max_depth,
    learn_mass_matrix            = spec$sampler$learn_mass_matrix,
    num_chains                   = spec$sampler$chains,
    num_groups                   = spec$data$num_groups,
    data_columnnames             = spec$data$data_columnnames,
    projection                   = spec$data$projection,
    num_categories               = spec$data$num_categories,
    is_ordinal_variable          = spec$variables$is_ordinal,
    group                        = sort(spec$data$group),
    pairwise_scaling_factors     = spec$prior$pairwise_scaling_factors
  )
}


# ==============================================================================
# print.bgm_spec()  — debugging summary
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
