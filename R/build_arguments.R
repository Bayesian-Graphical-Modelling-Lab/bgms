# build_arguments*() — convert a bgm_spec into the fit $arguments list
#
# Split out of bgm_spec.R (cleanup S4).


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

  arguments = if(mt == "ggm") {
    build_arguments_ggm(spec)
  } else if(mt == "omrf") {
    build_arguments_omrf(spec)
  } else if(mt == "mixed_mrf") {
    build_arguments_mixed_mrf(spec)
  } else {
    build_arguments_compare(spec)
  }
  c(arguments, build_arguments_sbm_provenance(spec))
}


# ==============================================================================
# build_arguments_sbm_provenance()  --- what a fit's SBM draws rest on
# ==============================================================================
#
# Recorded in every fit, NA without a Stochastic-Block prior on the
# indicators: the prior family on the number of components; whether the
# plain SBM update draws a new cluster's block probabilities from their
# conditional posterior (TRUE from bgms 0.2.1.0; NA under the joint
# precision-graph prior, which uses a separate update); and the convention of
# the stored posterior_num_blocks (shifted Poisson, B <= p, renormalized).
# ==============================================================================
build_arguments_sbm_provenance = function(spec) {
  sbm = if(spec$model_type == "compare") {
    identical(spec$prior$difference_prior, "Stochastic-Block")
  } else {
    identical(spec$prior$edge_prior, "Stochastic-Block")
  }
  joint = identical(spec$prior$precision_graph_prior, "joint")
  list(
    sbm_prior_family = if(sbm) "shifted Poisson" else NA_character_,
    sbm_sampler_fix = if(sbm && !joint) TRUE else NA,
    sbm_summary_convention = if(sbm) "shifted" else NA_character_
  )
}


build_arguments_ggm = function(spec) {
  list(
    num_variables = spec$data$num_variables,
    num_cases = spec$data$num_cases,
    na_impute = spec$missing$na_impute,
    variable_type = spec$variables$variable_type,
    iter = spec$sampler$iter,
    warmup = spec$sampler$warmup,
    edge_selection = spec$prior$edge_selection,
    edge_prior = spec$prior$edge_prior,
    precision_graph_prior = spec$prior$precision_graph_prior,
    inclusion_probability = spec$prior$inclusion_probability,
    beta_bernoulli_alpha = spec$prior$beta_bernoulli_alpha,
    beta_bernoulli_beta = spec$prior$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = spec$prior$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = spec$prior$beta_bernoulli_beta_between,
    dirichlet_alpha = spec$prior$dirichlet_alpha,
    lambda = spec$prior$lambda,
    na_action = spec$missing$na_action,
    version = packageVersion("bgms"),
    build_commit = bgms_build_commit(),
    update_method = spec$sampler$update_method,
    target_accept = spec$sampler$target_accept,
    num_chains = spec$sampler$chains,
    data_columnnames = spec$data$data_columnnames,
    no_variables = spec$data$num_variables,
    column_means = spec$data$column_means,
    is_continuous = TRUE,
    model_type = "ggm"
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
    build_commit                 = bgms_build_commit(),
    update_method                = spec$sampler$update_method,
    target_accept                = spec$sampler$target_accept,
    nuts_max_depth               = spec$sampler$nuts_max_depth,
    learn_mass_matrix            = spec$sampler$learn_mass_matrix,
    num_chains                   = spec$sampler$chains,
    num_categories               = spec$data$num_categories,
    category_levels              = spec$data$category_levels,
    blume_capel_shift            = spec$data$blume_capel_shift,
    data_columnnames             = spec$data$data_columnnames,
    baseline_category            = spec$variables$baseline_category,
    no_variables                 = spec$data$num_variables,
    model_type                   = "omrf"
  )
}


build_arguments_mixed_mrf = function(spec) {
  list(
    num_variables = spec$data$num_variables,
    num_discrete = spec$data$num_discrete,
    num_continuous = spec$data$num_continuous,
    num_cases = spec$data$num_cases,
    variable_type = spec$variables$variable_type,
    iter = spec$sampler$iter,
    warmup = spec$sampler$warmup,
    pairwise_scale = spec$prior$pairwise_scale,
    main_alpha = spec$prior$main_alpha,
    main_beta = spec$prior$main_beta,
    edge_selection = spec$prior$edge_selection,
    edge_prior = spec$prior$edge_prior,
    precision_graph_prior = spec$prior$precision_graph_prior,
    inclusion_probability = spec$prior$inclusion_probability,
    beta_bernoulli_alpha = spec$prior$beta_bernoulli_alpha,
    beta_bernoulli_beta = spec$prior$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = spec$prior$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = spec$prior$beta_bernoulli_beta_between,
    dirichlet_alpha = spec$prior$dirichlet_alpha,
    lambda = spec$prior$lambda,
    na_action = spec$missing$na_action,
    na_impute = spec$missing$na_impute,
    version = packageVersion("bgms"),
    build_commit = bgms_build_commit(),
    update_method = spec$sampler$update_method,
    target_accept = spec$sampler$target_accept,
    nuts_max_depth = spec$sampler$nuts_max_depth,
    num_chains = spec$sampler$chains,
    num_categories = spec$data$num_categories,
    category_levels = spec$data$category_levels,
    blume_capel_shift = spec$data$blume_capel_shift,
    data_columnnames = spec$data$data_columnnames,
    data_columnnames_discrete = spec$data$data_columnnames_discrete,
    data_columnnames_continuous = spec$data$data_columnnames_continuous,
    discrete_indices = spec$data$discrete_indices,
    continuous_indices = spec$data$continuous_indices,
    baseline_category = spec$variables$baseline_category,
    is_ordinal = spec$variables$is_ordinal,
    no_variables = spec$data$num_variables,
    is_mixed = TRUE,
    model_type = "mixed_mrf"
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
    build_commit                       = bgms_build_commit(),
    update_method                      = spec$sampler$update_method,
    target_accept                      = spec$sampler$target_accept,
    nuts_max_depth                     = spec$sampler$nuts_max_depth,
    learn_mass_matrix                  = spec$sampler$learn_mass_matrix,
    num_chains                         = spec$sampler$chains,
    num_groups                         = spec$data$num_groups,
    data_columnnames                   = spec$data$data_columnnames,
    projection                         = spec$data$projection,
    num_categories                     = spec$data$num_categories,
    category_levels                    = spec$data$category_levels,
    # Per-group counts on the final category codes: which group-by-category
    # cells are empty, and therefore which main-effect differences rest on the
    # prior rather than on data. See collapse_categories_across_groups().
    category_support                   = spec$data$category_support,
    blume_capel_shift                  = spec$data$blume_capel_shift,
    is_ordinal_variable                = spec$variables$is_ordinal,
    baseline_category                  = spec$variables$baseline_category,
    group                              = sort(spec$data$group),
    group_labels                       = spec$data$group_labels,
    model_type                         = "compare"
  )
}


# ==============================================================================
# bgms_build_commit()  --- git commit the installed bgms was built from
# ==============================================================================
#
# inst/BUILD_COMMIT holds "$Format:%H$", which git archive (and GitHub
# tarballs) replace with the commit hash through the export-subst attribute
# in .gitattributes. Installs through remotes/pak record RemoteSha instead.
# Returns NA when neither identifies the build (e.g. R CMD build run directly
# on a working tree, or devtools::load_all()).
# ==============================================================================
bgms_build_commit = function() {
  is_sha = function(x) {
    is.character(x) && length(x) == 1L && grepl("^[0-9a-f]{40}$", x)
  }
  path = system.file("BUILD_COMMIT", package = "bgms")
  if(nzchar(path)) {
    sha = trimws(readLines(path, n = 1L, warn = FALSE))
    if(is_sha(sha)) {
      return(sha)
    }
  }
  remote = utils::packageDescription("bgms")$RemoteSha
  if(is_sha(remote)) {
    return(remote)
  }
  NA_character_
}
