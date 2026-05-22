# ==============================================================================
# extract_prior_inclusion_probabilities() — parallel to the posterior extractor
# ==============================================================================
#
# Returns a symmetric p x p matrix of prior edge-inclusion probabilities for
# a bgms fit, in the same shape and orientation as
# extract_posterior_inclusion_probabilities().
#
# For most interaction priors the prior PIP is constant across edges,
# derived analytically from the edge prior:
#   - Bernoulli(p):              prior PIP = p
#   - Beta-Bernoulli(alpha, beta): prior PIP = alpha / (alpha + beta)
#
# For the joint-specification Graphical G-prior on GGM, the per-edge
# prior PIP is not analytic (the joint puts an implicit Z(Gamma)
# weighting on the marginal on Gamma) and must be estimated by running
# a short prior-only chain that samples (K, Gamma) under the same V_ij
# as the data fit. The estimate is cached on fit$cache so subsequent
# calls are free.
# ==============================================================================


# ------------------------------------------------------------------
# Accessors
# ------------------------------------------------------------------
get_fit_spec = function(fit) {
  if(inherits(fit, "S7_object")) {
    S7::prop(fit, ".bgm_spec")
  } else {
    .subset2(fit, ".bgm_spec")
  }
}


# ------------------------------------------------------------------
# run_prior_only_chain (internal)
# ------------------------------------------------------------------
# Runs a chain on the joint prior p(K, Gamma) using the fit's existing
# spec. Likelihood is muted via C++ `prior_only = TRUE`, which zeros
# n_ and suf_stat_ after V_ij has been cached. Returns the per-edge
# inclusion proportions in the same upper-triangle row-major order as
# posterior_mean_indicator.
# ------------------------------------------------------------------
run_prior_only_chain = function(spec, iter = 2000L, warmup = 500L,
                                seed_offset = 0L) {
  stopifnot(inherits(spec, "bgm_spec"))
  if(!identical(spec$model_type, "ggm")) {
    stop("Prior simulation is currently supported only for GGM fits ",
         "(continuous data). Got model_type = '", spec$model_type, "'.")
  }
  prior_spec = spec
  prior_spec$prior$prior_only = TRUE
  prior_spec$sampler$iter   = as.integer(iter)
  prior_spec$sampler$warmup = as.integer(warmup)
  prior_spec$sampler$chains = 1L
  prior_spec$sampler$cores  = 1L
  prior_spec$sampler$seed   = as.integer(spec$sampler$seed + seed_offset)
  prior_spec$sampler$progress_type = 0L
  prior_spec$sampler$display_progress = "none"
  prior_spec$sampler$update_method = "adaptive-metropolis"
  prior_spec$sampler$target_accept = 0.44

  raw = run_sampler(prior_spec)
  if(length(raw) == 0L) {
    stop("Prior simulation: empty chain output.")
  }
  ind = raw[[1L]]$indicator_samples
  if(is.null(ind)) {
    stop("Prior simulation: chain did not return indicator samples — ",
         "was edge_selection = TRUE on the fit?")
  }

  # ind is (p(p+1)/2) x n_iter with main-then-pair layout. Reduce to
  # off-diagonals (i < j) in row-major upper-triangle order.
  num_variables = spec$data$num_variables
  offdiag_idx = integer(num_variables * (num_variables - 1L) / 2L)
  pos = 0L; oi = 0L
  for(i in seq_len(num_variables)) {
    for(j in i:num_variables) {
      pos = pos + 1L
      if(i != j) {
        oi = oi + 1L
        offdiag_idx[oi] = pos
      }
    }
  }
  rowMeans(ind[offdiag_idx, , drop = FALSE])
}


# ------------------------------------------------------------------
# Analytic prior PIP (constant across edges)
# ------------------------------------------------------------------
# Reads the edge prior from the fit spec and returns the marginal
# P(gamma_ij = 1). Errors out on edge priors without a closed-form
# marginal (SBM and friends).
# ------------------------------------------------------------------
analytic_prior_pip = function(fit_spec) {
  p = fit_spec$prior
  switch(p$edge_prior,
    Bernoulli = mean(p$inclusion_probability),
    `Beta-Bernoulli` = p$beta_bernoulli_alpha /
      (p$beta_bernoulli_alpha + p$beta_bernoulli_beta),
    stop("extract_prior_inclusion_probabilities(): edge_prior = '",
         p$edge_prior, "' has no closed-form marginal. Currently ",
         "supported analytically: 'Bernoulli', 'Beta-Bernoulli'.")
  )
}


# ------------------------------------------------------------------
# Does this fit need the prior-only chain?
# ------------------------------------------------------------------
# Joint-spec Graphical G-prior is the only case where the marginal
# P(gamma_ij = 1) depends on V_ij and we can't read it off the edge
# prior. Everything else (Cauchy, Normal, beta-prime, ...) factorises
# and the prior PIP equals the edge-prior marginal.
# ------------------------------------------------------------------
needs_prior_simulation = function(fit_spec) {
  identical(fit_spec$prior$interaction_prior_type, "graphical_g")
}


# ------------------------------------------------------------------
# extract_prior_inclusion_probabilities() — main extractor
# ------------------------------------------------------------------

#' Extract Prior Inclusion Probabilities
#'
#' @description
#' Returns the prior inclusion probabilities for the edges of a model
#' fitted with [bgm()] (currently GGM only). The returned matrix has
#' the same shape and orientation as
#' [extract_posterior_inclusion_probabilities()] so the two can be
#' combined element-wise.
#'
#' For most interaction priors the prior PIP is constant across edges
#' and read analytically from the edge prior:
#' \itemize{
#'   \item Bernoulli(p): prior PIP = p
#'   \item Beta-Bernoulli(alpha, beta): prior PIP = alpha / (alpha + beta)
#' }
#'
#' For the joint-specification Graphical G-prior the marginal
#' P(gamma_ij = 1) depends on the data-defined V_ij and is not
#' analytic; it is estimated by running a short prior-only chain that
#' samples (K, Gamma) under the same V_ij as the data fit. The estimate
#' is cached on the fit so subsequent calls are free; pass
#' `recompute = TRUE` to re-run.
#'
#' @param bgms_object A fitted model object of class `bgms` from
#'   [bgm()] run with `edge_selection = TRUE`.
#' @param iter Integer. Number of post-warmup iterations for the
#'   prior-only chain. Used only when simulation is required (joint-spec
#'   Graphical G-prior). Default `2000`.
#' @param warmup Integer. Warmup length for the prior-only chain.
#'   Default `500`.
#' @param recompute Logical. If `TRUE`, recompute even if a cached
#'   result is present. Default `FALSE`.
#'
#' @return A symmetric p x p matrix of prior inclusion probabilities,
#'   with variable names as row and column names. Diagonal entries are
#'   `NA`. For non-Graphical-G priors all off-diagonal entries share
#'   the same value (the edge-prior marginal).
#'
#' @seealso [extract_posterior_inclusion_probabilities()], [bgm()]
#' @family extractors
#' @export
extract_prior_inclusion_probabilities = function(bgms_object,
                                                 iter = 2000L,
                                                 warmup = 500L,
                                                 recompute = FALSE) {
  UseMethod("extract_prior_inclusion_probabilities")
}


#' @inheritParams extract_prior_inclusion_probabilities
#' @exportS3Method
#' @noRd
extract_prior_inclusion_probabilities.bgms = function(bgms_object,
                                                      iter = 2000L,
                                                      warmup = 500L,
                                                      recompute = FALSE) {
  fit_spec = get_fit_spec(bgms_object)
  if(is.null(fit_spec)) {
    stop("extract_prior_inclusion_probabilities(): fit has no ",
         "embedded .bgm_spec. Re-fit with the current bgms version.")
  }
  if(!identical(fit_spec$model_type, "ggm")) {
    stop("extract_prior_inclusion_probabilities() currently supports ",
         "GGM fits only (continuous data). Got model_type = '",
         fit_spec$model_type, "'.")
  }
  if(!isTRUE(fit_spec$prior$edge_selection)) {
    stop("To extract prior inclusion probabilities, run bgm() with ",
         "edge_selection = TRUE.")
  }

  num_vars = fit_spec$data$num_variables
  data_columnnames = fit_spec$data$data_columnnames
  n_edges = num_vars * (num_vars - 1L) / 2L

  if(needs_prior_simulation(fit_spec)) {
    cache = get_fit_cache(bgms_object)
    edge_pips = if(!recompute && !is.null(cache) &&
                   !is.null(cache$prior_inclusion_pips) &&
                   length(cache$prior_inclusion_pips) == n_edges) {
      cache$prior_inclusion_pips
    } else {
      p = run_prior_only_chain(fit_spec, iter = iter, warmup = warmup)
      if(!is.null(cache)) cache$prior_inclusion_pips = p
      p
    }
  } else {
    edge_pips = rep(analytic_prior_pip(fit_spec), n_edges)
  }

  # Build symmetric matrix in the same orientation as
  # extract_posterior_inclusion_probabilities (lower.tri filled, then
  # mirrored). NA on the diagonal.
  pip_matrix = matrix(NA_real_, num_vars, num_vars)
  pip_matrix[lower.tri(pip_matrix)] = edge_pips
  pip_matrix[upper.tri(pip_matrix)] = t(pip_matrix)[upper.tri(pip_matrix)]
  colnames(pip_matrix) = data_columnnames
  rownames(pip_matrix) = data_columnnames
  pip_matrix
}
