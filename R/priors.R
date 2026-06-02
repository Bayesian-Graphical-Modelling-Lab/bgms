# ==============================================================================
# Prior specification classes for bgms
# ==============================================================================
#
# BAS-inspired prior class system. Each prior family has a named constructor
# that returns an S3 object of the appropriate class. These objects are passed
# to bgm() and bgmCompare() instead of loose numeric parameters.
#
# Three prior roles:
#   - bgms_parameter_prior   : prior on real-valued model parameters
#                               (interactions, thresholds, means)
#   - bgms_indicator_prior   : prior on inclusion indicators (edge selection,
#                               difference selection)
#
# Interaction vs threshold is determined by the bgm()/bgmCompare() argument
# slot (interaction_prior / threshold_prior) and the prior's `family` field,
# not by the prior's class.
# ==============================================================================


# ==============================================================================
# Parameter priors (interactions, thresholds, means)
# ==============================================================================

#' Cauchy Prior for Model Parameters
#'
#' @description
#' Specifies a Cauchy(0, scale) prior on model parameters.
#' This is the default prior for pairwise interactions in \code{\link{bgm}}
#' and produces heavy-tailed shrinkage toward zero.
#'
#' @param scale Positive numeric. Scale (half-width at half-maximum) of the
#'   Cauchy distribution. Default: \code{1}.
#'
#' @return An object of class \code{"bgms_parameter_prior"} with
#'   \code{family = "cauchy"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{normal_prior}}, \code{\link{beta_prime_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' cauchy_prior()
#' cauchy_prior(scale = 2.5)
#'
#' @export
cauchy_prior = function(scale = 1) {
  if(!is.numeric(scale) || length(scale) != 1L || is.na(scale)) {
    stop("'scale' must be a single positive number.")
  }
  if(scale <= 0) {
    stop("'scale' must be positive.")
  }
  if(!is.finite(scale)) {
    stop("'scale' must be finite.")
  }

  structure(
    list(
      family = "cauchy",
      hyper.parameters = list(scale = scale)
    ),
    class = "bgms_parameter_prior"
  )
}


#' Normal Prior for Model Parameters
#'
#' @description
#' Specifies a Normal(0, scale) prior on model parameters.
#' Produces lighter-tailed shrinkage than the Cauchy prior and is better
#' suited for simulation-based calibration (SBC) studies. Can be used for
#' interactions, thresholds, or continuous means.
#'
#' @param scale Positive numeric. Standard deviation of the normal
#'   distribution. Default: \code{1}.
#'
#' @return An object of class \code{"bgms_parameter_prior"} with
#'   \code{family = "normal"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{cauchy_prior}}, \code{\link{beta_prime_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' normal_prior()
#' normal_prior(scale = 0.5)
#'
#' @export
normal_prior = function(scale = 1) {
  if(!is.numeric(scale) || length(scale) != 1L || is.na(scale)) {
    stop("'scale' must be a single positive number.")
  }
  if(scale <= 0) {
    stop("'scale' must be positive.")
  }
  if(!is.finite(scale)) {
    stop("'scale' must be finite.")
  }

  structure(
    list(
      family = "normal",
      hyper.parameters = list(scale = scale)
    ),
    class = "bgms_parameter_prior"
  )
}


#' Beta-Prime Prior for Model Parameters
#'
#' @description
#' Specifies a beta-prime prior on model parameters.
#' The parameterization follows the logistic transformation:
#' \eqn{\sigma(\mu) \sim \textrm{Beta}(\alpha, \beta)}{sigma(mu) ~ Beta(alpha, beta)},
#' so \eqn{\mu = \textrm{logit}(Y)}{mu = logit(Y)} where
#' \eqn{Y \sim \textrm{Beta}(\alpha, \beta)}{Y ~ Beta(alpha, beta)}.
#'
#' @param alpha Positive numeric. First shape parameter. Default: \code{0.5}.
#' @param beta Positive numeric. Second shape parameter. Default: \code{0.5}.
#'
#' @return An object of class \code{"bgms_parameter_prior"} with
#'   \code{family = "beta-prime"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{cauchy_prior}}, \code{\link{normal_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' beta_prime_prior()
#' beta_prime_prior(alpha = 1, beta = 1)
#'
#' @export
beta_prime_prior = function(alpha = 0.5, beta = 0.5) {
  if(!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) {
    stop("'alpha' must be a single positive number.")
  }
  if(!is.numeric(beta) || length(beta) != 1L || is.na(beta)) {
    stop("'beta' must be a single positive number.")
  }
  if(alpha <= 0 || beta <= 0) {
    stop("'alpha' and 'beta' must be positive.")
  }
  if(!is.finite(alpha) || !is.finite(beta)) {
    stop("'alpha' and 'beta' must be finite.")
  }

  structure(
    list(
      family = "beta-prime",
      hyper.parameters = list(alpha = alpha, beta = beta)
    ),
    class = "bgms_parameter_prior"
  )
}


# ==============================================================================
# Scale priors (positive parameters: precision diagonal)
# ==============================================================================

#' Gamma Prior for Scale Parameters
#'
#' @description
#' Specifies a Gamma(shape, rate) prior for positive scale parameters such as
#' the diagonal elements of the precision matrix. The default
#' \code{gamma_prior(1, 1)} corresponds to an Exponential(1) distribution.
#'
#' @param shape Positive numeric. Shape parameter of the Gamma distribution.
#'   Default: \code{1}.
#' @param rate Positive numeric. Rate parameter of the Gamma distribution.
#'   Default: \code{1}.
#'
#' @return An object of class \code{"bgms_scale_prior"} with
#'   \code{family = "gamma"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{exponential_prior}}, \code{\link{bgm}}
#'
#' @examples
#' gamma_prior()
#' gamma_prior(shape = 2, rate = 0.5)
#'
#' @export
gamma_prior = function(shape = 1, rate = 1) {
  if(!is.numeric(shape) || length(shape) != 1L || is.na(shape)) {
    stop("'shape' must be a single positive number.")
  }
  if(!is.numeric(rate) || length(rate) != 1L || is.na(rate)) {
    stop("'rate' must be a single positive number.")
  }
  if(shape <= 0 || rate <= 0) {
    stop("'shape' and 'rate' must be positive.")
  }
  if(!is.finite(shape) || !is.finite(rate)) {
    stop("'shape' and 'rate' must be finite.")
  }

  structure(
    list(
      family = "gamma",
      hyper.parameters = list(shape = shape, rate = rate)
    ),
    class = "bgms_scale_prior"
  )
}


#' Exponential Prior for Scale Parameters
#'
#' @description
#' Specifies an Exponential(rate) prior for positive scale parameters.
#' This is a convenience function equivalent to
#' \code{gamma_prior(shape = 1, rate = rate)}.
#'
#' @param rate Positive numeric. Rate parameter of the Exponential
#'   distribution. Default: \code{1}.
#'
#' @return An object of class \code{"bgms_scale_prior"} with
#'   \code{family = "exponential"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{gamma_prior}}, \code{\link{bgm}}
#'
#' @examples
#' exponential_prior()
#' exponential_prior(rate = 2)
#'
#' @export
exponential_prior = function(rate = 1) {
  if(!is.numeric(rate) || length(rate) != 1L || is.na(rate)) {
    stop("'rate' must be a single positive number.")
  }
  if(rate <= 0) {
    stop("'rate' must be positive.")
  }
  if(!is.finite(rate)) {
    stop("'rate' must be finite.")
  }

  structure(
    list(
      family = "exponential",
      hyper.parameters = list(rate = rate)
    ),
    class = "bgms_scale_prior"
  )
}


# ==============================================================================
# Indicator priors (selection: edges, differences)
# ==============================================================================

#' Bernoulli Prior for Inclusion Indicators
#'
#' @description
#' Specifies a Bernoulli prior for inclusion indicators with a fixed
#' inclusion probability. Used for edge selection in \code{\link{bgm}} and
#' difference selection in \code{\link{bgmCompare}}.
#'
#' @param inclusion_probability Numeric scalar or symmetric matrix. Prior
#'   probability of each edge being included. A scalar applies to all edges;
#'   a matrix allows edge-specific probabilities. Must be in (0, 1).
#'   Default: \code{0.5}.
#'
#' @return An object of class \code{"bgms_indicator_prior"} with
#'   \code{family = "Bernoulli"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{beta_bernoulli_prior}}, \code{\link{sbm_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' bernoulli_prior()
#' bernoulli_prior(inclusion_probability = 0.25)
#'
#' @export
bernoulli_prior = function(inclusion_probability = 0.5) {
  structure(
    list(
      family = "Bernoulli",
      hyper.parameters = list(
        inclusion_probability = inclusion_probability
      )
    ),
    class = "bgms_indicator_prior"
  )
}


#' Beta-Bernoulli Prior for Inclusion Indicators
#'
#' @description
#' Specifies a Beta-Bernoulli prior for inclusion indicators. The inclusion
#' probability is drawn from a \eqn{\textrm{Beta}(\alpha, \beta)}{Beta(alpha, beta)}
#' distribution and shared across all edges.
#'
#' @param alpha Positive numeric. First shape parameter of the Beta
#'   distribution. Default: \code{1}.
#' @param beta Positive numeric. Second shape parameter of the Beta
#'   distribution. Default: \code{1}.
#'
#' @return An object of class \code{"bgms_indicator_prior"} with
#'   \code{family = "Beta-Bernoulli"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{bernoulli_prior}}, \code{\link{sbm_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' beta_bernoulli_prior()
#' beta_bernoulli_prior(alpha = 2, beta = 5)
#'
#' @export
beta_bernoulli_prior = function(alpha = 1, beta = 1) {
  if(!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) {
    stop("'alpha' must be a single positive number.")
  }
  if(!is.numeric(beta) || length(beta) != 1L || is.na(beta)) {
    stop("'beta' must be a single positive number.")
  }
  if(alpha <= 0 || beta <= 0) {
    stop("'alpha' and 'beta' must be positive.")
  }
  if(!is.finite(alpha) || !is.finite(beta)) {
    stop("'alpha' and 'beta' must be finite.")
  }

  structure(
    list(
      family = "Beta-Bernoulli",
      hyper.parameters = list(alpha = alpha, beta = beta)
    ),
    class = "bgms_indicator_prior"
  )
}


#' Stochastic Block Model Prior for Inclusion Indicators
#'
#' @description
#' Specifies a Stochastic Block Model (SBM) prior for inclusion indicators.
#' Variables are assigned to latent clusters, with separate Beta priors
#' on within-cluster and between-cluster inclusion probabilities.
#'
#' @param alpha Positive numeric. First shape parameter of the Beta
#'   distribution for within-cluster edges. Default: \code{1}.
#' @param beta Positive numeric. Second shape parameter of the Beta
#'   distribution for within-cluster edges. Default: \code{1}.
#' @param alpha_between Positive numeric. First shape parameter of the Beta
#'   distribution for between-cluster edges. Default: \code{1}.
#' @param beta_between Positive numeric. Second shape parameter of the Beta
#'   distribution for between-cluster edges. Default: \code{1}.
#' @param dirichlet_alpha Positive numeric. Concentration parameter of the
#'   Dirichlet prior on cluster assignments. Default: \code{1}.
#' @param lambda Positive numeric. Rate parameter of the zero-truncated
#'   Poisson prior on the number of clusters. Default: \code{1}.
#'
#' @return An object of class \code{"bgms_indicator_prior"} with
#'   \code{family = "Stochastic-Block"}.
#'
#' @family prior-constructors
#' @seealso \code{\link{bernoulli_prior}}, \code{\link{beta_bernoulli_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' sbm_prior()
#' sbm_prior(alpha = 2, beta = 1, alpha_between = 1, beta_between = 5)
#'
#' @export
sbm_prior = function(alpha = 1, beta = 1,
                     alpha_between = 1, beta_between = 1,
                     dirichlet_alpha = 1, lambda = 1) {
  params = list(
    alpha = alpha, beta = beta,
    alpha_between = alpha_between, beta_between = beta_between,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda
  )

  for(nm in names(params)) {
    val = params[[nm]]
    if(!is.numeric(val) || length(val) != 1L || is.na(val)) {
      stop(sprintf("'%s' must be a single positive number.", nm))
    }
    if(val <= 0) {
      stop(sprintf("'%s' must be positive.", nm))
    }
    if(!is.finite(val)) {
      stop(sprintf("'%s' must be finite.", nm))
    }
  }

  structure(
    list(
      family = "Stochastic-Block",
      hyper.parameters = params
    ),
    class = "bgms_indicator_prior"
  )
}


# ==============================================================================
# Print methods
# ==============================================================================

#' @export
print.bgms_parameter_prior = function(x, ...) {
  hp = x$hyper.parameters
  switch(x$family,
    "cauchy" = cat(sprintf("Parameter prior: Cauchy(0, %.4g)\n", hp$scale)),
    "normal" = cat(sprintf("Parameter prior: Normal(0, %.4g)\n", hp$scale)),
    "beta-prime" = cat(sprintf(
      "Parameter prior: Beta-prime(alpha = %.4g, beta = %.4g)\n",
      hp$alpha, hp$beta
    )),
    cat(sprintf("Parameter prior: %s\n", x$family))
  )
  invisible(x)
}

#' @export
print.bgms_scale_prior = function(x, ...) {
  hp = x$hyper.parameters
  switch(x$family,
    "gamma" = cat(sprintf(
      "Scale prior: Gamma(shape = %.4g, rate = %.4g)\n",
      hp$shape, hp$rate
    )),
    "exponential" = cat(sprintf("Scale prior: Exponential(rate = %.4g)\n", hp$rate)),
    cat(sprintf("Scale prior: %s\n", x$family))
  )
  invisible(x)
}

#' @export
print.bgms_indicator_prior = function(x, ...) {
  hp = x$hyper.parameters
  switch(x$family,
    "Bernoulli" = {
      ip = hp$inclusion_probability
      if(is.matrix(ip)) {
        cat("Edge prior: Bernoulli (variable-specific inclusion probabilities)\n")
      } else {
        cat(sprintf("Edge prior: Bernoulli(%.4g)\n", ip))
      }
    },
    "Beta-Bernoulli" = {
      cat(sprintf(
        "Edge prior: Beta-Bernoulli(alpha = %.4g, beta = %.4g)\n",
        hp$alpha, hp$beta
      ))
    },
    "Stochastic-Block" = {
      cat(sprintf(
        paste0(
          "Edge prior: Stochastic-Block\n",
          "  Within:    Beta(%.4g, %.4g)\n",
          "  Between:   Beta(%.4g, %.4g)\n",
          "  Dirichlet: %.4g, Lambda: %.4g\n"
        ),
        hp$alpha, hp$beta,
        hp$alpha_between, hp$beta_between,
        hp$dirichlet_alpha, hp$lambda
      ))
    },
    cat(sprintf("Edge prior: %s\n", x$family))
  )
  invisible(x)
}


# ==============================================================================
# Internal: extract prior parameters for spec / C++ interface
# ==============================================================================

#' Unpack a parameter prior into the flat parameters used by bgm_spec
#'
#' @param prior A \code{bgms_parameter_prior} object.
#'
#' @return A list with \code{prior_type} (character), \code{scale} (numeric),
#'   \code{alpha} (numeric), and \code{beta} (numeric). Unused hyperparameters
#'   are set to \code{NA_real_}.
#'
#' @keywords internal
unpack_parameter_prior = function(prior) {
  if(!inherits(prior, "bgms_parameter_prior")) {
    stop(
      "Prior must be a bgms_parameter_prior object.",
      " Use cauchy_prior(), normal_prior(), or beta_prime_prior()."
    )
  }
  hp = prior$hyper.parameters
  list(
    prior_type = prior$family,
    scale = hp$scale %||% NA_real_,
    alpha = hp$alpha %||% NA_real_,
    beta = hp$beta %||% NA_real_
  )
}


#' Unpack a scale prior into the flat parameters used by bgm_spec
#'
#' @param prior A \code{bgms_scale_prior} object.
#'
#' @return A list with \code{scale_prior_type} (character),
#'   \code{scale_shape} (numeric), and \code{scale_rate} (numeric).
#'
#' @keywords internal
unpack_scale_prior = function(prior) {
  if(!inherits(prior, "bgms_scale_prior")) {
    stop(
      "Prior must be a bgms_scale_prior object.",
      " Use gamma_prior() or exponential_prior()."
    )
  }
  hp = prior$hyper.parameters
  # Exponential is Gamma(1, rate)
  shape = if(prior$family == "exponential") 1 else hp$shape
  rate = hp$rate
  list(
    scale_prior_type = prior$family,
    scale_shape = shape,
    scale_rate = rate
  )
}


#' Unpack an interaction prior into the flat parameters used by bgm_spec
#'
#' @param prior A \code{bgms_parameter_prior} object.
#'
#' @return A list with \code{interaction_prior_type} (character) and
#'   \code{pairwise_scale} (numeric).
#'
#' @keywords internal
unpack_interaction_prior = function(prior) {
  pp = unpack_parameter_prior(prior)
  list(
    interaction_prior_type = pp$prior_type,
    pairwise_scale = pp$scale,
    interaction_alpha = pp$alpha,
    interaction_beta = pp$beta
  )
}


#' Unpack a threshold prior into the flat parameters used by bgm_spec
#'
#' @param prior A \code{bgms_parameter_prior} object.
#'
#' @return A list with \code{threshold_prior_type} (character),
#'   \code{main_alpha}, \code{main_beta}, and \code{threshold_scale}.
#'
#' @keywords internal
unpack_threshold_prior = function(prior) {
  pp = unpack_parameter_prior(prior)
  list(
    threshold_prior_type = pp$prior_type,
    main_alpha = pp$alpha,
    main_beta = pp$beta,
    threshold_scale = pp$scale
  )
}


#' Unpack an edge prior into the flat parameters used by bgm_spec
#'
#' @param prior A \code{bgms_indicator_prior} object.
#' @param num_variables Integer. Number of variables (for inclusion matrix).
#'
#' @return A list matching the fields expected by \code{validate_edge_prior}
#'   output.
#'
#' @keywords internal
unpack_indicator_prior = function(prior, num_variables) {
  if(!inherits(prior, "bgms_indicator_prior")) {
    stop(
      "'edge_prior' must be a bgms_indicator_prior object.",
      " Use bernoulli_prior(), beta_bernoulli_prior(), or sbm_prior()."
    )
  }

  hp = prior$hyper.parameters

  switch(prior$family,
    "Bernoulli" = {
      ip = validate_bernoulli_inclusion(
        probability = hp$inclusion_probability,
        num_variables = num_variables,
        include_diagonal = FALSE,
        context = ""
      )
      list(
        edge_selection = TRUE,
        edge_prior = "Bernoulli",
        inclusion_probability = ip,
        beta_bernoulli_alpha = 1,
        beta_bernoulli_beta = 1,
        beta_bernoulli_alpha_between = 1,
        beta_bernoulli_beta_between = 1,
        dirichlet_alpha = 1,
        lambda = 1
      )
    },
    "Beta-Bernoulli" = {
      list(
        edge_selection = TRUE,
        edge_prior = "Beta-Bernoulli",
        inclusion_probability = matrix(0.5, num_variables, num_variables),
        beta_bernoulli_alpha = hp$alpha,
        beta_bernoulli_beta = hp$beta,
        beta_bernoulli_alpha_between = 1,
        beta_bernoulli_beta_between = 1,
        dirichlet_alpha = 1,
        lambda = 1
      )
    },
    "Stochastic-Block" = {
      list(
        edge_selection = TRUE,
        edge_prior = "Stochastic-Block",
        inclusion_probability = matrix(0.5, num_variables, num_variables),
        beta_bernoulli_alpha = hp$alpha,
        beta_bernoulli_beta = hp$beta,
        beta_bernoulli_alpha_between = hp$alpha_between,
        beta_bernoulli_beta_between = hp$beta_between,
        dirichlet_alpha = hp$dirichlet_alpha,
        lambda = hp$lambda
      )
    }
  )
}
