check_positive_integer = function(value, name) {
  if(!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value <= 0) {
    stop(sprintf("Parameter `%s` must be a positive integer. Got: %s", name, value))
  }
}

# Helper function for validating non-negative integers
check_non_negative_integer = function(value, name) {
  if(!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value < 0) {
    stop(sprintf("Parameter `%s` must be a non-negative integer. Got: %s", name, value))
  }
}

# Helper function for validating logical inputs
check_logical = function(value, name) {

  value = as.logical(value)
  if(is.na(value)) {
    stop(sprintf("Parameter `%s` must be TRUE or FALSE. Got: %s", name, value))
  }
  return(value)
}

# Helper function for validating optional seed parameter
# Returns a valid integer seed (generates one if NULL)
check_seed = function(seed) {
  if(is.null(seed)) {
    return(sample.int(.Machine$integer.max, 1L))
  }
  if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
    stop("Argument 'seed' must be a single non-negative integer.")
  }
  as.integer(seed)
}

check_model = function(x,
                       variable_type,
                       baseline_category,
                       pairwise_scale = 2.5,
                       main_alpha = 0.5,
                       main_beta = 0.5,
                       edge_selection = TRUE,
                       edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
                       inclusion_probability = 0.5,
                       beta_bernoulli_alpha = 1,
                       beta_bernoulli_beta = 1,
                       beta_bernoulli_alpha_between = 1,
                       beta_bernoulli_beta_between = 1,
                       dirichlet_alpha = dirichlet_alpha,
                       lambda = lambda) {
  # Check variable type input ---------------------------------------------------
  vt <- validate_variable_types(
    variable_type    = variable_type,
    num_variables    = ncol(x),
    allow_continuous = TRUE,
    caller           = "bgm"
  )
  variable_bool <- vt$variable_bool
  is_continuous <- vt$is_continuous

  # Check Blume-Capel variable input --------------------------------------------
  baseline_category <- validate_baseline_category(
    baseline_category          = baseline_category,
    baseline_category_provided = hasArg("baseline_category"),
    x                          = x,
    variable_bool              = variable_bool
  )

  # Check prior set-up for the interaction parameters ---------------------------
  if(pairwise_scale <= 0 || is.na(pairwise_scale) || is.infinite(pairwise_scale)) {
    stop("The scale of the Cauchy prior needs to be positive.")
  }

  # Check prior set-up for the threshold parameters -----------------------------
  if(main_alpha <= 0 | !is.finite(main_alpha)) {
    stop("Parameter main_alpha needs to be positive.")
  }
  if(main_beta <= 0 | !is.finite(main_beta)) {
    stop("Parameter main_beta needs to be positive.")
  }

  # Check set-up for the Bayesian edge selection model --------------------------
  ep <- validate_edge_prior(
    edge_selection               = edge_selection,
    edge_prior                   = edge_prior,
    inclusion_probability        = inclusion_probability,
    num_variables                = ncol(x),
    beta_bernoulli_alpha         = beta_bernoulli_alpha,
    beta_bernoulli_beta          = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between  = beta_bernoulli_beta_between,
    dirichlet_alpha              = dirichlet_alpha,
    lambda                       = lambda
  )
  edge_selection <- ep$edge_selection
  edge_prior     <- ep$edge_prior
  theta          <- ep$inclusion_probability

  return(list(
    variable_bool = variable_bool,
    baseline_category = baseline_category,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = theta,
    is_continuous = is_continuous
  ))
}


check_compare_model = function(
  x,
  y,
  group_indicator,
  difference_selection,
  variable_type,
  baseline_category,
  difference_scale = 2.5,
  difference_prior = c("Bernoulli", "Beta-Bernoulli"),
  difference_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  pairwise_scale = 2.5,
  main_alpha = 0.5,
  main_beta = 0.5
) {
  if(!is.null(group_indicator)) {
    unique_g = unique(group_indicator)
    if(length(unique_g) == 0) {
      stop(paste0(
        "The bgmCompare function expects at least two groups, but the input group_indicator contains\n",
        "no group value."
      ))
    }
    if(length(unique_g) == 1) {
      stop(paste0(
        "The bgmCompare function expects at least two groups, but the input group_indicator contains\n",
        "only one group value."
      ))
    }
    if(length(unique_g) == length(group_indicator)) {
      stop("The input group_indicator contains only unique group values.")
    }

    group = group_indicator
    for(u in unique_g) {
      group[group_indicator == u] = which(unique_g == u)
    }
    tab = tabulate(group)

    if(any(tab < 2)) {
      stop("One or more groups only had one member in the input group_indicator.")
    }
  } else {
    group = c(rep.int(1, times = nrow(x)), rep.int(2, times = nrow(y)))
    x = rbind(x, y)
  }

  # Check variable type input ---------------------------------------------------
  vt <- validate_variable_types(
    variable_type    = variable_type,
    num_variables    = ncol(x),
    allow_continuous = FALSE,
    caller           = "bgmCompare"
  )
  variable_bool <- vt$variable_bool

  # Check Blume-Capel variable input --------------------------------------------
  baseline_category <- validate_baseline_category(
    baseline_category          = baseline_category,
    baseline_category_provided = hasArg("baseline_category"),
    x                          = x,
    variable_bool              = variable_bool
  )

  # Check prior set-up for the interaction parameters ---------------------------
  if(pairwise_scale <= 0 || is.na(pairwise_scale) || is.infinite(pairwise_scale)) {
    stop("The scale of the Cauchy prior for the interactions needs to be positive.")
  }

  # Check prior set-up for the interaction differences --------------------------
  if(difference_scale <= 0 || is.na(difference_scale) || is.infinite(difference_scale)) {
    stop("The scale of the Cauchy prior for the differences needs to be positive.")
  }

  # Check prior set-up for the threshold parameters -----------------------------
  if(main_alpha <= 0 | !is.finite(main_alpha)) {
    stop("Parameter main_alpha needs to be positive.")
  }
  if(main_beta <= 0 | !is.finite(main_beta)) {
    stop("Parameter main_beta needs to be positive.")
  }

  # Check set-up for the Bayesian difference selection model --------------------
  difference_selection = as.logical(difference_selection)
  if(is.na(difference_selection)) {
    stop("The parameter difference_selection needs to be TRUE or FALSE.")
  }
  if(difference_selection == TRUE) {
    inclusion_probability_difference = matrix(0,
      nrow = ncol(x),
      ncol = ncol(x)
    )

    difference_prior = match.arg(difference_prior)
    if(difference_prior == "Bernoulli") {
      if(length(difference_probability) == 1) {
        difference_inclusion_probability = difference_probability[1]
        if(is.na(difference_inclusion_probability) || is.null(difference_inclusion_probability)) {
          stop("There is no value specified for the inclusion probability for the differences.")
        }
        if(difference_inclusion_probability <= 0) {
          stop("The inclusion probability for differences needs to be positive.")
        }
        if(difference_inclusion_probability >= 1) {
          stop("The inclusion probability for differences cannot equal or exceed the value one.")
        }

        inclusion_probability_difference = matrix(difference_probability,
          nrow = ncol(x),
          ncol = ncol(x)
        )
      } else {
        if(!inherits(difference_probability, what = "matrix") &&
          !inherits(difference_probability, what = "data.frame")) {
          stop("The input for the inclusion probability argument for differences needs to be a single number, matrix, or dataframe.")
        }

        if(inherits(difference_probability, what = "data.frame")) {
          inclusion_probability_difference = data.matrix(difference_probability)
        } else {
          inclusion_probability_difference = difference_probability
        }

        if(!isSymmetric(inclusion_probability_difference)) {
          stop("The inclusion probability matrix needs to be symmetric.")
        }
        if(ncol(inclusion_probability_difference) != ncol(x)) {
          stop(paste0(
            "The inclusion probability matrix needs to have as many rows (columns) as there\n",
            " are variables in the data."
          ))
        }

        if(anyNA(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)]) ||
          any(is.null(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)]))) {
          stop("One or more inclusion probabilities for differences are not specified.")
        }
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)] <= 0)) {
          stop("One or more inclusion probabilities for differences are negative or zero.")
        }
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)] >= 1)) {
          stop("One or more inclusion probabilities for differences are one or larger.")
        }
      }
    } else {
      inclusion_probability_difference = matrix(0.5,
        nrow = ncol(x),
        ncol = ncol(x)
      )
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0) {
        stop("The scale parameters of the beta distribution for the differences need to be positive.")
      }
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta)) {
        stop("The scale parameters of the beta distribution for the differences need to be finite.")
      }
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
        is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta)) {
        stop("The scale parameters of the beta distribution for the differences need to be specified.")
      }
    }
  } else {
    difference_prior = "Not applicable"
    inclusion_probability_difference = matrix(0.5, 1, 1)
  }

  return(
    list(
      x = x,
      group_indicator = group,
      variable_bool = variable_bool,
      baseline_category = baseline_category,
      difference_prior = difference_prior,
      inclusion_probability_difference = inclusion_probability_difference
    )
  )
}

progress_type_from_display_progress <- function(display_progress = c("per-chain", "total", "none")) {
  if(is.logical(display_progress) && length(display_progress) == 1) {
    if(is.na(display_progress)) {
      stop("The display_progress argument must be a single logical value, but not NA.")
    }
    display_progress = if(display_progress) "per-chain" else "none"
  } else {
    display_progress = match.arg(display_progress)
  }
  return(if(display_progress == "per-chain") 2L else if(display_progress == "total") 1L else 0L)
}
