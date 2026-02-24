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
  dp <- validate_difference_prior(
    difference_selection  = difference_selection,
    difference_prior      = difference_prior,
    difference_probability = difference_probability,
    num_variables         = ncol(x),
    beta_bernoulli_alpha  = beta_bernoulli_alpha,
    beta_bernoulli_beta   = beta_bernoulli_beta
  )

  return(
    list(
      x = x,
      group_indicator = group,
      variable_bool = variable_bool,
      baseline_category = baseline_category,
      difference_prior = dp$difference_prior,
      inclusion_probability_difference = dp$inclusion_probability_difference
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
