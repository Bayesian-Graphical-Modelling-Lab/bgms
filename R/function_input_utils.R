check_positive_integer = function(value, name) {
  if (!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value <= 0) {
    stop(sprintf("Parameter `%s` must be a positive integer. Got: %s", name, value))
  }
}

# Helper function for validating non-negative integers
check_non_negative_integer = function(value, name) {
  if (!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value < 0) {
    stop(sprintf("Parameter `%s` must be a non-negative integer. Got: %s", name, value))
  }
}

# Helper function for validating logical inputs
check_logical = function(value, name) {
  value = as.logical(value)
  if (is.na(value)) {
    stop(sprintf("Parameter `%s` must be TRUE or FALSE. Got: %s", name, value))
  }
  return(value)
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

  #Check variable type input ---------------------------------------------------
  if(length(variable_type) == 1) {
    variable_input = variable_type
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel")),
                        silent = TRUE)
    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  variable_input, "."))
    variable_bool = (variable_type == "ordinal")
    variable_bool = rep(variable_bool, ncol(x))
  } else {
    if(length(variable_type) != ncol(x))
      stop(paste0("The variable type vector variable_type should be either a single character\n",
                  "string or a vector of character strings of length p."))

    variable_input = unique(variable_type)
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel"),
                                  several.ok = TRUE), silent = TRUE)

    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input, collapse = ", "), "."))

    num_types = sapply(variable_input, function(type) {
      tmp = try(match.arg(arg = type,
                          choices = c("ordinal", "blume-capel")),
                silent = TRUE)
      inherits(tmp, what = "try-error")
    })

    if(length(variable_type) != ncol(x))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input[num_types], collapse = ", "), "."))

    variable_bool = (variable_type == "ordinal")
  }

  #Check Blume-Capel variable input --------------------------------------------
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    if(!hasArg("baseline_category"))
      stop("The argument baseline_category is required for Blume-Capel variables.")

    if(length(baseline_category) != ncol(x) && length(baseline_category) != 1)
      stop(paste0("The argument baseline_category for the Blume-Capel model needs to be a \n",
                  "single integer or a vector of integers of length p."))

    if(length(baseline_category) == 1) {
      #Check if the input is integer -------------------------------------------
      integer_check = try(as.integer(baseline_category), silent = TRUE)
      if(is.na(integer_check))
        stop(paste0("The baseline_category argument for the Blume-Capel model contains either \n",
                    "a missing value or a value that could not be forced into an integer value."))
      integer_check = baseline_category - round(baseline_category)
      if(integer_check > .Machine$double.eps)
        stop("Reference category needs to an integer value or a vector of integers of length p.")
      baseline_category = rep.int(baseline_category, times = ncol(x))
    }

    #Check if the input is integer -------------------------------------------
    blume_capel_variables = which(!variable_bool)
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    integer_check = try(as.integer(baseline_category[blume_capel_variables]),
                        silent = TRUE)
    if(anyNA(integer_check))
      stop(paste0("The baseline_category argument for the Blume-Capel model contains either \n",
                  "missing values or values that could not be forced into an integer value."))

    integer_check = baseline_category[blume_capel_variables] -
      round(baseline_category[blume_capel_variables])

    if(any(integer_check > .Machine$double.eps)) {
      non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
      if(length(non_integers) > 1) {
        stop(paste0("The entries in baseline_category for variables ",
                    paste0(non_integers, collapse = ", "), " need to be integer."))
      } else {
        stop(paste0("The entry in baseline_category for variable ",
                    non_integers, " needs to be an integer."))
      }
    }

    variable_lower = apply(x, 2, min, na.rm = TRUE)
    variable_upper = apply(x, 2, max, na.rm = TRUE)

    if(any(baseline_category < variable_lower) | any(baseline_category > variable_upper)) {
      out_of_range = which(baseline_category < variable_lower | baseline_category > variable_upper)
      stop(paste0("The Blume-Capel model assumes that the reference category is within the range \n",
                  "of the observed category scores. This was not the case for variable(s) \n",
                  paste0(out_of_range, collapse =", "),
                  "."))
    }

  } else {
    baseline_category = rep.int(0, times = ncol(x))
  }

  #Check prior set-up for the interaction parameters ---------------------------
  if(pairwise_scale <= 0 || is.na(pairwise_scale) || is.infinite(pairwise_scale))
    stop("The scale of the Cauchy prior needs to be positive.")

  #Check prior set-up for the threshold parameters -----------------------------
  if(main_alpha <= 0 | !is.finite(main_alpha))
    stop("Parameter main_alpha needs to be positive.")
  if(main_beta <= 0 | !is.finite(main_beta))
    stop("Parameter main_beta needs to be positive.")

  #Check set-up for the Bayesian edge selection model --------------------------
  edge_selection = as.logical(edge_selection)
  if(is.na(edge_selection))
    stop("The parameter edge_selection needs to be TRUE or FALSE.")
  if(edge_selection == TRUE) {
    #Check prior set-up for the edge indicators --------------------------------
    edge_prior = match.arg(edge_prior)
    if(edge_prior == "Bernoulli") {
      if(length(inclusion_probability) == 1) {
        theta = inclusion_probability[1]
        if(is.na(theta) || is.null(theta))
          stop("There is no value specified for the inclusion probability.")
        if(theta <= 0)
          stop("The inclusion probability needs to be positive.")
        if(theta > 1)
          stop("The inclusion probability cannot exceed the value one.")
        if(theta == 1)
          stop("The inclusion probability cannot equal one.")

        theta = matrix(theta, nrow = ncol(x), ncol = ncol(x))
      } else {
        if(!inherits(inclusion_probability, what = "matrix") &&
           !inherits(inclusion_probability, what = "data.frame"))
          stop("The input for the inclusion probability argument needs to be a single number, matrix, or dataframe.")

        if(inherits(inclusion_probability, what = "data.frame")) {
          theta = data.matrix(inclusion_probability)
        } else {
          theta = inclusion_probability
        }
        if(!isSymmetric(theta))
          stop("The inclusion probability matrix needs to be symmetric.")
        if(ncol(theta) != ncol(x))
          stop("The inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

        if(anyNA(theta[lower.tri(theta)]) ||
           any(is.null(theta[lower.tri(theta)])))
          stop("One or more elements of the elements in inclusion probability matrix are not specified.")
        if(any(theta[lower.tri(theta)] <= 0))
          stop(paste0("The inclusion probability matrix contains negative or zero values;\n",
                      "inclusion probabilities need to be positive."))
        if(any(theta[lower.tri(theta)] >= 1))
          stop(paste0("The inclusion probability matrix contains values greater than or equal to one;\n",
                      "inclusion probabilities cannot exceed or equal the value one."))
      }
    }
    if(edge_prior == "Beta-Bernoulli") {
      theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0)
        stop("The scale parameters of the beta distribution need to be positive.")
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta))
        stop("The scale parameters of the beta distribution need to be finite.")
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
         is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta))
        stop("Values for both scale parameters of the beta distribution need to be specified.")
    }

    if(edge_prior == "Stochastic-Block") {
      theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))

      # Check that all beta parameters are provided
      if (is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
          is.null(beta_bernoulli_alpha_between) || is.null(beta_bernoulli_beta_between)) {
        stop("The Stochastic-Block prior requires all four beta parameters: ",
             "beta_bernoulli_alpha, beta_bernoulli_beta, ",
             "beta_bernoulli_alpha_between, and beta_bernoulli_beta_between.")
      }

      # Check that all beta parameters are positive
      if (beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0 ||
          beta_bernoulli_alpha_between <= 0 || beta_bernoulli_beta_between <= 0 ||
          dirichlet_alpha <= 0 || lambda <= 0) {
        stop("The parameters of the beta and Dirichlet distributions need to be positive.")
      }

      # Check that all beta parameters are finite
      if (!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta) ||
          !is.finite(beta_bernoulli_alpha_between) || !is.finite(beta_bernoulli_beta_between) ||
          !is.finite(dirichlet_alpha) || !is.finite(lambda)) {
        stop("The shape parameters of the beta distribution, the concentration parameter of the Dirichlet distribution, ",
             "and the rate parameter of the Poisson distribution need to be finite.")
      }

      # Check for NAs
      if (is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
          is.na(beta_bernoulli_alpha_between) || is.na(beta_bernoulli_beta_between) ||
          is.na(dirichlet_alpha) || is.na(lambda)) {
        stop("Values for all shape parameters of the beta distribution, the concentration parameter of the Dirichlet distribution, ",
             "and the rate parameter of the Poisson distribution cannot be NA.")
      }
    }
 }else {
    theta = matrix(0.5, nrow = 1, ncol = 1)
    edge_prior = "Not Applicable"
  }

  return(list(variable_bool = variable_bool,
              baseline_category = baseline_category,
              edge_selection = edge_selection,
              edge_prior = edge_prior,
              inclusion_probability = theta))
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
    if(length(unique_g) == 0)
      stop(paste0("The bgmCompare function expects at least two groups, but the input group_indicator contains\n",
                  "no group value."))
    if(length(unique_g) == 1)
      stop(paste0("The bgmCompare function expects at least two groups, but the input group_indicator contains\n",
                  "only one group value."))
    if(length(unique_g) == length(group_indicator))
      stop("The input group_indicator contains only unique group values.")

    group = group_indicator
    for(u in unique_g) {
      group[group_indicator == u] = which(unique_g == u)
    }
    tab = tabulate(group)

    if(any(tab < 2))
      stop("One or more groups only had one member in the input group_indicator.")
  } else {
    group = c(rep.int(1, times = nrow(x)), rep.int(2, times = nrow(y)))
    x = rbind(x, y)
  }

  #Check variable type input ---------------------------------------------------
  if(length(variable_type) == 1) {
    variable_input = variable_type
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel")),
                        silent = TRUE)
    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgmCompare function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  variable_input, "."))
    variable_bool = (variable_type == "ordinal")
    variable_bool = rep(variable_bool, ncol(x))
  } else {
    if(length(variable_type) != ncol(x))
      stop(paste0("The variable type vector variable_type should be either a single character\n",
                  "string or a vector of character strings of length p."))

    variable_input = unique(variable_type)
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel"),
                                  several.ok = TRUE), silent = TRUE)

    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgmCompare function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input, collapse = ", "), "."))

    num_types = sapply(variable_input, function(type) {
      tmp = try(match.arg(arg = type,
                          choices = c("ordinal", "blume-capel")),
                silent = TRUE)
      inherits(tmp, what = "try-error")
    })

    if(length(variable_type) != ncol(x))
      stop(paste0("The bgmCompare function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input[num_types], collapse = ", "), "."))

    variable_bool = (variable_type == "ordinal")
  }

  #Check Blume-Capel variable input --------------------------------------------
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    if(!hasArg("baseline_category"))
      stop("The argument baseline_category is required for Blume-Capel variables.")

    if(length(baseline_category) != ncol(x) && length(baseline_category) != 1)
      stop(paste0("The argument baseline_category for the Blume-Capel model needs to be a \n",
                  "single integer or a vector of integers of length p."))

    if(length(baseline_category) == 1) {
      #Check if the input is integer -------------------------------------------
      integer_check = try(as.integer(baseline_category), silent = TRUE)
      if(is.na(integer_check))
        stop(paste0("The baseline_category argument for the Blume-Capel model contains either \n",
                    "a missing value or a value that could not be forced into an integer value."))
      integer_check = baseline_category - round(baseline_category)
      if(integer_check > .Machine$double.eps)
        stop("Reference category needs to an integer value or a vector of integers of length p.")
      baseline_category = rep.int(baseline_category, times = ncol(x))
    }

    #Check if the input is integer -------------------------------------------
    blume_capel_variables = which(!variable_bool)
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    integer_check = try(as.integer(baseline_category[blume_capel_variables]),
                        silent = TRUE)
    if(anyNA(integer_check))
      stop(paste0("The baseline_category argument for the Blume-Capel model contains either \n",
                  "missing values or values that could not be forced into an integer value."))

    integer_check = baseline_category[blume_capel_variables] -
      round(baseline_category[blume_capel_variables])

    if(any(integer_check > .Machine$double.eps)) {
      non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
      if(length(non_integers) > 1) {
        stop(paste0("The entries in baseline_category for variables ",
                    paste0(non_integers, collapse = ", "), " need to be integer."))
      } else {
        stop(paste0("The entry in baseline_category for variable ",
                    non_integers, " needs to be an integer."))
      }
    }

    variable_lower = apply(x, 2, min, na.rm = TRUE)
    variable_upper = apply(x, 2, max, na.rm = TRUE)

    if(any(baseline_category < variable_lower) | any(baseline_category > variable_upper)) {
      out_of_range = which(baseline_category < variable_lower | baseline_category > variable_upper)
      stop(paste0("The Blume-Capel model assumes that the reference category is within the range \n",
                  "of the observed category scores. This was not the case for variable(s) \n",
                  paste0(out_of_range, collapse =", "),
                  "."))
    }
  } else {
    baseline_category = rep.int(0, times = ncol(x))
  }

  #Check prior set-up for the interaction parameters ---------------------------
  if(pairwise_scale <= 0 || is.na(pairwise_scale) || is.infinite(pairwise_scale))
    stop("The scale of the Cauchy prior for the interactions needs to be positive.")

  #Check prior set-up for the interaction differences --------------------------
  if(difference_scale <= 0 || is.na(difference_scale) || is.infinite(difference_scale))
    stop("The scale of the Cauchy prior for the differences needs to be positive.")

  #Check prior set-up for the threshold parameters -----------------------------
  if(main_alpha <= 0 | !is.finite(main_alpha))
    stop("Parameter main_alpha needs to be positive.")
  if(main_beta <= 0 | !is.finite(main_beta))
    stop("Parameter main_beta needs to be positive.")

  #Check set-up for the Bayesian difference selection model --------------------
  difference_selection = as.logical(difference_selection)
  if(is.na(difference_selection))
    stop("The parameter difference_selection needs to be TRUE or FALSE.")
  if(difference_selection == TRUE) {
    inclusion_probability_difference = matrix(0,
                                              nrow = ncol(x),
                                              ncol = ncol(x))

    difference_prior = match.arg(difference_prior)
    if(difference_prior == "Bernoulli") {
      if(length(difference_probability) == 1) {
        difference_inclusion_probability = difference_probability[1]
        if(is.na(difference_inclusion_probability) || is.null(difference_inclusion_probability))
          stop("There is no value specified for the inclusion probability for the differences.")
        if(difference_inclusion_probability <= 0)
          stop("The inclusion probability for differences needs to be positive.")
        if(difference_inclusion_probability >= 1)
          stop("The inclusion probability for differences cannot equal or exceed the value one.")

       inclusion_probability_difference = matrix(difference_probability,
                                                 nrow = ncol(x),
                                                 ncol = ncol(x))

      } else {
        if(!inherits(difference_probability, what = "matrix") &&
           !inherits(difference_probability, what = "data.frame"))
          stop("The input for the inclusion probability argument for differences needs to be a single number, matrix, or dataframe.")

        if(inherits(difference_probability, what = "data.frame")) {
          inclusion_probability_difference = data.matrix(difference_probability)
        } else {
          inclusion_probability_difference = difference_probability
        }

        if(!isSymmetric(inclusion_probability_difference))
          stop("The inclusion probability matrix needs to be symmetric.")
        if(ncol(inclusion_probability_difference) != ncol(x))
          stop(paste0("The inclusion probability matrix needs to have as many rows (columns) as there\n",
                      " are variables in the data."))

        if(anyNA(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)]) ||
           any(is.null(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)])))
          stop("One or more inclusion probabilities for differences are not specified.")
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)] <= 0))
          stop("One or more inclusion probabilities for differences are negative or zero.")
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference, diag = TRUE)] >= 1))
          stop("One or more inclusion probabilities for differences are one or larger.")
      }
    } else {
      inclusion_probability_difference = matrix(0.5,
                                                nrow = ncol(x),
                                                ncol = ncol(x))
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0)
        stop("The scale parameters of the beta distribution for the differences need to be positive.")
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta))
        stop("The scale parameters of the beta distribution for the differences need to be finite.")
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
         is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta))
        stop("The scale parameters of the beta distribution for the differences need to be specified.")
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
  if (is.logical(display_progress) && length(display_progress) == 1) {
    if (is.na(display_progress))
      stop("The display_progress argument must be a single logical value, but not NA.")
    display_progress = if (display_progress) "per-chain" else "none"
  } else {
    display_progress = match.arg(display_progress)
  }
  return(if (display_progress == "per-chain") 2L else if (display_progress == "total") 1L else 0L)
}
