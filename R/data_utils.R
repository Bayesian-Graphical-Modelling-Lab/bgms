reformat_data = function(x,
                         na_action,
                         variable_bool,
                         baseline_category) {
  # Handle missing data --------------------------------------------------------
  md <- validate_missing_data(
    x             = x,
    na_action     = na_action,
    is_continuous = FALSE
  )
  x             <- md$x
  na_impute     <- md$na_impute
  missing_index <- md$missing_index

  check_fail_zero = FALSE
  num_variables = ncol(x)
  num_categories = vector(length = num_variables)
  for(node in 1:num_variables) {
    unq_vls = sort(unique(x[, node]))
    mx_vl = max(unq_vls)

    # Check if observed responses are not all unique ---------------------------
    if(mx_vl == nrow(x)) {
      stop(paste0(
        "Only unique responses observed for variable ",
        node,
        ". We expect >= 1 observations per category."
      ))
    }

    # Recode data --------------------------------------------------------------
    if(variable_bool[node]) { # Regular ordinal variable
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      if(length(unq_vls) != mx_vl + 1 || any(unq_vls != 0:mx_vl)) {
        y = x[, node]
        cntr = 0
        for(value in unq_vls) {
          x[y == value, node] = cntr
          cntr = cntr + 1
        }
      }
    } else { # Blume-Capel ordinal variable
      # Check if observations are integer or can be recoded --------------------
      if(any(abs(unq_vls - round(unq_vls)) > .Machine$double.eps)) {
        int_unq_vls = unique(as.integer(unq_vls))
        if(anyNA(int_unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers, but \n",
            "the category scores for node ", node, " were not integer. An attempt to recode \n",
            "them to integer failed. Please inspect the documentation for the base R \n",
            "function as.integer(), which bgm uses for recoding category scores."
          ))
        }

        if(length(int_unq_vls) != length(unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers. The \n",
            "category scores of the observations for node ", node, " were not integers. An \n",
            "attempt to recode these observations as integers failed because, after rounding, \n",
            "a single integer value was used for several observed score categories."
          ))
        }
        x[, node] = as.integer(x[, node])

        if(baseline_category[node] < 0 | baseline_category[node] > max(x[, node])) {
          stop(paste0(
            "The reference category for the Blume-Capel variable ", node, "is outside its \n",
            "range of observations."
          ))
        }
      }

      # Check if observations start at zero and recode otherwise ---------------
      if(min(x[, node]) != 0) {
        baseline_category[node] = baseline_category[node] - min(x[, node])
        x[, node] = x[, node] - min(x[, node])

        if(check_fail_zero == FALSE) {
          check_fail_zero = TRUE
          failed_zeroes = c(node)
        } else {
          failed_zeroes = c(failed_zeroes, node)
        }
      }

      check_range = length(unique(x[, node]))
      if(check_range < 3) {
        stop(paste0(
          "The Blume-Capel is only available for variables with more than one category \n",
          "observed. There two or less categories observed for variable ",
          node,
          "."
        ))
      }
    }

    # Warn that maximum category value is large --------------------------------
    num_categories[node] = max(x[, node])
    if(!variable_bool[node] & num_categories[node] > 10) {
      warning(
        "Blume-Capel variable ", node, " has ", num_categories[node], " categories. ",
        "This may slow computation. Empty categories are not collapsed.",
        call. = FALSE
      )
    }

    # Check to see if not all responses are in one category --------------------
    if(num_categories[node] == 0) {
      stop(paste0(
        "Only one value [",
        unq_vls,
        "] was observed for variable ",
        node,
        "."
      ))
    }
  }

  if(check_fail_zero == TRUE && isTRUE(getOption("bgms.verbose", TRUE))) {
    nodes_str <- paste(failed_zeroes, collapse = ", ")
    message(
      "Variable", if(length(failed_zeroes) > 1) "s" else "", " ", nodes_str,
      " recoded to start at 0 (baseline categor",
      if(length(failed_zeroes) > 1) "ies" else "y", " adjusted)."
    )
  }

  return(list(
    x = x,
    num_categories = num_categories,
    baseline_category = baseline_category,
    missing_index = missing_index,
    na_impute = na_impute
  ))
}

# Helper function for data checks
data_check = function(data, name) {
  if(!inherits(data, c("matrix", "data.frame"))) {
    stop(paste(name, "must be a matrix or data.frame."))
  }
  if(inherits(data, "data.frame")) {
    data = data.matrix(data)
  }
  if(nrow(data) < 2 || ncol(data) < 2) {
    stop(paste(name, "must have at least 2 rows and 2 columns."))
  }
  return(data)
}

# Helper function for computing `counts_per_category`
compute_counts_per_category = function(x, num_categories, group = NULL) {
  counts_per_category = list()
  for(g in unique(group)) {
    counts_per_category_gr = matrix(0, nrow = max(num_categories), ncol = ncol(x))
    for(variable in seq_len(ncol(x))) {
      for(category in seq_len(num_categories[variable])) {
        counts_per_category_gr[category, variable] = sum(x[group == g, variable] == category)
      }
    }
    counts_per_category[[length(counts_per_category) + 1]] = counts_per_category_gr
  }
  return(counts_per_category)
}

# Helper function for computing sufficient statistics for Blume-Capel variables
compute_blume_capel_stats = function(x, baseline_category, ordinal_variable, group = NULL) {
  if(is.null(group)) { # One-group design
    sufficient_stats = matrix(0, nrow = 2, ncol = ncol(x))
    bc_vars = which(!ordinal_variable)
    for (i in bc_vars) {
      sufficient_stats[1, i] = sum(x[, i] - baseline_category[i])
      sufficient_stats[2, i] = sum((x[, i] - baseline_category[i]) ^ 2)
    }
    return(sufficient_stats)
  } else { # Multi-group design
    sufficient_stats = list()
    for(g in unique(group)) {
      sufficient_stats_gr = matrix(0, nrow = 2, ncol = ncol(x))
      bc_vars = which(!ordinal_variable)
      for (i in bc_vars) {
        sufficient_stats_gr[1, i] = sum(x[group == g, i] - baseline_category[i])
        sufficient_stats_gr[2, i] = sum((x[group == g, i] - baseline_category[i]) ^ 2)
      }
      sufficient_stats[[length(sufficient_stats) + 1]] = sufficient_stats_gr
    }
    return(sufficient_stats)
  }
}

# Helper function for computing sufficient statistics for pairwise interactions
compute_pairwise_stats <- function(x, group) {
  result <- list()

  for(g in unique(group)) {
    obs <- x[group == g, , drop = FALSE]
    # cross-product: gives number of co-occurrences of categories
    result[[length(result) + 1]] <- t(obs) %*% obs
  }

  result
}


compare_reformat_data = function(
  x,
  group,
  na_action,
  variable_bool,
  baseline_category
) {
  # Handle missing data --------------------------------------------------------
  md <- validate_missing_data(
    x             = x,
    na_action     = na_action,
    is_continuous = FALSE,
    group         = group
  )
  x             <- md$x
  na_impute     <- md$na_impute
  missing_index <- md$missing_index
  group         <- md$group

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

  check_fail_zero = FALSE
  num_variables = ncol(x)
  num_categories = vector(length = num_variables)

  for(node in 1:num_variables) {
    unq_vls = sort(unique(x[, node]))
    mx_vls = length(unq_vls)

    # Check if observed responses are not all unique ---------------------------
    if(mx_vls == nrow(x)) {
      stop(paste0(
        "Only unique responses observed for variable ",
        node,
        " in the matrix x (group 1). We expect >= 1 observations per category."
      ))
    }

    # Recode data --------------------------------------------------------------
    if(variable_bool[node]) { # Regular ordinal variable
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      observed_scores = matrix(NA,
        nrow = mx_vls,
        ncol = max(group)
      )

      for(value in unq_vls) {
        unique_g = unique(group)
        for(g in unique_g) {
          observed_scores[which(unq_vls == value), g] =
            any(x[group == g, node] == value) * 1
        }
      }

      xx = x[, node]
      cntr = -1
      for(value in unq_vls) {
        # Collapse categories when not observed in one or more groups.
        if(sum(observed_scores[which(unq_vls == value), ]) == max(group)) {
          cntr = cntr + 1 # increment score if category observed in all groups
        }
        x[xx == value, node] = max(0, cntr)
      }
    } else { # Blume-Capel ordinal variable
      # Check if observations are integer or can be recoded --------------------
      if(any(abs(unq_vls - round(unq_vls)) > .Machine$double.eps)) {
        int_unq_vls = unique(as.integer(unq_vls))
        if(anyNA(int_unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers, but \n",
            "the category scores for node ", node, " were not integer. An attempt to recode \n",
            "them to integer failed. Please inspect the documentation for the base R function \n",
            "as.integer(), which bgmCompare uses for recoding category scores."
          ))
        }

        if(length(int_unq_vls) != length(unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers. The \n",
            "category scores of the observations for node ", node, " were not integers. An \n",
            "attempt to recode these observations as integers failed because, after rounding,\n",
            "a single integer value was used for several observed score categories."
          ))
        }
        x[, node] = as.integer(x[, node])
      }

      mi = min(x[, node])

      ma = max(x[, node])

      if(baseline_category[node] < mi | baseline_category[node] > ma) {
        stop(paste0(
          "The reference category for the Blume-Capel variable ", node, "is outside its \n",
          "range of observations in the matrices x (and y)."
        ))
      }

      # Check if observations start at zero and recode otherwise ---------------
      if(mi != 0) {
        baseline_category[node] = baseline_category[node] - mi
        x[, node] = x[, node] - mi

        if(check_fail_zero == FALSE) {
          check_fail_zero = TRUE
          failed_zeroes = c(node)
        } else {
          failed_zeroes = c(failed_zeroes, node)
        }
      }

      check_range = length(unique(x[, node]))

      if(check_range < 3) {
        stop(paste0(
          "The Blume-Capel is only available for variables with more than two categories \n",
          "observed. There are two or less categories observed for variable ",
          node,
          "."
        ))
      }
    }


    # Warn that maximum category value is large --------------------------------

    num_categories[node] = max(x[, node])

    if(!variable_bool[node] & max(num_categories[node]) > 10) {
      warning(
        "Blume-Capel variable ", node, " has ", max(num_categories[node]), " categories. ",
        "This may slow computation. Empty categories are not collapsed.",
        call. = FALSE
      )
    }


    # Check to see if not all responses are in one category --------------------
    if(any(num_categories[node] == 0)) {
      stop(paste0("Only one value was observed for variable ", node, "."))
    }
  }


  if(check_fail_zero == TRUE && isTRUE(getOption("bgms.verbose", TRUE))) {
    nodes_str <- paste(failed_zeroes, collapse = ", ")
    message(
      "Variable", if(length(failed_zeroes) > 1) "s" else "", " ", nodes_str,
      " recoded to start at 0 (baseline categor",
      if(length(failed_zeroes) > 1) "ies" else "y", " adjusted)."
    )
  }

  return(list(
    x = x,
    group = group,
    num_categories = num_categories,
    baseline_category = baseline_category,
    missing_index = missing_index,
    na_impute = na_impute
  ))
}
