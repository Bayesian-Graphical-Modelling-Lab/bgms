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

  # Recode ordinal / Blume-Capel variables -------------------------------------
  ord <- reformat_ordinal_data(
    x                 = x,
    is_ordinal        = variable_bool,
    baseline_category = baseline_category
  )
  x                 <- ord$x
  num_categories    <- ord$num_categories
  baseline_category <- ord$baseline_category

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

  # Recode ordinal / Blume-Capel variables (same as single-group) --------------
  ord <- reformat_ordinal_data(
    x                 = x,
    is_ordinal        = variable_bool,
    baseline_category = baseline_category
  )
  x                 <- ord$x
  num_categories    <- ord$num_categories
  baseline_category <- ord$baseline_category

  # Collapse categories not observed in all groups (compare-specific) ----------
  col <- collapse_categories_across_groups(
    x                 = x,
    group             = group,
    is_ordinal        = variable_bool,
    num_categories    = num_categories,
    baseline_category = baseline_category
  )
  x                 <- col$x
  num_categories    <- col$num_categories
  baseline_category <- col$baseline_category

  return(list(
    x = x,
    group = group,
    num_categories = num_categories,
    baseline_category = baseline_category,
    missing_index = missing_index,
    na_impute = na_impute
  ))
}
