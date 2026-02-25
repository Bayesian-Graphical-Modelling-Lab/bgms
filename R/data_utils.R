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
