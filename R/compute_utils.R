# ==============================================================================
# Compute utilities
# ==============================================================================
#
# Pure helper functions for computing derived quantities from data
# (scaling factors, category counts, sufficient statistics). Each
# function is pure: input -> output (no side-effects).
# ==============================================================================


# ------------------------------------------------------------------------------
# compute_scaling_factors
# ------------------------------------------------------------------------------
#
# Computes pairwise scaling factors for standardized Cauchy priors.
# When standardize = FALSE the result is a ones matrix (no scaling).
# When standardize = TRUE we compute scaling factors based on the
# maximum product of category ranges for each pair of variables.
#
# @param num_variables  Integer: number of variables.
# @param is_ordinal     Logical vector of length num_variables:
#   TRUE = regular ordinal (range 0..M), FALSE = Blume-Capel (range -b..M-b).
# @param num_categories Integer vector: number of response categories per
#   variable (i.e. max score M, where responses are 0..M).
# @param baseline_category Integer vector: baseline category for Blume-Capel
#   variables. Only used for entries where is_ordinal is FALSE.
# @param standardize    Logical scalar: whether to apply standardization.
# @param varnames       Character vector of variable names (used for row/col
#   names). NULL produces default names "Variable 1", "Variable 2", ...
#
# Returns:
#   A num_variables x num_variables matrix of scaling factors with named
#   rows and columns.
# ------------------------------------------------------------------------------
compute_scaling_factors = function(num_variables,
                                   is_ordinal,
                                   num_categories,
                                   baseline_category,
                                   standardize,
                                   varnames = NULL) {
  pairwise_scaling_factors = matrix(1,
    nrow = num_variables,
    ncol = num_variables
  )

  if(standardize) {
    for(v1 in seq_len(num_variables - 1)) {
      for(v2 in seq(v1 + 1, num_variables)) {
        if(is_ordinal[v1] && is_ordinal[v2]) {
          # Both ordinal: M_i * M_j (range 0..M)
          pairwise_scaling_factors[v1, v2] =
            num_categories[v1] * num_categories[v2]
        } else if(!is_ordinal[v1] && !is_ordinal[v2]) {
          # Both Blume-Capel: max of absolute endpoint products
          b1 = baseline_category[v1]
          b2 = baseline_category[v2]
          m1 = num_categories[v1]
          m2 = num_categories[v2]
          endpoints1 = c(-b1, m1 - b1)
          endpoints2 = c(-b2, m2 - b2)
          all_products = abs(outer(endpoints1, endpoints2))
          pairwise_scaling_factors[v1, v2] = max(all_products)
        } else {
          # Mixed: one ordinal, one Blume-Capel
          if(is_ordinal[v1]) {
            m1 = num_categories[v1]
            b2 = baseline_category[v2]
            m2 = num_categories[v2]
            endpoints1 = c(0, m1)
            endpoints2 = c(-b2, m2 - b2)
          } else {
            b1 = baseline_category[v1]
            m1 = num_categories[v1]
            m2 = num_categories[v2]
            endpoints1 = c(-b1, m1 - b1)
            endpoints2 = c(0, m2)
          }
          all_products = abs(outer(endpoints1, endpoints2))
          pairwise_scaling_factors[v1, v2] = max(all_products)
        }
        pairwise_scaling_factors[v2, v1] = pairwise_scaling_factors[v1, v2]
      }
    }
  }

  # Label rows and columns
  if(is.null(varnames)) {
    varnames = paste0("Variable ", seq_len(num_variables))
  }
  rownames(pairwise_scaling_factors) = varnames
  colnames(pairwise_scaling_factors) = varnames

  pairwise_scaling_factors
}


# ------------------------------------------------------------------------------
# compute_counts_per_category
# ------------------------------------------------------------------------------
#
# Compute per-group category counts for each variable. Used to build
# precomputed structures for the bgmCompare C++ backend.
#
# @param x  Numeric matrix: the recoded data.
# @param num_categories  Integer vector: max category per variable.
# @param group  Integer vector: group membership.
#
# Returns: list of matrices (one per group), each max_cat x num_variables.
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# compute_blume_capel_stats
# ------------------------------------------------------------------------------
#
# Compute sufficient statistics for Blume-Capel variables (linear and
# quadratic deviations from baseline). Used to build precomputed
# structures for the bgmCompare C++ backend.
#
# @param x  Numeric matrix: the recoded data.
# @param baseline_category  Integer vector: baseline categories.
# @param ordinal_variable  Logical vector: TRUE = ordinal, FALSE = BC.
# @param group  Integer vector or NULL: group membership.
#
# Returns: matrix (one-group) or list of matrices (multi-group),
#   each 2 x num_variables (row 1 = linear, row 2 = quadratic).
# ------------------------------------------------------------------------------
compute_blume_capel_stats = function(x, baseline_category, ordinal_variable, group = NULL) {
  if(is.null(group)) { # One-group design
    sufficient_stats = matrix(0, nrow = 2, ncol = ncol(x))
    bc_vars = which(!ordinal_variable)
    for(i in bc_vars) {
      sufficient_stats[1, i] = sum(x[, i] - baseline_category[i])
      sufficient_stats[2, i] = sum((x[, i] - baseline_category[i])^2)
    }
    return(sufficient_stats)
  } else { # Multi-group design
    sufficient_stats = list()
    for(g in unique(group)) {
      sufficient_stats_gr = matrix(0, nrow = 2, ncol = ncol(x))
      bc_vars = which(!ordinal_variable)
      for(i in bc_vars) {
        sufficient_stats_gr[1, i] = sum(x[group == g, i] - baseline_category[i])
        sufficient_stats_gr[2, i] = sum((x[group == g, i] - baseline_category[i])^2)
      }
      sufficient_stats[[length(sufficient_stats) + 1]] = sufficient_stats_gr
    }
    return(sufficient_stats)
  }
}


# ------------------------------------------------------------------------------
# compute_pairwise_stats
# ------------------------------------------------------------------------------
#
# Compute sufficient statistics for pairwise interactions (cross-product
# of observations per group). Used to build precomputed structures for
# the bgmCompare C++ backend.
#
# @param x  Numeric matrix: the centered data.
# @param group  Integer vector: group membership.
#
# Returns: list of p x p cross-product matrices (one per group).
# ------------------------------------------------------------------------------
compute_pairwise_stats = function(x, group) {
  result = list()

  for(g in unique(group)) {
    obs = x[group == g, , drop = FALSE]
    # cross-product: gives number of co-occurrences of categories
    result[[length(result) + 1]] = t(obs) %*% obs
  }

  result
}
