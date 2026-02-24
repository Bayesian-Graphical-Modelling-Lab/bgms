# ==============================================================================
# Compute utilities
# ==============================================================================
#
# Extracted from inline code in bgm.R and bgmCompare.R as part of
# Phase A.8 of the R scaffolding refactor.
#
# Each function is pure: input -> output (no side-effects).
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
compute_scaling_factors <- function(num_variables,
                                    is_ordinal,
                                    num_categories,
                                    baseline_category,
                                    standardize,
                                    varnames = NULL) {

  pairwise_scaling_factors <- matrix(1,
                                     nrow = num_variables,
                                     ncol = num_variables)

  if (standardize) {
    for (v1 in seq_len(num_variables - 1)) {
      for (v2 in seq(v1 + 1, num_variables)) {
        if (is_ordinal[v1] && is_ordinal[v2]) {
          # Both ordinal: M_i * M_j (range 0..M)
          pairwise_scaling_factors[v1, v2] <-
            num_categories[v1] * num_categories[v2]

        } else if (!is_ordinal[v1] && !is_ordinal[v2]) {
          # Both Blume-Capel: max of absolute endpoint products
          b1 <- baseline_category[v1]
          b2 <- baseline_category[v2]
          m1 <- num_categories[v1]
          m2 <- num_categories[v2]
          endpoints1 <- c(-b1, m1 - b1)
          endpoints2 <- c(-b2, m2 - b2)
          all_products <- abs(outer(endpoints1, endpoints2))
          pairwise_scaling_factors[v1, v2] <- max(all_products)

        } else {
          # Mixed: one ordinal, one Blume-Capel
          if (is_ordinal[v1]) {
            m1 <- num_categories[v1]
            b2 <- baseline_category[v2]
            m2 <- num_categories[v2]
            endpoints1 <- c(0, m1)
            endpoints2 <- c(-b2, m2 - b2)
          } else {
            b1 <- baseline_category[v1]
            m1 <- num_categories[v1]
            m2 <- num_categories[v2]
            endpoints1 <- c(-b1, m1 - b1)
            endpoints2 <- c(0, m2)
          }
          all_products <- abs(outer(endpoints1, endpoints2))
          pairwise_scaling_factors[v1, v2] <- max(all_products)
        }
        pairwise_scaling_factors[v2, v1] <- pairwise_scaling_factors[v1, v2]
      }
    }
  }

  # Label rows and columns
  if (is.null(varnames)) {
    varnames <- paste0("Variable ", seq_len(num_variables))
  }
  rownames(pairwise_scaling_factors) <- varnames
  colnames(pairwise_scaling_factors) <- varnames

  pairwise_scaling_factors
}
