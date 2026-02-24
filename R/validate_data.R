# ==============================================================================
# Data validation functions
# ==============================================================================
#
# Extracted from reformat_data(), compare_reformat_data(), and inline GGM
# handling in bgm.R as part of the R scaffolding refactor
# (dev/scaffolding/plan.md, Phase A).
#
# Each validator is a pure function: input -> validated output (or error).
# ==============================================================================


# ------------------------------------------------------------------------------
# validate_missing_data
# ------------------------------------------------------------------------------
#
# Handles missing data for all model types (OMRF, GGM, bgmCompare).
# Shared by reformat_data(), compare_reformat_data(), and the GGM path
# in bgm.R.
#
# @param x  Numeric matrix: the data.
# @param na_action  Character: "listwise" or "impute".
# @param is_continuous  Logical: TRUE for GGM (continuous) models.
#   When TRUE and na_action == "impute", an error is raised
#   (imputation not yet supported for GGM).
# @param group  Optional integer vector: group indicators for bgmCompare.
#   If provided, listwise deletion also filters the group vector.
#   NULL for bgm() calls.
#
# Returns:
#   list(x, na_impute, missing_index)
#   - x: data matrix (rows removed for listwise, NAs imputed for impute)
#   - na_impute: logical
#   - missing_index: matrix(NA, 1, 1) if no imputation,
#       otherwise Nx2 matrix of 0-based (row, col) indices
#   - group: filtered group vector (only present when group was non-NULL)
#   - n_removed: integer count of rows removed (listwise only, 0 otherwise)
#
# Replaces:
#   - reformat_data()         lines 5-68  (data_utils.R)
#   - compare_reformat_data() lines 268-350 (data_utils.R)
#   - bgm.R                   lines 606-622 (inline GGM path)
# ------------------------------------------------------------------------------
validate_missing_data <- function(x,
                                  na_action,
                                  is_continuous = FALSE,
                                  group = NULL) {
  # --- GGM + impute guard ---
  if (is_continuous && na_action == "impute") {
    stop(
      "Imputation is not yet supported for the Gaussian model. ",
      "Use na_action = 'listwise'."
    )
  }

  if (na_action == "listwise") {
    return(handle_listwise(x, group))
  }

  # --- impute path ---
  handle_impute(x, group)
}


# ------------------------------------------------------------------------------
# handle_listwise (internal helper)
# ------------------------------------------------------------------------------
handle_listwise <- function(x, group = NULL) {
  missing_rows <- apply(x, 1, anyNA)

  if (all(missing_rows)) {
    stop(paste0(
      "All rows in x contain at least one missing response.\n",
      "You could try option na_action = impute."
    ))
  }

  n_removed <- sum(missing_rows)
  if (n_removed > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
    n_remaining <- nrow(x) - n_removed
    message(
      n_removed, " row", if (n_removed > 1) "s" else "",
      " with missing values excluded (n = ", n_remaining, " remaining).\n",
      "To impute missing values instead, use na_action = \"impute\"."
    )
  }

  x <- x[!missing_rows, , drop = FALSE]

  if (is.null(ncol(x)) || ncol(x) < 2) {
    stop(paste0(
      "After removing missing observations from the input matrix x,\n",
      "there were less than two columns left in x."
    ))
  }
  if (is.null(nrow(x)) || nrow(x) < 2) {
    stop(paste0(
      "After removing missing observations from the input matrix x,\n",
      "there were less than two rows left in x."
    ))
  }

  result <- list(
    x             = x,
    na_impute     = FALSE,
    missing_index = matrix(NA, nrow = 1, ncol = 1),
    n_removed     = n_removed
  )

  if (!is.null(group)) {
    result$group <- group[!missing_rows]
  }

  result
}


# ------------------------------------------------------------------------------
# handle_impute (internal helper)
# ------------------------------------------------------------------------------
handle_impute <- function(x, group = NULL) {
  num_missings <- sum(is.na(x))

  if (num_missings == 0) {
    result <- list(
      x             = x,
      na_impute     = FALSE,
      missing_index = matrix(NA, nrow = 1, ncol = 1),
      n_removed     = 0L
    )
    if (!is.null(group)) result$group <- group
    return(result)
  }

  num_variables <- ncol(x)
  missing_index <- matrix(0, nrow = num_missings, ncol = 2)
  cntr <- 0
  for (node in seq_len(num_variables)) {
    mis <- which(is.na(x[, node]))
    if (length(mis) > 0) {
      for (i in seq_along(mis)) {
        cntr <- cntr + 1
        missing_index[cntr, 1] <- mis[i] - 1  # C++ 0-based index
        missing_index[cntr, 2] <- node - 1     # C++ 0-based index
        x[mis[i], node] <- sample(x[-mis, node], size = 1)
      }
    }
  }

  result <- list(
    x             = x,
    na_impute     = TRUE,
    missing_index = missing_index,
    n_removed     = 0L
  )
  if (!is.null(group)) result$group <- group

  result
}
