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
# data_check
# ------------------------------------------------------------------------------
#
# Coerce user data to a numeric matrix and perform basic dimension checks.
# Originally in data_utils.R; moved here in Phase D.2.
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# reformat_ordinal_data
# ------------------------------------------------------------------------------
#
# Per-variable recoding of ordinal and Blume-Capel data to 0-based
# contiguous categories. Single-group only (no group-conditional
# collapsing — see collapse_categories_across_groups() for that).
#
# Extracted from the per-variable for-loop in reformat_data()
# (data_utils.R) as part of Phase A.6.
#
# @param x  Numeric matrix: the data (after missing-data handling).
# @param is_ordinal  Logical vector of length ncol(x): TRUE = regular
#   ordinal variable, FALSE = Blume-Capel variable.
# @param baseline_category  Integer vector of length ncol(x): baseline
#   (reference) categories for Blume-Capel variables.
#
# Returns:
#   list(x, num_categories, baseline_category)
#   - x: matrix with recoded values
#   - num_categories: integer vector (max observed value per variable)
#   - baseline_category: possibly adjusted baseline categories
# ------------------------------------------------------------------------------
reformat_ordinal_data <- function(x, is_ordinal, baseline_category) {
  check_fail_zero <- FALSE
  num_variables <- ncol(x)
  num_categories <- vector(length = num_variables)

  for (node in 1:num_variables) {
    unq_vls <- sort(unique(x[, node]))
    mx_vl <- max(unq_vls)

    # Check if observed responses are not all unique ---------------------------
    if (mx_vl == nrow(x)) {
      stop(paste0(
        "Only unique responses observed for variable ",
        node,
        ". We expect >= 1 observations per category."
      ))
    }

    # Recode data --------------------------------------------------------------
    if (is_ordinal[node]) { # Regular ordinal variable
      if (length(unq_vls) != mx_vl + 1 || any(unq_vls != 0:mx_vl)) {
        y <- x[, node]
        cntr <- 0
        for (value in unq_vls) {
          x[y == value, node] <- cntr
          cntr <- cntr + 1
        }
      }
    } else { # Blume-Capel ordinal variable
      # Check if observations are integer or can be recoded --------------------
      if (any(abs(unq_vls - round(unq_vls)) > .Machine$double.eps)) {
        int_unq_vls <- unique(as.integer(unq_vls))
        if (anyNA(int_unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers, but \n",
            "the category scores for node ", node, " were not integer. An attempt to recode \n",
            "them to integer failed. Please inspect the documentation for the base R \n",
            "function as.integer(), which bgm uses for recoding category scores."
          ))
        }

        if (length(int_unq_vls) != length(unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers. The \n",
            "category scores of the observations for node ", node, " were not integers. An \n",
            "attempt to recode these observations as integers failed because, after rounding, \n",
            "a single integer value was used for several observed score categories."
          ))
        }
        x[, node] <- as.integer(x[, node])

        if (baseline_category[node] < 0 | baseline_category[node] > max(x[, node])) {
          stop(paste0(
            "The reference category for the Blume-Capel variable ", node, "is outside its \n",
            "range of observations."
          ))
        }
      }

      # Check if observations start at zero and recode otherwise ---------------
      if (min(x[, node]) != 0) {
        baseline_category[node] <- baseline_category[node] - min(x[, node])
        x[, node] <- x[, node] - min(x[, node])

        if (check_fail_zero == FALSE) {
          check_fail_zero <- TRUE
          failed_zeroes <- c(node)
        } else {
          failed_zeroes <- c(failed_zeroes, node)
        }
      }

      check_range <- length(unique(x[, node]))
      if (check_range < 3) {
        stop(paste0(
          "The Blume-Capel is only available for variables with more than one category \n",
          "observed. There two or less categories observed for variable ",
          node,
          "."
        ))
      }
    }

    # Warn that maximum category value is large --------------------------------
    num_categories[node] <- max(x[, node])
    if (!is_ordinal[node] & num_categories[node] > 10) {
      warning(
        "Blume-Capel variable ", node, " has ", num_categories[node], " categories. ",
        "This may slow computation. Empty categories are not collapsed.",
        call. = FALSE
      )
    }

    # Check to see if not all responses are in one category --------------------
    if (num_categories[node] == 0) {
      stop(paste0(
        "Only one value [",
        unq_vls,
        "] was observed for variable ",
        node,
        "."
      ))
    }
  }

  if (check_fail_zero == TRUE && isTRUE(getOption("bgms.verbose", TRUE))) {
    nodes_str <- paste(failed_zeroes, collapse = ", ")
    message(
      "Variable", if (length(failed_zeroes) > 1) "s" else "", " ", nodes_str,
      " recoded to start at 0 (baseline categor",
      if (length(failed_zeroes) > 1) "ies" else "y", " adjusted)."
    )
  }

  list(
    x                 = x,
    num_categories    = num_categories,
    baseline_category = baseline_category
  )
}


# ------------------------------------------------------------------------------
# collapse_categories_across_groups
# ------------------------------------------------------------------------------
#
# For bgmCompare: collapses ordinal categories that are not observed in
# *all* groups, then renumbers the remaining categories contiguously
# (0-based). Blume-Capel variables are left unchanged.
#
# Called immediately after reformat_ordinal_data() in the compare path.
#
# Extracted from the group-aware ordinal recoding loop in
# compare_reformat_data() (data_utils.R) as part of Phase A.6b.
#
# @param x  Numeric matrix: data already recoded by reformat_ordinal_data().
# @param group  Integer vector of length nrow(x): group membership (1:K).
# @param is_ordinal  Logical vector of length ncol(x).
# @param num_categories  Integer vector from reformat_ordinal_data().
# @param baseline_category  Integer vector from reformat_ordinal_data().
#
# Returns:
#   list(x, num_categories, baseline_category)
# ------------------------------------------------------------------------------
collapse_categories_across_groups <- function(x,
                                              group,
                                              is_ordinal,
                                              num_categories,
                                              baseline_category) {
  num_variables <- ncol(x)
  num_groups <- max(group)

  for (node in seq_len(num_variables)) {
    if (!is_ordinal[node]) next  # BC variables: no group collapsing

    unq_vls <- sort(unique(x[, node]))
    n_unique <- length(unq_vls)

    # Build observed_scores matrix: which categories appear in which groups
    observed_scores <- matrix(NA, nrow = n_unique, ncol = num_groups)
    for (i in seq_along(unq_vls)) {
      for (g in seq_len(num_groups)) {
        observed_scores[i, g] <-
          as.integer(any(x[group == g, node] == unq_vls[i]))
      }
    }

    # Recode: keep only categories observed in ALL groups
    original <- x[, node]
    cntr <- -1L
    for (i in seq_along(unq_vls)) {
      if (sum(observed_scores[i, ]) == num_groups) {
        cntr <- cntr + 1L
      }
      x[original == unq_vls[i], node] <- max(0L, cntr)
    }

    num_categories[node] <- max(x[, node])

    # After collapsing, check that at least two categories remain
    if (num_categories[node] == 0) {
      stop(paste0("Only one value was observed for variable ", node, "."))
    }
  }

  list(
    x                 = x,
    num_categories    = num_categories,
    baseline_category = baseline_category
  )
}
