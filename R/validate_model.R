# ==============================================================================
# Model validation functions
# ==============================================================================
#
# Extracted from check_model() and check_compare_model() as part of the
# R scaffolding refactor (dev/scaffolding/plan.md, Phase A).
#
# Each validator is a pure function: input -> validated output (or error).
# Validators do NOT read from parent environments or use hasArg().
# ==============================================================================

# ------------------------------------------------------------------------------
# validate_variable_types
# ------------------------------------------------------------------------------
#
# Parses and validates the variable_type argument for both bgm() and
# bgmCompare(). Returns a canonical list with:
#   $variable_type  - character vector of length num_variables
#   $variable_bool  - logical vector (TRUE = ordinal, FALSE = blume-capel)
#   $is_continuous  - scalar logical
#
# @param variable_type Character: single string or vector of variable types.
# @param num_variables Integer: number of columns in the data matrix.
# @param allow_continuous Logical: whether "continuous" is a valid type.
#   TRUE for bgm(), FALSE for bgmCompare() (until GGM-Compare is added).
# @param caller Character: name of the calling function, used in error messages.
#
# Replaces:
#   - check_model()         lines 53-133  (function_input_utils.R)
#   - check_compare_model() lines 405-468 (function_input_utils.R)
# ------------------------------------------------------------------------------
validate_variable_types <- function(variable_type,
                                    num_variables,
                                    allow_continuous = TRUE,
                                    caller = "bgm") {
  valid_choices <- if (allow_continuous) {
    c("ordinal", "blume-capel", "continuous")
  } else {
    c("ordinal", "blume-capel")
  }

  supported_str <- paste(valid_choices, collapse = ", ")

  is_continuous <- FALSE

  if (length(variable_type) == 1) {
    # --- Single string: replicate to all variables ---
    variable_input <- variable_type
    variable_type <- try(
      match.arg(arg = variable_type, choices = valid_choices),
      silent = TRUE
    )
    if (inherits(variable_type, what = "try-error")) {
      stop(paste0(
        "The ", caller, " function supports variables of type ", supported_str,
        ", but not of type ", variable_input, "."
      ))
    }

    if (variable_type == "continuous") {
      is_continuous <- TRUE
      variable_bool <- rep(TRUE, num_variables)
      variable_type <- rep("continuous", num_variables)
    } else {
      variable_bool <- (variable_type == "ordinal")
      variable_bool <- rep(variable_bool, num_variables)
      variable_type <- rep(variable_type, num_variables)
    }
  } else {
    # --- Vector of types: validate each element ---
    if (length(variable_type) != num_variables) {
      stop(paste0(
        "The variable type vector variable_type should be either a single character\n",
        "string or a vector of character strings of length p."
      ))
    }

    has_continuous <- any(variable_type == "continuous")
    if (has_continuous && !all(variable_type == "continuous")) {
      stop(paste0(
        "When using continuous variables, all variables must be of type ",
        "'continuous'. Mixtures of continuous and ordinal/blume-capel ",
        "variables are not supported."
      ))
    }

    if (has_continuous) {
      if (!allow_continuous) {
        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type continuous."
        ))
      }
      is_continuous <- TRUE
      variable_bool <- rep(TRUE, num_variables)
    } else {
      variable_input <- unique(variable_type)
      non_continuous_choices <- c("ordinal", "blume-capel")

      variable_type_checked <- try(
        match.arg(
          arg = variable_type,
          choices = non_continuous_choices,
          several.ok = TRUE
        ),
        silent = TRUE
      )

      if (inherits(variable_type_checked, what = "try-error")) {
        # Identify which types are invalid
        num_types <- sapply(variable_input, function(type) {
          tmp <- try(
            match.arg(arg = type, choices = non_continuous_choices),
            silent = TRUE
          )
          inherits(tmp, what = "try-error")
        })

        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type ",
          paste0(variable_input[num_types], collapse = ", "), "."
        ))
      }

      # The match.arg with several.ok may have normalized partial matches
      variable_type <- variable_type_checked

      # Re-check length after match.arg (defensive)
      if (length(variable_type) != num_variables) {
        num_types <- sapply(variable_input, function(type) {
          tmp <- try(
            match.arg(arg = type, choices = non_continuous_choices),
            silent = TRUE
          )
          inherits(tmp, what = "try-error")
        })

        stop(paste0(
          "The ", caller, " function supports variables of type ", supported_str,
          ", but not of type ",
          paste0(variable_input[num_types], collapse = ", "), "."
        ))
      }

      variable_bool <- (variable_type == "ordinal")
    }
  }

  list(
    variable_type  = variable_type,
    variable_bool  = variable_bool,
    is_continuous  = is_continuous
  )
}


# ------------------------------------------------------------------------------
# validate_baseline_category
# ------------------------------------------------------------------------------
#
# Validates and normalizes the baseline_category argument for Blume-Capel
# variables. Shared by both bgm() and bgmCompare().
#
# @param baseline_category  The user-supplied baseline_category value.
# @param baseline_category_provided  Logical: whether the user actually
#   supplied the argument (i.e., `hasArg("baseline_category")` from the
#   calling scope). Needed because baseline_category has no default.
# @param x  Numeric matrix: the (validated) data.
# @param variable_bool  Logical vector: TRUE = ordinal, FALSE = blume-capel.
#
# Returns:
#   Integer vector of length ncol(x). For ordinal-only models, all zeros.
#
# Replaces:
#   - check_model()         lines 63-139  (function_input_utils.R)
#   - check_compare_model() lines 340-418 (function_input_utils.R)
# ------------------------------------------------------------------------------
validate_baseline_category <- function(baseline_category,
                                       baseline_category_provided,
                                       x,
                                       variable_bool) {
  num_variables <- ncol(x)

  # If all ordinal (no Blume-Capel), return zeros
  if (!any(!variable_bool)) {
    return(rep.int(0, times = num_variables))
  }

  # --- Blume-Capel variables present ---

  if (!baseline_category_provided) {
    stop("The argument baseline_category is required for Blume-Capel variables.")
  }

  if (length(baseline_category) != num_variables && length(baseline_category) != 1) {
    stop(paste0(
      "The argument baseline_category for the Blume-Capel model needs to be a \n",
      "single integer or a vector of integers of length p."
    ))
  }

  # Scalar: validate then replicate
  if (length(baseline_category) == 1) {
    integer_check <- try(as.integer(baseline_category), silent = TRUE)
    if (is.na(integer_check)) {
      stop(paste0(
        "The baseline_category argument for the Blume-Capel model contains either \n",
        "a missing value or a value that could not be forced into an integer value."
      ))
    }
    integer_check <- abs(baseline_category - round(baseline_category))
    if (integer_check > .Machine$double.eps) {
      stop("Reference category needs to an integer value or a vector of integers of length p.")
    }
    baseline_category <- rep.int(baseline_category, times = num_variables)
  }

  # Validate integer-ness for Blume-Capel variables
  blume_capel_variables <- which(!variable_bool)

  integer_check <- try(as.integer(baseline_category[blume_capel_variables]),
    silent = TRUE
  )
  if (anyNA(integer_check)) {
    stop(paste0(
      "The baseline_category argument for the Blume-Capel model contains either \n",
      "missing values or values that could not be forced into an integer value."
    ))
  }

  integer_check <- abs(baseline_category[blume_capel_variables] -
    round(baseline_category[blume_capel_variables]))

  if (any(integer_check > .Machine$double.eps)) {
    non_integers <- blume_capel_variables[integer_check > .Machine$double.eps]
    if (length(non_integers) > 1) {
      stop(paste0(
        "The entries in baseline_category for variables ",
        paste0(non_integers, collapse = ", "), " need to be integer."
      ))
    } else {
      stop(paste0(
        "The entry in baseline_category for variable ",
        non_integers, " needs to be an integer."
      ))
    }
  }

  # Validate within observed data range
  variable_lower <- apply(x, 2, min, na.rm = TRUE)
  variable_upper <- apply(x, 2, max, na.rm = TRUE)

  if (any(baseline_category < variable_lower) | any(baseline_category > variable_upper)) {
    out_of_range <- which(baseline_category < variable_lower | baseline_category > variable_upper)
    stop(paste0(
      "The Blume-Capel model assumes that the reference category is within the range \n",
      "of the observed category scores. This was not the case for variable(s) \n",
      paste0(out_of_range, collapse = ", "),
      "."
    ))
  }

  baseline_category
}
