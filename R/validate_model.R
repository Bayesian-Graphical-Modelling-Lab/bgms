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
