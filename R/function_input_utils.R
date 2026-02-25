check_positive_integer = function(value, name) {
  if(!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value <= 0) {
    stop(sprintf("Parameter `%s` must be a positive integer. Got: %s", name, value))
  }
}

# Helper function for validating non-negative integers
check_non_negative_integer = function(value, name) {
  if(!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value < 0) {
    stop(sprintf("Parameter `%s` must be a non-negative integer. Got: %s", name, value))
  }
}

# Helper function for validating logical inputs
check_logical = function(value, name) {

  value = as.logical(value)
  if(is.na(value)) {
    stop(sprintf("Parameter `%s` must be TRUE or FALSE. Got: %s", name, value))
  }
  return(value)
}

# Helper function for validating optional seed parameter
# Returns a valid integer seed (generates one if NULL)
check_seed = function(seed) {
  if(is.null(seed)) {
    return(sample.int(.Machine$integer.max, 1L))
  }
  if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
    stop("Argument 'seed' must be a single non-negative integer.")
  }
  as.integer(seed)
}

progress_type_from_display_progress <- function(display_progress = c("per-chain", "total", "none")) {
  if(is.logical(display_progress) && length(display_progress) == 1) {
    if(is.na(display_progress)) {
      stop("The display_progress argument must be a single logical value, but not NA.")
    }
    display_progress = if(display_progress) "per-chain" else "none"
  } else {
    display_progress = match.arg(display_progress)
  }
  return(if(display_progress == "per-chain") 2L else if(display_progress == "total") 1L else 0L)
}
