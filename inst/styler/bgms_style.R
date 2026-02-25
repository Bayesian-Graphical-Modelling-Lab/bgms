bgms_style <- function() {
  style <- styler::tidyverse_style()
  # style$space$spacing_before_comments()
  # "if (condition)" -> "if(condition)"
  style$space$add_space_after_for_if_while <- function(pd_flat) {
    paren_after <- pd_flat$token %in% c("FOR", "IF", "WHILE")
    if (!any(paren_after)) {
      return(pd_flat)
    }
    pd_flat$spaces[paren_after & (pd_flat$newlines == 0L)] <- 0L
    pd_flat
  }

  # "a <- 1" -> "a = 1"
  #
  # Preserves <- inside function-call arguments, where = would be interpreted

  # as a named argument (e.g. expect_message(result <- foo(), "pattern")).
  # Detection: the LHS expr's child node has token_before == "'('" when the
  # assignment sits inside a call's argument list.
  style$token$force_assignment_op <- function(pd) {
    to_replace <- pd$token %in% c("EQ_ASSIGN", "LEFT_ASSIGN")
    if(any(to_replace & pd$token == "LEFT_ASSIGN")) {
      for(i in which(to_replace & pd$token == "LEFT_ASSIGN")) {
        lhs_row <- i - 1
        if(lhs_row >= 1 && pd$token[lhs_row] == "expr") {
          child <- pd$child[[lhs_row]]
          if(!is.null(child) && nrow(child) > 0) {
            tb <- child$token_before[1]
            if(!is.na(tb) && tb == "'('") {
              to_replace[i] <- FALSE
            }
          }
        }
      }
    }
    pd$token[to_replace] <- "EQ_ASSIGN"
    pd$text[to_replace] <- "="
    pd
  }
  style$transformers_drop$token$force_assignment_op <-
    c("EQ_ASSIGN", "LEFT_ASSIGN")
  style
}

# run the styler locally via
# note that this cannot be undone! make sure to commit your changes beforehand
# source("inst/styler/bgms_style.R")
# styler::style_pkg(
#   style = bgms_style,
#   exclude_files = c("R/RcppExports\\.R", "R/cpp11\\.R", "R/import-standalone.*\\.R", "\\.Rprofile\\.R")
# )
