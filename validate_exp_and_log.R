
values <- matrix(c(
  -1e6, 1e6,
  -1e4, 1e4,
  -1e3, 1e3,
  -1e2, 1e2,
  -1e1, 1e1,
  -1e-3, 1e-3
), ncol = 2, byrow = TRUE)
size <- 1E6L

for (i in seq_len(nrow(values))) {
  from <- values[i, 1]
  to   <- values[i, 2]
  cat("i = ", i, ", Testing range: ", from, " to ", to, "\n")
  test_values <- seq(from, to, length.out = size)
  y_R <- exp(test_values)
  y_R2 <- bgms:::rcpp_ieee754_exp(test_values)
  result <- all.equal(y_R, y_R2)

  abs_error = max(abs(y_R - y_R2),         na.rm = TRUE)
  rel_error = max(abs((y_R - y_R2) / y_R), na.rm = TRUE)

  cat("i = ", i, ", All.equal = ", isTRUE(result),
      " max abs error:", sprintf("%.8g", abs_error),
      " max rel error:", sprintf("%.8g", rel_error),
      "\n")
}

values <- matrix(c(
  1e10, 1e30,
  1e04, 1e09,
  1e01, 1e03,
  1e-2, 1e01,
  1e-5, 1e-3,
  1e-20, 1e-6
), ncol = 2, byrow = TRUE)
size <- 1E6L

for (i in seq_len(nrow(values))) {
  from <- values[i, 1]
  to   <- values[i, 2]
  cat("i = ", i, ", Testing range: ", from, " to ", to, "\n")
  test_values <- seq(from, to, length.out = size)
  y_R <- log(test_values)
  y_R2 <- bgms:::rcpp_ieee754_log(test_values)
  result <- all.equal(y_R, y_R2)

  abs_error = max(abs(y_R - y_R2),         na.rm = TRUE)
  rel_error = max(abs((y_R - y_R2) / y_R), na.rm = TRUE)

  cat("i = ", i, ", All.equal = ", isTRUE(result),
      " max abs error:", sprintf("%.8g", abs_error),
      " max rel error:", sprintf("%.8g", rel_error),
      "\n")
}

