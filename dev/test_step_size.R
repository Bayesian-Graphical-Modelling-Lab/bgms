devtools::load_all()

set.seed(42)
n = 500L
x = matrix(sample(0:1, n * 2, replace = TRUE), n, 2)
y = matrix(rnorm(n * 2), n, 2)

input_list = list(
  discrete_observations = x,
  continuous_observations = y,
  num_categories = as.integer(c(1, 1)),
  is_ordinal_variable = as.integer(c(1, 1)),
  baseline_category = as.integer(c(0, 0)),
  main_alpha = 0.5,
  main_beta = 0.5,
  pairwise_scale = 2.5,
  pseudolikelihood = "full"
)

# p=2 binary: mux=2, Kxx=1, muy=2, Kxy=4 => dim=9
theta0 = rep(0, 9)

cat("--- Gradient test ---\n")
r = test_mixed_gradient(input_list, theta0)
cat("logp:", r$logp, "\n")
cat("grad length:", length(r$gradient), "\n")
cat("grad:", round(r$gradient, 2), "\n")
cat("grad norm:", round(sqrt(sum(r$gradient^2)), 2), "\n")

cat("\n--- Step-size heuristic ---\n")
r2 = test_mixed_step_size(input_list, theta0)
cat("Step size:", r2$step_size, "\n")
cat("Elapsed ms:", r2$elapsed_ms, "\n")
cat("logp at theta:", r2$logp_at_theta, "\n")
cat("Grad norm:", round(r2$grad_norm, 2), "\n")
