# --------------------------------------------------------------------------- #
# Exactness of the MFM-SBM allocation sweep (hierarchical graph prior).
#
# GGM data with large n pin the graph down, and under the hierarchical graph
# prior p(z | G, Y) = p(z | G). The sampled number of blocks must therefore
# match the exact partition posterior, enumerated over all partitions and
# averaged over the graphs the chain visits. This guards the new-cluster
# parameter draw: a node that opens a new cluster (weighted by the integrated
# marginal) must draw the new between-cluster row from its conditional
# posterior given that node's edges, not from the prior (Neal 2000, Alg. 2).
# The prior draw over-weights one block here (P(T = 1) ~0.36 vs exact ~0.33).
# --------------------------------------------------------------------------- #

enumerate_partitions = function(p) {
  out = list()
  extend = function(z) {
    if(length(z) == p) {
      out[[length(out) + 1L]] <<- z
    } else {
      for(j in seq_len(max(z) + 1L)) extend(c(z, j))
    }
  }
  extend(1L)
  out
}

# log V_p(t) for K - 1 ~ Poisson(lambda) components, symmetric Dirichlet(eta).
mfm_log_v = function(p, eta = 1, lambda = 1, k_max = 500L) {
  vapply(seq_len(p), function(t) {
    k = seq.int(t, k_max)
    terms = dpois(k - 1, lambda, log = TRUE) + lgamma(k + 1) -
      lgamma(k - t + 1) + lgamma(k * eta) - lgamma(k * eta + p)
    m = max(terms)
    m + log(sum(exp(terms - m)))
  }, numeric(1))
}

# Exact p(z | G) over all partitions; block probabilities integrated out under
# Beta(a, b) within and Beta(a_between, b_between) between clusters.
exact_partition_posterior = function(g, partitions, log_v, eta = 1,
                                     a = 1, b = 1,
                                     a_between = 1, b_between = 1) {
  log_mass = vapply(partitions, function(z) {
    out = log_v[max(z)] + sum(lgamma(tabulate(z) + eta) - lgamma(eta))
    for(r in seq_len(max(z))) {
      for(s in seq.int(r, max(z))) {
        use = upper.tri(g) & (outer(z == r, z == s) | outer(z == s, z == r))
        edges = sum(g[use])
        ar = if(r == s) a else a_between
        br = if(r == s) b else b_between
        out = out + lbeta(ar + edges, br + sum(use) - edges) - lbeta(ar, br)
      }
    }
    out
  }, numeric(1))
  exp(log_mass - max(log_mass)) / sum(exp(log_mass - max(log_mass)))
}

test_that("sbm allocation sweep returns the exact partition posterior", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the SBM allocation exactness check"
  )

  # Two triangles joined by one edge; n = 3000 leaves the graph nearly fixed.
  p = 6
  edges = rbind(c(1, 2), c(1, 3), c(2, 3), c(4, 5), c(4, 6), c(5, 6), c(3, 4))
  omega = diag(p)
  omega[edges] = -0.3
  omega[edges[, 2:1]] = -0.3
  set.seed(20261023)
  x = matrix(rnorm(3000 * p), 3000, p) %*% chol(solve(omega))

  fit = bgm(x,
    variable_type = "continuous", edge_prior = sbm_prior(),
    iter = 100000, warmup = 2000, chains = 2, cores = 2, seed = 20261023,
    display_progress = "none", verbose = FALSE
  )

  partitions = enumerate_partitions(p)
  num_blocks = vapply(partitions, max, integer(1))
  log_v = mfm_log_v(p)
  pairs = which(upper.tri(diag(p)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]

  indicator = do.call(rbind, fit$raw_samples$indicator)
  allocations = do.call(rbind, fit$raw_samples$allocations)
  graph_key = apply(indicator, 1, paste, collapse = "")
  exact = numeric(p)
  for(k in unique(graph_key)) {
    g = matrix(0, p, p)
    g[pairs] = indicator[match(k, graph_key), ]
    g = g + t(g)
    prob = exact_partition_posterior(g, partitions, log_v)
    exact = exact + mean(graph_key == k) *
      vapply(seq_len(p), function(t) sum(prob[num_blocks == t]), numeric(1))
  }
  sampled = tabulate(
    apply(allocations, 1, function(z) length(unique(z))), nbins = p
  ) / nrow(allocations)

  expect_lt(max(abs(sampled - exact)), 0.02)
})
