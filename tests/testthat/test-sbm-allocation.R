# --------------------------------------------------------------------------- #
# Exactness of the MFM-SBM allocation sweep (hierarchical graph prior).
#
# The sweep must leave p(z, theta | G) invariant. This guards the new-cluster
# parameter draw: a node that opens a new cluster (weighted by the integrated
# marginal) must draw the new between-cluster row from its conditional
# posterior given that node's edges, not from the prior (Neal 2000, Alg. 2).
# The prior draw over-weights one block on two triangles joined by an edge
# (P(T = 1) ~0.35 vs exact ~0.30) and is off by up to 0.08 in P(T = t) on
# graphs that disagree with an informative block prior.
#
# Checks, all against exact enumeration over partitions:
#   T0  the plain update on one fixed graph (test_sbm_allocation_chain);
#   T1  the same across graphs, block priors and count priors at p = 7; the
#       prior-only chain against the joint prior of (G, z) at p = 5; and GGM
#       data with large n, averaged over the graphs the chain visits;
#   T2  every route that uses the plain update (ordinal, binary, mixed,
#       bgmCompare), through the Rao-Blackwell identity
#       E[1{T = t} | X] = E[p(T = t | G) | X], valid because p(z | G, X) =
#       p(z | G) under the hierarchical specification.
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

# Node pairs in row-major upper-triangle order (1-2, 1-3, ..., (p-1)-p).
upper_pairs = function(p) {
  pr = which(upper.tri(diag(p)), arr.ind = TRUE)
  pr[order(pr[, 1L], pr[, 2L]), , drop = FALSE]
}

# log p(G | z) with block probabilities integrated out under Beta(a, b) within
# and Beta(a_between, b_between) between clusters, for every row of `graphs`
# (0/1, one column per pair in upper_pairs() order).
log_graph_given_partition = function(graphs, z, a, b, a_between, b_between) {
  pr = upper_pairs(length(z))
  k = max(z)
  lo = pmin(z[pr[, 1L]], z[pr[, 2L]])
  hi = pmax(z[pr[, 1L]], z[pr[, 2L]])
  cell = (lo - 1L) * k + hi
  cells = sort(unique(cell))
  onehot = outer(cell, cells, "==") * 1
  edges = graphs %*% onehot
  dyads = colSums(onehot)
  within = (cells - 1L) %/% k + 1L == (cells - 1L) %% k + 1L
  ar = ifelse(within, a, a_between)
  br = ifelse(within, b, b_between)
  n = nrow(edges)
  rowSums(lbeta(
    matrix(ar, n, length(cells), byrow = TRUE) + edges,
    matrix(br + dyads, n, length(cells), byrow = TRUE) - edges
  )) - sum(lbeta(ar, br))
}

# log p(z) for a partition in canonical labels: V_p(T) prod_c eta^(|c|).
log_partition_prior = function(z, log_v, eta) {
  log_v[max(z)] + sum(lgamma(tabulate(z) + eta) - lgamma(eta))
}

# Exact p(z | G) over all partitions.
exact_partition_posterior = function(g, partitions, log_v, eta = 1,
                                     a = 1, b = 1,
                                     a_between = 1, b_between = 1) {
  graphs = matrix(g[upper_pairs(nrow(g))], nrow = 1L)
  log_mass = vapply(partitions, function(z) {
    log_partition_prior(z, log_v, eta) +
      log_graph_given_partition(graphs, z, a, b, a_between, b_between)
  }, numeric(1))
  exp(log_mass - max(log_mass)) / sum(exp(log_mass - max(log_mass)))
}

# Batch-means Monte Carlo standard error of each column mean of x.
batch_mcse = function(x, n_batches = 100L) {
  x = as.matrix(x)
  size = floor(nrow(x) / n_batches)
  x = x[seq_len(size * n_batches), , drop = FALSE]
  means = rowsum(x, rep(seq_len(n_batches), each = size)) / size
  apply(means, 2, sd) / sqrt(n_batches)
}

num_clusters = function(allocations) {
  apply(allocations, 1, function(z) length(unique(z)))
}

# Sampled allocations against an exact partition law: partition-level total
# variation, and P(T = t) with z-scores in units of its MCSE.
compare_with_exact = function(allocations, partitions, prob) {
  p = length(partitions[[1L]])
  keys = vapply(partitions, paste, character(1), collapse = ",")
  sampled_keys = apply(allocations, 1, function(z) {
    paste(match(z, unique(z)), collapse = ",")
  })
  idx = match(sampled_keys, keys)
  t_of = vapply(partitions, max, integer(1))
  t_exact = vapply(seq_len(p), function(t) sum(prob[t_of == t]), numeric(1))
  ind = outer(num_clusters(allocations), seq_len(p), "==") * 1
  t_sampled = colMeans(ind)
  mcse = batch_mcse(ind)
  keep = t_exact > 1e-4 & mcse > 0
  list(
    tv_partition = 0.5 * sum(abs(tabulate(idx, length(keys)) / length(idx) - prob)),
    t_exact = t_exact,
    t_sampled = t_sampled,
    max_abs_z = max(abs(t_sampled[keep] - t_exact[keep]) / mcse[keep])
  )
}

skip_unless_slow_sbm = function() {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the SBM allocation exactness checks"
  )
}

two_triangles_bridge = function() {
  g = matrix(0L, 6, 6)
  g[rbind(c(1, 2), c(1, 3), c(2, 3), c(4, 5), c(4, 6), c(5, 6), c(3, 4))] = 1L
  g + t(g)
}


# --- T0 --------------------------------------------------------------------- #

test_that("sbm allocation update returns the exact partition posterior on a fixed graph", {
  g = two_triangles_bridge()
  partitions = enumerate_partitions(6)
  exact = exact_partition_posterior(g, partitions, mfm_log_v(6))

  out = test_sbm_allocation_chain(g, 80500L, 1, 1, 1, 1, 1, 1, seed = 20261101L)
  cmp = compare_with_exact(out$allocations[-seq_len(500L), ], partitions, exact)

  # Over 50 seeds the conditional draw gave TV <= 0.023 and |dP(T = t)| <=
  # 0.014; the prior draw gave TV >= 0.095 and |dP(T = t)| >= 0.039.
  expect_lt(cmp$tv_partition, 0.04)
  expect_lt(max(abs(cmp$t_sampled - cmp$t_exact)), 0.02)
})


# --- T1 --------------------------------------------------------------------- #

test_that("sbm allocation update is exact across graphs, block priors and count priors", {
  skip_unless_slow_sbm()

  p = 7
  pairs = upper_pairs(p)
  as_graph = function(keep) {
    g = matrix(0L, p, p)
    g[pairs[keep, , drop = FALSE]] = 1L
    g + t(g)
  }
  in_block = function(i, j) (i <= 4 & j <= 4) | (i >= 5 & j >= 5)
  graphs = list(
    two_blocks = as_graph(in_block(pairs[, 1], pairs[, 2]) |
      (pairs[, 1] == 4 & pairs[, 2] == 5)),
    one_block = as_graph((pairs[, 2] - pairs[, 1]) %% 7 %in% c(1, 2, 5, 6)),
    empty = as_graph(rep(FALSE, nrow(pairs))),
    near_complete = as_graph(!(pairs[, 1] == 1 & pairs[, 2] == 7) &
      !(pairs[, 1] == 3 & pairs[, 2] == 5))
  )
  block_priors = list(
    default = c(1, 1, 1, 1),
    weak_asymmetric = c(2, 1, 1, 2),
    strong = c(8, 1, 1, 8)
  )
  count_priors = list(c(lambda = 1, eta = 1), c(lambda = 3, eta = 0.5))
  partitions = enumerate_partitions(p)

  seed = 20261101L
  for(gn in names(graphs)) {
    for(bn in names(block_priors)) {
      for(cp in count_priors) {
        bp = block_priors[[bn]]
        exact = exact_partition_posterior(
          graphs[[gn]], partitions, mfm_log_v(p, cp[["eta"]], cp[["lambda"]]),
          cp[["eta"]], bp[1], bp[2], bp[3], bp[4]
        )
        seed = seed + 1L
        out = test_sbm_allocation_chain(
          graphs[[gn]], 101000L, bp[1], bp[2], bp[3], bp[4],
          cp[["eta"]], cp[["lambda"]],
          seed = seed
        )
        cmp = compare_with_exact(out$allocations[-seq_len(1000L), ], partitions, exact)
        info = sprintf("%s, %s, lambda = %g, eta = %g", gn, bn, cp[["lambda"]], cp[["eta"]])
        expect_lt(cmp$max_abs_z, 4.5, label = info)
        expect_lt(cmp$tv_partition, 0.05, label = info)
      }
    }
  }
})

test_that("prior-only sbm chain returns the joint prior of graph and partition", {
  skip_unless_slow_sbm()

  # Alternating the update with an exact draw G | z, theta targets the joint
  # p(G, z, theta) of the hierarchical specification. Enumerating all 52
  # partitions and 1024 graphs at p = 5 gives the exact law of the number of
  # clusters T, the number of edges E and the number of within-cluster edges
  # W, which ties the graph to the partition.
  p = 5
  pairs = upper_pairs(p)
  graphs = t(vapply(0:1023, function(g) {
    as.integer(intToBits(g))[seq_len(nrow(pairs))]
  }, integer(nrow(pairs))))
  partitions = enumerate_partitions(p)
  stat_table = function(z, g) {
    same = z[, pairs[, 1]] == z[, pairs[, 2]]
    cbind(T = num_clusters(z), E = rowSums(g), W = rowSums(g * same))
  }
  z_all = do.call(rbind, partitions)[rep(seq_along(partitions), nrow(graphs)), ]
  g_all = graphs[rep(seq_len(nrow(graphs)), each = length(partitions)), ]
  stats_all = stat_table(z_all, g_all)

  configs = list(
    c(a = 1, b = 1, ab = 1, bb = 1, lambda = 1, eta = 1),
    c(a = 2, b = 1, ab = 1, bb = 2, lambda = 3, eta = 0.5),
    c(a = 8, b = 1, ab = 1, bb = 8, lambda = 1, eta = 1)
  )
  for(i in seq_along(configs)) {
    cf = configs[[i]]
    log_v = mfm_log_v(p, cf[["eta"]], cf[["lambda"]])
    log_joint = vapply(partitions, function(z) {
      log_partition_prior(z, log_v, cf[["eta"]]) +
        log_graph_given_partition(graphs, z, cf[["a"]], cf[["b"]], cf[["ab"]], cf[["bb"]])
    }, numeric(nrow(graphs)))
    weight = as.numeric(t(exp(log_joint)))
    expect_equal(sum(weight), 1, tolerance = 1e-10)

    out = test_sbm_allocation_chain(
      matrix(0L, p, p), 201000L, cf[["a"]], cf[["b"]], cf[["ab"]], cf[["bb"]],
      cf[["eta"]], cf[["lambda"]],
      seed = 20261200L + i, resample_graph = TRUE
    )
    keep = -seq_len(1000L)
    sampled = stat_table(out$allocations[keep, ], out$graphs[keep, ])
    for(s in colnames(sampled)) {
      levels = sort(unique(stats_all[, s]))
      exact = vapply(levels, function(v) sum(weight[stats_all[, s] == v]), numeric(1))
      ind = outer(sampled[, s], levels, "==") * 1
      mcse = batch_mcse(ind, 50L)
      use = exact > 1e-4 & mcse > 0
      z = (colMeans(ind)[use] - exact[use]) / mcse[use]
      expect_lt(max(abs(z)), 4.5, label = sprintf("config %d, statistic %s", i, s))
    }
  }
})

test_that("sbm allocation sweep returns the exact partition posterior", {
  skip_unless_slow_sbm()

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
  pairs = upper_pairs(p)

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
  sampled = tabulate(num_clusters(allocations), nbins = p) / nrow(allocations)

  expect_lt(max(abs(sampled - exact)), 0.02)
})


# --- T2 --------------------------------------------------------------------- #

# Rao-Blackwell check on a fitted model: sampled P(T = t) against the average
# of the exact p(T = t | G) over the graphs the chains visit. `pair_names`
# label the indicator columns ("A-B"); `keep` selects the columns that carry
# the SBM prior.
expect_sbm_rao_blackwell = function(fit, pair_names, keep, node_names) {
  p = length(node_names)
  partitions = enumerate_partitions(p)
  t_of = vapply(partitions, max, integer(1))
  log_v = mfm_log_v(p)
  ends = do.call(rbind, strsplit(pair_names[keep], "-", fixed = TRUE))
  idx = cbind(match(ends[, 1], node_names), match(ends[, 2], node_names))
  expect_false(anyNA(idx))

  cache = list()
  diffs = list()
  rb_sum = 0
  for(c in seq_along(fit$raw_samples$allocations)) {
    ind = fit$raw_samples$indicator[[c]][, keep, drop = FALSE]
    keys = apply(ind, 1, paste, collapse = "")
    for(k in setdiff(unique(keys), names(cache))) {
      g = matrix(0, p, p)
      g[idx] = ind[match(k, keys), ]
      g = pmax(g, t(g))
      prob = exact_partition_posterior(g, partitions, log_v)
      cache[[k]] = vapply(seq_len(p), function(t) sum(prob[t_of == t]), numeric(1))
    }
    rb = do.call(rbind, cache[keys])
    sampled = outer(num_clusters(fit$raw_samples$allocations[[c]]), seq_len(p), "==") * 1
    diffs[[c]] = sampled - rb
    rb_sum = rb_sum + colSums(rb)
  }
  d = do.call(rbind, diffs)
  # Chains are pooled: batch means within each chain, combined.
  mcse = sqrt(Reduce(`+`, lapply(diffs, function(x) batch_mcse(x, 50L)^2))) /
    length(diffs)
  # A count with negligible mass is never visited, so its difference series
  # is nearly constant and the MCSE degenerates; score counts with real mass.
  use = mcse > 0 & rb_sum / nrow(d) > 1e-3
  expect_lt(max(abs(colMeans(d)[use]) / mcse[use]), 4.5)
  expect_lt(max(abs(colMeans(d))), 0.03)
}

sbm_route_blocks = function(within, bridge) {
  w = matrix(0, 6, 6)
  w[1:3, 1:3] = within
  w[4:6, 4:6] = within
  w[3, 4] = w[4, 3] = bridge
  diag(w) = 0
  w
}

# The route seeds below are those of certification/sbm_allocation_cert.R
# (20261101 + route index: ordinal, binary, ggm, mixed, compare), so these
# tests reproduce certified runs; GGM is checked by the T1 sweep test above.
# With the fix reverted, the certified ordinal, mixed and compare runs gave
# max |z| of 6.6 to 8.7, but the binary run only 4.2, below the 4.5 threshold:
# the binary case is a smoke test of the route, and the power comes from the
# other routes and from the T1 grid.

test_that("ordinal, binary and mixed bgm() fits sample the exact partition posterior", {
  skip_on_cran()
  skip_unless_certification()

  vn = paste0("V", 1:6)
  ordinal = simulate_mrf(1000, 6, 2, sbm_route_blocks(0.5, 0.25),
    matrix(-0.5, 6, 2),
    iter = 1000, seed = 101
  )
  binary = simulate_mrf(1500, 6, 1, sbm_route_blocks(0.8, 0.4),
    matrix(-1, 6, 1),
    iter = 1000, seed = 102
  )
  omega = diag(6)
  omega[sbm_route_blocks(1, 1) != 0] = -0.3
  set.seed(104)
  mixed = as.data.frame(matrix(rnorm(2000 * 6), 2000, 6) %*% chol(solve(omega)))
  for(j in c(1, 2, 5)) {
    mixed[[j]] = findInterval(mixed[[j]], quantile(mixed[[j]], c(1, 2) / 3))
  }
  mixed_type = rep("continuous", 6)
  mixed_type[c(1, 2, 5)] = "ordinal"

  cases = list(
    list(x = ordinal, type = "ordinal"),
    list(x = binary, type = "ordinal"),
    list(x = mixed, type = mixed_type)
  )
  for(i in seq_along(cases)) {
    x = cases[[i]]$x
    colnames(x) = vn
    fit = bgm(x,
      variable_type = cases[[i]]$type, edge_prior = sbm_prior(),
      iter = 20000, warmup = 2000, chains = 2, cores = 2,
      seed = c(20261102L, 20261103L, 20261105L)[i],
      display_progress = "none", verbose = FALSE
    )
    nm = fit$raw_samples$parameter_names$indicator
    expect_sbm_rao_blackwell(fit, nm, rep(TRUE, length(nm)), vn)
  }
})

test_that("bgmCompare() with an sbm difference prior samples the exact partition posterior", {
  skip_on_cran()
  skip_unless_certification()

  vn = paste0("V", 1:6)
  w1 = sbm_route_blocks(0.5, 0.25)
  w2 = w1
  w2[1:3, 1:3] = 0.05
  diag(w2) = 0
  # Thresholds that keep every category observed in both groups.
  main = matrix(c(-1.5, -4), 6, 2, byrow = TRUE)
  x = rbind(
    simulate_mrf(1000, 6, 2, w1, main, iter = 1000, seed = 105),
    simulate_mrf(1000, 6, 2, w2, main, iter = 1000, seed = 106)
  )
  colnames(x) = vn
  fit = bgmCompare(x,
    group_indicator = rep(1:2, each = 1000), variable_type = "ordinal",
    difference_prior = sbm_prior(), iter = 20000, warmup = 2000,
    chains = 2, cores = 2, seed = 20261106L,
    display_progress = "none", verbose = FALSE
  )
  nm = fit$raw_samples$parameter_names$indicator
  keep = grepl(" (pairwise)", nm, fixed = TRUE)
  expect_sbm_rao_blackwell(fit, sub(" (pairwise)", "", nm, fixed = TRUE), keep, vn)
})
