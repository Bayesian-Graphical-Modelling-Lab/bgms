# SBM allocation certification.
#
# Certifies that the plain (hierarchical) MFM-SBM update samples the partition
# law the model defines, against exact enumeration over partitions. It guards
# the new-cluster parameter draw: a node that opens a new cluster must draw the
# new between-cluster row from its conditional posterior given that node's
# edges (Neal 2000, Alg. 2), not from the prior, as bgms <= 0.2.0.0 did.
#
# Three arms, each scored as total variation (TV) and as z-scores in units of
# the batch-means Monte Carlo standard error (MCSE):
#
#   1. Fixed graphs. test_sbm_allocation_chain() on four p = 7 graphs (two
#      blocks, one block, empty, near complete) x three block priors (1,1,1,1),
#      (2,1,1,2), (8,1,1,8) x four (lambda, eta), plus three p = 6 graphs from
#      the audit, against the exact p(z | G) over all 877 (203) partitions.
#   2. Prior only. The same update alternated with an exact draw
#      G | z, theta ~ Bernoulli, and the production GGM hierarchical prior
#      chain (sample_ggm_prior(spec = "hierarchical")), against the exact
#      joint law of (G, z) at p = 5 (52 partitions x 1024 graphs). I.i.d.
#      ancestral draws check the enumeration itself.
#   3. Routes. Fitted ordinal, binary, GGM, mixed and bgmCompare models,
#      through the Rao-Blackwell identity E[1{T = t} | X] = E[p(T = t | G) | X],
#      valid because p(z | G, X) = p(z | G) under the hierarchical
#      specification. z uses the MCSE of d_s = 1{T_s = t} - p(T = t | G_s).
#      The fixed-build run also exercises summary(), print() and the
#      extractors on each fit.
#
# The precision_graph_prior = "joint" kernel (block_allocations_mfm_sbm_
# corrected) is not covered: it is a separate, approximate update.
#
# A control build with the fix reverted shows the arms have power:
#   git archive <commit> | tar -x -C ctrl
#   git show a12a36f -- src/priors/sbm_edge_prior.cpp | (cd ctrl && patch -R -p1)
#   R CMD INSTALL -l ctrl_lib ctrl
#
# Usage, from the repo root:
#   Rscript certification/sbm_allocation_cert.R --lib=<library> --label=fixed \
#     --tree=<git tree hash of the source> --parts=grid,prior,routes
# writes certification/sbm_allocation_cert_<label>.rds.

opt = function(name, default) {
  a = grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
  if(length(a) == 0L) default else sub(paste0("^--", name, "="), "", a[1])
}
lib = opt("lib", .libPaths()[1])
label = opt("label", "fixed")
tree = opt("tree", NA_character_)
parts = strsplit(opt("parts", "grid,prior,routes"), ",")[[1]]
out = opt("out", file.path("certification", paste0("sbm_allocation_cert_", label, ".rds")))

suppressMessages(library(bgms, lib.loc = lib))
stopifnot(normalizePath(dirname(find.package("bgms"))) == normalizePath(lib))


# ==============================================================================
# Exact enumeration
# ==============================================================================
# Hierarchy: B - 1 ~ Poisson(lambda); delta | B ~ Dir_B(eps); z_i ~ Cat(delta);
# theta_rr ~ Beta(a, b), theta_rs ~ Beta(ab, bb); G_ij ~ Bern(theta_{z_i z_j}).

enumerate_partitions = function(p) {
  rows = list(1L)
  for(k in seq_len(p - 1L)) {
    rows = unlist(lapply(rows, function(z) {
      lapply(seq_len(max(z) + 1L), function(j) c(z, j))
    }), recursive = FALSE)
  }
  do.call(rbind, rows)
}

canon_key = function(Z) {
  apply(Z, 1, function(z) paste(match(z, unique(z)), collapse = ""))
}

upper_pairs = function(p) {
  pr = which(upper.tri(diag(p)), arr.ind = TRUE)
  pr[order(pr[, 1L], pr[, 2L]), , drop = FALSE]
}

# log V_p(t), t = 1..p, with the omitted tail checked.
log_vp = function(p, eps, lambda) {
  b = seq_len(max(600L, ceiling(lambda + 60 * sqrt(lambda) + 100)))
  base = dpois(b - 1, lambda, log = TRUE) + lgamma(b * eps) - lgamma(b * eps + p)
  vapply(seq_len(p), function(t) {
    k = b[b >= t]
    w = base[b >= t] + lgamma(k + 1) - lgamma(k - t + 1)
    m = max(w)
    stopifnot(w[length(w)] - m < -40)
    m + log(sum(exp(w - m)))
  }, numeric(1))
}

log_partition_prior = function(P, eps, lambda) {
  lv = log_vp(ncol(P), eps, lambda)
  lv[apply(P, 1, max)] +
    apply(P, 1, function(z) sum(lgamma(tabulate(z) + eps) - lgamma(eps)))
}

# log p(G | z) for each partition (row of P); block probabilities integrated.
log_graph_given_partition = function(G, P, a, b, ab, bb) {
  diag(G) = 0
  apply(P, 1, function(z) {
    Zo = outer(z, seq_len(max(z)), "==") * 1
    E = crossprod(Zo, G %*% Zo)
    n = colSums(Zo)
    N = outer(n, n)
    diag(E) = diag(E) / 2
    diag(N) = n * (n - 1) / 2
    up = upper.tri(E)
    sum(lbeta(a + diag(E), b + diag(N) - diag(E)) - lbeta(a, b)) +
      sum(lbeta(ab + E[up], bb + N[up] - E[up]) - lbeta(ab, bb))
  })
}

normalize_log = function(lw) {
  w = exp(lw - max(lw))
  w / sum(w)
}

exact_partition_posterior = function(G, P, pr) {
  normalize_log(log_partition_prior(P, pr$eps, pr$lambda) +
                  log_graph_given_partition(G, P, pr$a, pr$b, pr$ab, pr$bb))
}

t_marginal = function(prob, P) {
  t = apply(P, 1, max)
  vapply(seq_len(ncol(P)), function(k) sum(prob[t == k]), numeric(1))
}


# ==============================================================================
# Monte Carlo summaries
# ==============================================================================

batch_mcse = function(x, n_batches = 100L) {
  x = as.matrix(x)
  size = floor(nrow(x) / n_batches)
  x = x[seq_len(size * n_batches), , drop = FALSE]
  means = rowsum(x, rep(seq_len(n_batches), each = size)) / size
  apply(means, 2, sd) / sqrt(n_batches)
}

tv = function(p, q) 0.5 * sum(abs(p - q))

compare_partitions = function(alloc, P, exact) {
  p = ncol(P)
  idx = match(canon_key(alloc), canon_key(P))
  stopifnot(!anyNA(idx))
  emp = tabulate(idx, nbins = nrow(P)) / length(idx)
  ind_t = outer(apply(alloc, 1, function(z) length(unique(z))), seq_len(p), "==") * 1
  t_emp = colMeans(ind_t)
  t_exact = t_marginal(exact, P)
  t_mcse = batch_mcse(ind_t)
  z_t = ifelse(t_mcse > 0, (t_emp - t_exact) / t_mcse, NA_real_)
  big = which(exact >= 0.01)
  z_p = (emp[big] - exact[big]) / batch_mcse(outer(idx, big, "==") * 1)
  list(t_exact = t_exact, t_emp = t_emp, t_mcse = t_mcse, z_t = z_t,
       tv_t = tv(t_emp, t_exact), tv_partition = tv(emp, exact),
       max_abs_z_t = max(abs(z_t[t_exact > 1e-4]), na.rm = TRUE),
       max_abs_z_partition = max(abs(z_p)))
}


# ==============================================================================
# Arm 1: fixed graphs
# ==============================================================================

arm_grid = function() {
  make_graph = function(p, keep) {
    pr = upper_pairs(p)
    G = matrix(0L, p, p)
    G[pr[keep, , drop = FALSE]] = 1L
    G + t(G)
  }
  pr7 = upper_pairs(7)
  in_block = (pr7[, 1] <= 4 & pr7[, 2] <= 4) | (pr7[, 1] >= 5 & pr7[, 2] >= 5)
  graphs = list(
    two_blocks = make_graph(7, in_block | (pr7[, 1] == 4 & pr7[, 2] == 5)),
    one_block = make_graph(7, (pr7[, 2] - pr7[, 1]) %% 7 %in% c(1, 2, 5, 6)),
    empty = make_graph(7, rep(FALSE, nrow(pr7))),
    near_complete = make_graph(7, !(pr7[, 1] == 1 & pr7[, 2] == 7) &
                                 !(pr7[, 1] == 3 & pr7[, 2] == 5))
  )
  pr6 = upper_pairs(6)
  e6 = function(edges) {
    G = matrix(0L, 6, 6)
    G[edges] = 1L
    pmax(G, t(G))
  }
  audit = list(
    two_triangles_bridge = e6(rbind(c(1, 2), c(1, 3), c(2, 3), c(4, 5),
                                    c(4, 6), c(5, 6), c(3, 4))),
    star = e6(cbind(1L, 2:6)),
    complete = e6(pr6)
  )
  block_priors = list(
    default = c(a = 1, b = 1, ab = 1, bb = 1),
    weak_asymmetric = c(a = 2, b = 1, ab = 1, bb = 2),
    strong = c(a = 8, b = 1, ab = 1, bb = 8)
  )
  count_priors = list(c(lambda = 1, eps = 1), c(lambda = 3, eps = 1),
                      c(lambda = 1, eps = 0.5), c(lambda = 0.5, eps = 2))
  cells = list()
  for(g in names(graphs)) for(bp in names(block_priors)) for(cp in count_priors) {
    cells[[length(cells) + 1L]] = list(g = g, G = graphs[[g]], bp = bp, cp = cp)
  }
  for(g in names(audit)) {
    cells[[length(cells) + 1L]] = list(g = g, G = audit[[g]], bp = "default",
                                       cp = c(lambda = 1, eps = 1))
  }
  parts = list(`6` = enumerate_partitions(6), `7` = enumerate_partitions(7))
  n_iter = 401000L
  burn = 1000L
  rows = lapply(seq_along(cells), function(i) {
    cl = cells[[i]]
    b = block_priors[[cl$bp]]
    P = parts[[as.character(nrow(cl$G))]]
    exact = exact_partition_posterior(cl$G, P, as.list(c(b, cl$cp)))
    seed = 20261101L + i
    r = bgms:::test_sbm_allocation_chain(cl$G, n_iter, b[["a"]], b[["b"]], b[["ab"]],
                                  b[["bb"]], cl$cp[["eps"]], cl$cp[["lambda"]], seed)
    cmp = compare_partitions(r$allocations[-seq_len(burn), , drop = FALSE], P, exact)
    cat(sprintf("grid %-21s %-15s lambda=%.1f eps=%.1f TV_T=%.4f TV_part=%.4f max|z_T|=%.1f\n",
                cl$g, cl$bp, cl$cp[["lambda"]], cl$cp[["eps"]], cmp$tv_t,
                cmp$tv_partition, cmp$max_abs_z_t))
    data.frame(graph = cl$g, p = nrow(cl$G), block_prior = cl$bp,
               lambda = cl$cp[["lambda"]], eps = cl$cp[["eps"]], seed = seed,
               n_draws = n_iter - burn, tv_t = cmp$tv_t,
               tv_partition = cmp$tv_partition, max_abs_z_t = cmp$max_abs_z_t,
               max_abs_z_partition = cmp$max_abs_z_partition,
               t_exact = I(list(cmp$t_exact)), t_sampled = I(list(cmp$t_emp)),
               z_t = I(list(cmp$z_t)))
  })
  do.call(rbind, rows)
}


# ==============================================================================
# Arm 2: prior only, joint (G, z) at p = 5
# ==============================================================================

arm_prior = function() {
  p = 5L
  P = enumerate_partitions(p)
  pairs = upper_pairs(p)
  np = nrow(pairs)
  bits = t(sapply(0:(2^np - 1), function(g) as.integer(intToBits(g))[seq_len(np)]))
  stats_of = function(z, g) {
    same = z[, pairs[, 1], drop = FALSE] == z[, pairs[, 2], drop = FALSE]
    cbind(T = apply(z, 1, function(r) length(unique(r))),
          E = rowSums(g), W = rowSums(g * same))
  }
  st_all = stats_of(P[rep(seq_len(nrow(P)), times = nrow(bits)), ],
                    bits[rep(seq_len(nrow(bits)), each = nrow(P)), ])
  exact_joint = function(pr) {
    lpc = log_partition_prior(P, pr$eps, pr$lambda)
    lg = sapply(seq_len(nrow(bits)), function(k) {
      G = matrix(0, p, p)
      G[pairs] = bits[k, ]
      log_graph_given_partition(G + t(G), P, pr$a, pr$b, pr$ab, pr$bb)
    })
    joint = exp(lpc + lg)
    stopifnot(abs(sum(joint) - 1) < 1e-10)
    list(partition = rowSums(joint), graph = colSums(joint), w = as.numeric(joint))
  }
  compare_stats = function(z, g, ex, n_batches) {
    st = stats_of(z, g)
    res = list()
    for(s in c("T", "E", "W")) {
      lv = sort(unique(st_all[, s]))
      exact = vapply(lv, function(v) sum(ex$w[st_all[, s] == v]), numeric(1))
      ind = outer(st[, s], lv, "==") * 1
      se = batch_mcse(ind, n_batches)
      use = se > 0 & exact > 1e-4
      res[[paste0("tv_", s)]] = tv(colMeans(ind), exact)
      res[[paste0("max_abs_z_", s)]] = max(abs((colMeans(ind)[use] - exact[use]) / se[use]))
    }
    key_all = paste(st_all[, "T"], st_all[, "E"])
    lev = sort(unique(key_all))
    exact_te = tapply(ex$w, factor(key_all, lev), sum)
    key = factor(paste(st[, "T"], st[, "E"]), lev)
    emp_te = tabulate(as.integer(key), nbins = length(lev)) / nrow(st)
    big = which(exact_te >= 0.01)
    res$tv_TE = tv(emp_te, exact_te)
    res$max_abs_z_TE = max(abs((emp_te[big] - exact_te[big]) /
                                 batch_mcse(outer(as.integer(key), big, "==") * 1, n_batches)))
    pk = match(canon_key(z), canon_key(P))
    res$tv_partition = tv(tabulate(pk, nrow(P)) / length(pk), ex$partition)
    gk = as.integer(g %*% (2^(seq_len(np) - 1))) + 1L
    res$tv_graph = tv(tabulate(gk, nrow(bits)) / length(gk), ex$graph)
    res
  }
  configs = list(
    default = list(a = 1, b = 1, ab = 1, bb = 1, lambda = 1, eps = 1),
    weak_asymmetric = list(a = 2, b = 1, ab = 1, bb = 2, lambda = 3, eps = 0.5),
    strong = list(a = 8, b = 1, ab = 1, bb = 8, lambda = 1, eps = 1)
  )
  rows = list()
  for(nm in names(configs)) {
    pr = configs[[nm]]
    ex = exact_joint(pr)
    seed = 7000L + match(nm, names(configs))
    ch = bgms:::test_sbm_allocation_chain(matrix(0L, p, p), 401000L, pr$a, pr$b, pr$ab,
                                   pr$bb, pr$eps, pr$lambda, seed, resample_graph = TRUE)
    keep = -seq_len(1000L)
    src = list(kernel = list(seed = seed, n = 400000L, r = compare_stats(
      ch$allocations[keep, ], ch$graphs[keep, ], ex, 100L)))
    if(label == "fixed") {
      set.seed(seed + 100L)
      n_anc = 100000L
      za = matrix(0L, n_anc, p)
      ga = matrix(0L, n_anc, np)
      for(s in seq_len(n_anc)) {
        zz = bgms:::ancestral_mfm_sbm_partition(p, pr$lambda, pr$eps)
        pp = bgms:::ancestral_sbm_pair_probabilities(zz, pr$a, pr$b, pr$ab, pr$bb)
        za[s, ] = zz
        ga[s, ] = as.integer(runif(np) < pp[pairs])
      }
      src$ancestral = list(seed = seed + 100L, n = n_anc,
                           r = compare_stats(za, ga, ex, 100L))
    }
    gd = sample_ggm_prior(
      p = p, n_samples = 400000L, n_warmup = 5000L, spec = "hierarchical",
      edge_prior = sbm_prior(alpha = pr$a, beta = pr$b, alpha_between = pr$ab,
                             beta_between = pr$bb, dirichlet_alpha = pr$eps,
                             lambda = pr$lambda),
      seed = seed + 200L, verbose = FALSE, zratio_diagnostics = FALSE
    )
    # Longer batches: the GGM chain mixes more slowly in the edge count.
    src$ggm = list(seed = seed + 200L, n = 400000L,
                   r = compare_stats(gd$allocations, gd$edge_indicators, ex, 50L))
    for(sn in names(src)) {
      r = src[[sn]]$r
      cat(sprintf("prior %-15s %-9s TV: T=%.4f E=%.4f W=%.4f TE=%.4f | max|z|: T=%.1f E=%.1f W=%.1f TE=%.1f\n",
                  nm, sn, r$tv_T, r$tv_E, r$tv_W, r$tv_TE, r$max_abs_z_T,
                  r$max_abs_z_E, r$max_abs_z_W, r$max_abs_z_TE))
      rows[[length(rows) + 1L]] = data.frame(config = nm, source = sn,
                                             seed = src[[sn]]$seed, n_draws = src[[sn]]$n, r)
    }
  }
  do.call(rbind, rows)
}


# ==============================================================================
# Arm 3: routes
# ==============================================================================

arm_routes = function() {
  p = 6L
  P = enumerate_partitions(p)
  pr = list(a = 1, b = 1, ab = 1, bb = 1, eps = 1, lambda = 1)
  vn = paste0("V", seq_len(p))
  blocks = function(within, bridge) {
    w = matrix(0, p, p)
    w[1:3, 1:3] = within
    w[4:6, 4:6] = within
    w[3, 4] = w[4, 3] = bridge
    diag(w) = 0
    w
  }
  omega = diag(p)
  omega[blocks(1, 1) != 0] = -0.3
  gen = list(
    ordinal = function() {
      list(x = simulate_mrf(1000, p, 2, blocks(0.5, 0.25), matrix(-0.5, p, 2),
                            iter = 1000, seed = 101), type = "ordinal")
    },
    binary = function() {
      list(x = simulate_mrf(1500, p, 1, blocks(0.8, 0.4), matrix(-1, p, 1),
                            iter = 1000, seed = 102), type = "ordinal")
    },
    ggm = function() {
      set.seed(103)
      list(x = matrix(rnorm(3000 * p), 3000, p) %*% chol(solve(omega)),
           type = "continuous")
    },
    mixed = function() {
      set.seed(104)
      x = as.data.frame(matrix(rnorm(2000 * p), 2000, p) %*% chol(solve(omega)))
      for(j in c(1, 2, 5)) x[[j]] = findInterval(x[[j]], quantile(x[[j]], c(1, 2) / 3))
      type = rep("continuous", p)
      type[c(1, 2, 5)] = "ordinal"
      list(x = x, type = type)
    },
    compare = function() {
      w2 = blocks(0.5, 0.25)
      w2[1:3, 1:3] = 0.05
      diag(w2) = 0
      main = matrix(c(-1.5, -4), p, 2, byrow = TRUE)
      x = rbind(simulate_mrf(1000, p, 2, blocks(0.5, 0.25), main, iter = 1000, seed = 105),
                simulate_mrf(1000, p, 2, w2, main, iter = 1000, seed = 106))
      list(x = x, type = "ordinal", group = rep(1:2, each = 1000))
    }
  )
  rb_check = function(ind_list, alloc_list, pidx) {
    cache = new.env()
    d_list = list()
    emp = 0
    rb = 0
    n_tot = 0
    for(c in seq_along(ind_list)) {
      ind = ind_list[[c]]
      keys = apply(ind, 1, paste, collapse = "")
      uk = unique(keys)
      tab = t(vapply(uk, function(k) {
        if(is.null(cache[[k]])) {
          G = matrix(0, p, p)
          G[pidx] = ind[match(k, keys), ]
          assign(k, t_marginal(exact_partition_posterior(pmax(G, t(G)), P, pr), P),
                 envir = cache)
        }
        cache[[k]]
      }, numeric(p)))
      pt = tab[match(keys, uk), , drop = FALSE]
      tt = outer(apply(alloc_list[[c]], 1, function(z) length(unique(z))),
                 seq_len(p), "==") * 1
      d_list[[c]] = tt - pt
      emp = emp + colSums(tt)
      rb = rb + colSums(pt)
      n_tot = n_tot + nrow(tt)
    }
    se = sqrt(Reduce(`+`, lapply(d_list, function(d) batch_mcse(d, 50L)^2))) / length(d_list)
    z = ifelse(se > 0, (emp / n_tot - rb / n_tot) / se, NA_real_)
    list(t_sampled = emp / n_tot, t_rb = rb / n_tot, z = z, n_draws = n_tot,
         n_unique_graphs = length(ls(cache)),
         max_abs_z = max(abs(z[rb / n_tot > 1e-3]), na.rm = TRUE))
  }
  smoke = function(fit, is_compare) {
    calls = list(
      summary = function() capture.output(print(summary(fit))),
      print = function() capture.output(print(fit)),
      coef = function() coef(fit),
      extract_sbm = function() extract_sbm(fit),
      extract_arguments = function() extract_arguments(fit),
      extract_indicators = function() extract_indicators(fit),
      extract_posterior_inclusion_probabilities = function() {
        extract_posterior_inclusion_probabilities(fit)
      }
    )
    if(is_compare) {
      calls$extract_group_params = function() extract_group_params(fit)
    } else {
      calls$extract_pairwise_interactions = function() extract_pairwise_interactions(fit)
      calls$extract_inclusion_bf = function() extract_inclusion_bf(fit)
    }
    vapply(calls, function(f) tryCatch({ f(); "ok" }, error = conditionMessage),
           character(1))
  }
  rows = list()
  for(r in names(gen)) {
    d = gen[[r]]()
    colnames(d$x) = vn
    seed = 20261101L + match(r, names(gen))
    t0 = Sys.time()
    if(r == "compare") {
      fit = bgmCompare(d$x, group_indicator = d$group, variable_type = d$type,
                       difference_prior = sbm_prior(), iter = 20000, warmup = 2000,
                       chains = 2, cores = 2, seed = seed,
                       display_progress = "none", verbose = FALSE)
      nm = fit$raw_samples$parameter_names$indicator
      keep = grepl(" (pairwise)", nm, fixed = TRUE)
      nm = sub(" (pairwise)", "", nm[keep], fixed = TRUE)
      ind_list = lapply(fit$raw_samples$indicator, function(m) m[, keep, drop = FALSE])
    } else {
      fit = bgm(d$x, variable_type = d$type, edge_prior = sbm_prior(),
                iter = 20000, warmup = 2000, chains = 2, cores = 2, seed = seed,
                display_progress = "none", verbose = FALSE)
      nm = fit$raw_samples$parameter_names$indicator
      ind_list = fit$raw_samples$indicator
    }
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ends = do.call(rbind, strsplit(nm, "-", fixed = TRUE))
    pidx = cbind(match(ends[, 1], vn), match(ends[, 2], vn))
    stopifnot(!anyNA(pidx))
    rb = rb_check(ind_list, fit$raw_samples$allocations, pidx)
    a = extract_arguments(fit)
    surface = if(label == "fixed") smoke(fit, r == "compare") else NULL
    cat(sprintf("route %-8s fit %4.0fs graphs=%5d | sampled T: %s | RB-exact: %s | max|z|=%.1f%s\n",
                r, secs, rb$n_unique_graphs,
                paste(sprintf("%.3f", rb$t_sampled[1:4]), collapse = " "),
                paste(sprintf("%.3f", rb$t_rb[1:4]), collapse = " "), rb$max_abs_z,
                if(is.null(surface)) "" else sprintf(" | surface ok %d/%d",
                                                     sum(surface == "ok"), length(surface))))
    rows[[r]] = list(route = r, seed = seed, fit_seconds = secs,
                     precision_graph_prior = if(is.null(a$precision_graph_prior)) NA else a$precision_graph_prior,
                     rb = rb, surface = surface)
  }
  rows
}


# ==============================================================================
# Run
# ==============================================================================

res = list(meta = list(
  label = label, lib = lib, tree = tree,
  build_commit = readLines(file.path(lib, "bgms", "BUILD_COMMIT"), n = 1L),
  version = format(packageVersion("bgms", lib.loc = lib)),
  date = format(Sys.time(), tz = "UTC", usetz = TRUE),
  session = sessionInfo()
))
if("grid" %in% parts) res$grid = arm_grid()
if("prior" %in% parts) res$prior = arm_prior()
if("routes" %in% parts) res$routes = arm_routes()
saveRDS(res, out)
cat("wrote", out, "\n")
