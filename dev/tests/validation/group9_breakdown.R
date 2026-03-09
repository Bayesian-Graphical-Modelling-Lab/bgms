# ===========================================================================
# Targeted timing breakdown for mixed MRF hotspot identification
# ===========================================================================
library(bgms)
source("dev/tests/validation/helpers.R")

cat("===================================================================\n")
cat("  Mixed MRF Timing Breakdown\n")
cat("===================================================================\n\n")

# ---------------------------------------------------------------------------
# Setup: p=8, q=6 problem
# ---------------------------------------------------------------------------
p = 8; q = 6; n = 300
n_cat = rep(2L, p)
net = make_network(p = p, q = q, n_cat = n_cat, density = 0.4, seed = 42)
dat = generate_data(net, n = n, seed = 42)
vt = c(rep("ordinal", p), rep("continuous", q))

n_iter = 200; n_warmup = 100

cat(sprintf("Config: p=%d, q=%d, n=%d, iter=%d, warmup=%d\n\n", p, q, n, n_iter, n_warmup))

# ---------------------------------------------------------------------------
# A. Component breakdown: estimation-only vs edge-selection
# ---------------------------------------------------------------------------
cat("--- A. Edge selection overhead ---\n")

t_est = system.time(suppressWarnings(
  bgm(dat, variable_type = vt, iter = n_iter, warmup = n_warmup,
      chains = 1, cores = 1, edge_selection = FALSE,
      display_progress = FALSE)))["elapsed"]

t_es = system.time(suppressWarnings(
  bgm(dat, variable_type = vt, iter = n_iter, warmup = n_warmup,
      chains = 1, cores = 1, edge_selection = TRUE,
      display_progress = FALSE)))["elapsed"]

cat(sprintf("  Estimation only:  %6.2f s (%5.2f ms/iter)\n", t_est, t_est/n_iter*1000))
cat(sprintf("  + Edge selection: %6.2f s (%5.2f ms/iter)\n", t_es, t_es/n_iter*1000))
cat(sprintf("  Edge selection adds: %.2f ms/iter (%.0f%%)\n\n",
            (t_es - t_est)/n_iter*1000, (t_es - t_est)/t_es*100))

# ---------------------------------------------------------------------------
# B. Conditional vs Marginal PL at varying p
# ---------------------------------------------------------------------------
cat("--- B. Conditional vs Marginal PL scaling ---\n")

for(pp in c(4, 8, 12)) {
  qq = round(pp * 0.75)
  nc = rep(2L, pp)
  net2 = make_network(p = pp, q = qq, n_cat = nc, density = 0.4, seed = 42)
  dat2 = generate_data(net2, n = n, seed = 42)
  vt2 = c(rep("ordinal", pp), rep("continuous", qq))

  tc = system.time(suppressWarnings(
    bgm(dat2, variable_type = vt2, iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        display_progress = FALSE)))["elapsed"]

  tm = system.time(suppressWarnings(
    bgm(dat2, variable_type = vt2, iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        pseudolikelihood = "marginal",
        display_progress = FALSE)))["elapsed"]

  cat(sprintf("  p=%2d q=%2d : cond=%5.2f, marg=%5.2f ms/iter (%.1fx)\n",
              pp, qq, tc/n_iter*1000, tm/n_iter*1000, tm/tc))
}

# ---------------------------------------------------------------------------
# C. Pure sub-model comparison (isolate cross-term cost)
# ---------------------------------------------------------------------------
cat("\n--- C. Cross-model comparison ---\n")

for(sz in c(4, 8, 12)) {
  nc = rep(2L, sz)
  net3 = make_network(p = sz, q = sz, n_cat = nc, density = 0.4, seed = 42)
  dat3 = generate_data(net3, n = n, seed = 42)
  vt3 = c(rep("ordinal", sz), rep("continuous", sz))

  # Mixed
  t_mix = system.time(suppressWarnings(
    bgm(dat3, variable_type = vt3, iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        display_progress = FALSE)))["elapsed"]

  # OMRF only (discrete columns)
  dat_o = dat3[, 1:sz, drop = FALSE]
  t_omrf = system.time(suppressWarnings(
    bgm(dat_o, variable_type = rep("ordinal", sz), iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        display_progress = FALSE)))["elapsed"]

  # GGM only (continuous columns)
  dat_g = dat3[, (sz+1):(2*sz), drop = FALSE]
  t_ggm = system.time(suppressWarnings(
    bgm(dat_g, variable_type = rep("continuous", sz), iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        display_progress = FALSE)))["elapsed"]

  omrf_ms = t_omrf/n_iter*1000
  ggm_ms = t_ggm/n_iter*1000
  mix_ms = t_mix/n_iter*1000
  cross_ms = mix_ms - omrf_ms - ggm_ms

  cat(sprintf("  size=%2d : mixed=%5.2f  omrf=%5.2f  ggm=%5.2f  cross(residual)=%5.2f ms/iter\n",
              sz, mix_ms, omrf_ms, ggm_ms, cross_ms))
}

# ---------------------------------------------------------------------------
# D. Scaling: isolate p-scaling (fix q=4) and q-scaling (fix p=4)
# ---------------------------------------------------------------------------
cat("\n--- D. Scaling curves ---\n")

cat("  Fix q=4, vary p:\n")
for(pp in c(3, 6, 9, 12)) {
  nc = rep(2L, pp)
  net4 = make_network(p = pp, q = 4, n_cat = nc, density = 0.4, seed = 42)
  dat4 = generate_data(net4, n = n, seed = 42)
  vt4 = c(rep("ordinal", pp), rep("continuous", 4))

  tt = system.time(suppressWarnings(
    bgm(dat4, variable_type = vt4, iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        display_progress = FALSE)))["elapsed"]

  cat(sprintf("    p=%2d q=4 : %5.2f ms/iter\n", pp, tt/n_iter*1000))
}

cat("  Fix p=4, vary q:\n")
for(qq in c(2, 4, 6, 8, 10)) {
  nc = rep(2L, 4)
  net5 = make_network(p = 4, q = qq, n_cat = nc, density = 0.4, seed = 42)
  dat5 = generate_data(net5, n = n, seed = 42)
  vt5 = c(rep("ordinal", 4), rep("continuous", qq))

  tt = system.time(suppressWarnings(
    bgm(dat5, variable_type = vt5, iter = n_iter, warmup = n_warmup,
        chains = 1, cores = 1, edge_selection = FALSE,
        display_progress = FALSE)))["elapsed"]

  cat(sprintf("    p=4 q=%2d : %5.2f ms/iter\n", qq, tt/n_iter*1000))
}

# ---------------------------------------------------------------------------
# E. Operation counts and expected complexity
# ---------------------------------------------------------------------------
cat("\n--- E. Complexity analysis ---\n")
cat("Per-iteration MH proposals for p=8, q=6:\n")

mux_ops = sum(rep(2L, p))  # binary ordinal: 2 categories → 1 threshold each
muy_ops = q
kxx_ops = p * (p - 1) / 2
kyy_ops = q * (q - 1) / 2 + q  # off-diag + diag
kxy_ops = p * q

cat(sprintf("  Main effects (mux):  %3d proposals\n", mux_ops))
cat(sprintf("  Continuous means:    %3d proposals\n", muy_ops))
cat(sprintf("  Kxx updates:         %3d proposals\n", kxx_ops))
cat(sprintf("  Kyy updates:         %3d proposals (off-diag: %d, diag: %d)\n",
            kyy_ops, q*(q-1)/2, q))
cat(sprintf("  Kxy updates:         %3d proposals\n", kxy_ops))
cat(sprintf("  TOTAL:               %3d proposals/iter\n",
            mux_ops + muy_ops + kxx_ops + kyy_ops + kxy_ops))

# Edge selection adds:
es_ops = kxx_ops + q*(q-1)/2 + kxy_ops
cat(sprintf("  Edge selection:      %3d proposals/iter (Kxx:%d + Kyy:%d + Kxy:%d)\n",
            es_ops, kxx_ops, q*(q-1)/2, kxy_ops))

cat("\n===================================================================\n")
cat("  Done\n")
cat("===================================================================\n")
