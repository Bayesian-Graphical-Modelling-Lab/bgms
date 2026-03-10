# ===========================================================================
# Shared helpers for the mixed MRF validation test suite
# ===========================================================================
# Provides:
#   - make_network()       : generate reproducible true parameter sets
#   - generate_data()      : simulate data from true parameters via bgms or
#                            mixedGM
#   - extract_bgms_blocks(): pull (mux, muy, Kxx, Kxy, Kyy) from bgms fit
#   - extract_mgm_blocks() : pull (mux, muy, Kxx, Kxy, Kyy) from mixedGM fit
#   - flatten_params()     : flatten all blocks to a single named vector
#   - recovery_table()     : compare estimated vs true as a data.frame
#   - recovery_scatter()   : scatterplot of estimated vs true
#   - summarise_recovery() : one-line correlation + bias + RMSE summary
# ===========================================================================

# ------------------------------------------------------------------
# make_network
# ------------------------------------------------------------------
# Build a reproducible mixed MRF parameter set.
#
# @param p              Number of discrete variables.
# @param q              Number of continuous variables.
# @param n_cat          Integer vector of length p: number of categories per
#                       discrete variable (bgms convention = max index,
#                       e.g. binary = 1).
# @param variable_type  Character vector of length p: "ordinal" or
#                       "blume-capel" per discrete variable. Default: all
#                       ordinal.
# @param baseline_category  Integer vector of length p: baseline category
#                       for Blume-Capel variables. Ignored for ordinal.
#                       Default: all zeros.
# @param density        Approximate fraction of non-zero edges.
# @param seed           Random seed for reproducibility.
#
# Returns: named list with mux, muy, Kxx, Kxy, Kyy, n_cat, p, q,
#   variable_type, baseline_category.
# ------------------------------------------------------------------
make_network = function(p, q, n_cat, variable_type = rep("ordinal", p),
                        baseline_category = rep(0L, p),
                        density = 0.5, seed = 42) {
  set.seed(seed)

  max_cat = max(n_cat)
  # BC variables use 2 columns (alpha, beta); ensure mux is wide enough
  mux_cols = max(max_cat, 2L)

  # --- Main effects (mux) ---
  # Ordinal: threshold parameters, NA-padded.
  # Blume-Capel: alpha (linear) and beta (quadratic), rest NA.
  mux = matrix(NA_real_, nrow = p, ncol = mux_cols)
  for(i in seq_len(p)) {
    if(variable_type[i] == "blume-capel") {
      mux[i, 1] = round(runif(1, -0.3, 0.3), 2) # alpha
      mux[i, 2] = round(runif(1, -0.4, -0.05), 2) # beta (negative = peaked)
    } else {
      vals = sort(round(seq(-0.5, 0.5, length.out = n_cat[i]), 2))
      if(n_cat[i] == 1) vals = round(runif(1, -0.3, 0.3), 2)
      mux[i, seq_len(n_cat[i])] = vals
    }
  }

  # --- Continuous means ---
  muy = round(runif(q, -0.5, 0.5), 2)

  # --- Kxx: ordinal-ordinal interactions (symmetric, zero diagonal) ---
  n_edges_xx = p * (p - 1) / 2
  mask_xx = rbinom(n_edges_xx, 1, density)
  vals_xx = mask_xx * round(runif(n_edges_xx, 0.15, 0.4) * sample(c(-1, 1), n_edges_xx, replace = TRUE), 2)
  Kxx = matrix(0, p, p)
  Kxx[upper.tri(Kxx)] = vals_xx
  Kxx = Kxx + t(Kxx)

  # --- Kyy: continuous precision (positive definite, sparse off-diag with diagonal dominance) ---
  n_edges_yy = q * (q - 1) / 2
  mask_yy = rbinom(n_edges_yy, 1, density)
  vals_yy = mask_yy * round(runif(n_edges_yy, 0.05, 0.2) * sample(c(-1, 1), n_edges_yy, replace = TRUE), 2)
  Kyy = matrix(0, q, q)
  Kyy[upper.tri(Kyy)] = vals_yy
  Kyy = Kyy + t(Kyy)
  diag(Kyy) = abs(rowSums(Kyy)) + runif(q, 1.2, 1.8)

  # --- Kxy: ordinal-continuous cross interactions ---
  n_cross = p * q
  mask_xy = rbinom(n_cross, 1, density)
  vals_xy = mask_xy * round(runif(n_cross, 0.1, 0.3) * sample(c(-1, 1), n_cross, replace = TRUE), 2)
  Kxy = matrix(vals_xy, p, q)

  list(
    mux = mux, muy = muy,
    Kxx = Kxx, Kxy = Kxy, Kyy = Kyy,
    n_cat = n_cat, p = p, q = q,
    variable_type = variable_type,
    baseline_category = as.integer(baseline_category)
  )
}

# ------------------------------------------------------------------
# generate_data
# ------------------------------------------------------------------
# Simulate data from true network parameters.
#
# @param net        Output from make_network().
# @param n          Number of observations.
# @param source     "bgms" or "mixedGM".
# @param iter       Gibbs burn-in iterations.
# @param seed       Random seed.
#
# Returns: data.frame (n x (p+q)) with ordinal columns first, then continuous.
# ------------------------------------------------------------------
generate_data = function(net, n, source = "bgms", iter = 1000L, seed = 1) {
  if(source == "bgms") {
    # simulate_mrf() does not support mixed types; use the C++ Gibbs
    # sampler directly via sample_mixed_mrf_gibbs().
    sim = bgms:::sample_mixed_mrf_gibbs(
      num_states = as.integer(n),
      Kxx_r = net$Kxx,
      Kxy_r = net$Kxy,
      Kyy_r = net$Kyy,
      mux_r = net$mux,
      muy_r = net$muy,
      num_categories_r = as.integer(net$n_cat),
      variable_type_r = net$variable_type,
      baseline_category_r = net$baseline_category,
      iter = as.integer(iter),
      seed = as.integer(seed)
    )
    df = as.data.frame(cbind(sim$x, sim$y))
    names(df) = c(paste0("X", seq_len(net$p)), paste0("Y", seq_len(net$q)))
    df
  } else if(source == "mixedGM") {
    sim = mixedGM::mixed_gibbs_generate(
      n = n,
      Kxx = net$Kxx, Kxy = net$Kxy, Kyy = net$Kyy,
      mux = net$mux, muy = net$muy,
      num_categories = net$n_cat + 1L,
      n_burnin = iter
    )
    df = as.data.frame(cbind(sim$x, sim$y))
    names(df) = c(paste0("X", seq_len(net$p)), paste0("Y", seq_len(net$q)))
    df
  } else {
    stop('generate_data(): source must be "bgms" or "mixedGM", not "', source, '".')
  }
}

# ------------------------------------------------------------------
# extract_bgms_blocks
# ------------------------------------------------------------------
# Extract parameter blocks from a bgms fit object.
#
# @param fit  A bgms object.
# @param net  The true network (for dimension reference).
#
# Returns: list(mux, muy, Kxx, Kxy, Kyy).
# ------------------------------------------------------------------
extract_bgms_blocks = function(fit, net) {
  pm = coef(fit)
  mux = pm$main$discrete # p x max_cat
  muy_vec = pm$main$continuous[, "mean"] # length q
  Kyy_diag = pm$main$continuous[, "precision"] # length q

  pw = pm$pairwise # (p+q) x (p+q)
  p = net$p
  q = net$q
  Kxx = pw[seq_len(p), seq_len(p)]
  Kxy = pw[seq_len(p), p + seq_len(q)]
  Kyy_off = pw[p + seq_len(q), p + seq_len(q)]
  Kyy = Kyy_off
  diag(Kyy) = Kyy_diag

  list(mux = mux, muy = muy_vec, Kxx = Kxx, Kxy = Kxy, Kyy = Kyy)
}

# ------------------------------------------------------------------
# extract_mgm_blocks
# ------------------------------------------------------------------
# Extract posterior means from a mixedGM fit object.
#
# @param fit    Output of mixedGM::mixed_sampler().
# @param n_cat  Integer vector of length p: number of categories per ordinal
#               variable (bgms convention). Used to NA-pad unused threshold
#               slots so flatten_params() produces matching lengths.
#
# Returns: list(mux, muy, Kxx, Kxy, Kyy).
# ------------------------------------------------------------------
extract_mgm_blocks = function(fit, n_cat) {
  mux = apply(fit$samples$mux, c(2, 3), mean)
  for(i in seq_along(n_cat)) {
    if(n_cat[i] < ncol(mux)) {
      mux[i, (n_cat[i] + 1):ncol(mux)] = NA_real_
    }
  }
  list(
    mux = mux,
    muy = colMeans(fit$samples$muy),
    Kxx = apply(fit$samples$Kxx, c(2, 3), mean),
    Kxy = apply(fit$samples$Kxy, c(2, 3), mean),
    Kyy = apply(fit$samples$Kyy, c(2, 3), mean)
  )
}

# ------------------------------------------------------------------
# flatten_params
# ------------------------------------------------------------------
# Flatten parameter blocks to a single named vector for comparison.
# Excludes NA entries in mux.
#
# @param blocks  list(mux, muy, Kxx, Kxy, Kyy).
# @param prefix  Optional prefix for names.
#
# Returns: named numeric vector.
# ------------------------------------------------------------------
flatten_params = function(blocks, prefix = "") {
  mux_vals = as.vector(t(blocks$mux))
  mux_keep = !is.na(mux_vals)
  mux_named = mux_vals[mux_keep]
  names(mux_named) = paste0(prefix, "mux_", which(mux_keep))

  muy_named = blocks$muy
  names(muy_named) = paste0(prefix, "muy_", seq_along(muy_named))

  # Kxx upper triangle
  kxx_ut = blocks$Kxx[upper.tri(blocks$Kxx)]
  names(kxx_ut) = paste0(prefix, "Kxx_", seq_along(kxx_ut))

  # Kxy full
  kxy_vals = as.vector(blocks$Kxy)
  names(kxy_vals) = paste0(prefix, "Kxy_", seq_along(kxy_vals))

  # Kyy upper triangle (includes diagonal)
  kyy_ut = blocks$Kyy[upper.tri(blocks$Kyy, diag = TRUE)]
  names(kyy_ut) = paste0(prefix, "Kyy_", seq_along(kyy_ut))

  c(mux_named, muy_named, kxx_ut, kxy_vals, kyy_ut)
}

# ------------------------------------------------------------------
# recovery_table
# ------------------------------------------------------------------
# Build a data.frame comparing true vs estimated parameters.
#
# @param true_blocks  list(mux, muy, Kxx, Kxy, Kyy).
# @param est_blocks   list(mux, muy, Kxx, Kxy, Kyy).
# @param label        Character label for the method.
#
# Returns: data.frame with columns: parameter, true, estimate, diff, block, method.
# ------------------------------------------------------------------
recovery_table = function(true_blocks, est_blocks, label = "estimate") {
  true_flat = flatten_params(true_blocks)
  est_flat = flatten_params(est_blocks)

  # Determine block membership
  block = ifelse(grepl("^mux", names(true_flat)), "mux",
    ifelse(grepl("^muy", names(true_flat)), "muy",
      ifelse(grepl("^Kxx", names(true_flat)), "Kxx",
        ifelse(grepl("^Kxy", names(true_flat)), "Kxy", "Kyy")
      )
    )
  )

  data.frame(
    parameter = names(true_flat),
    true = true_flat,
    estimate = est_flat,
    diff = est_flat - true_flat,
    block = block,
    method = label,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ------------------------------------------------------------------
# recovery_scatter
# ------------------------------------------------------------------
# Scatterplot of estimated vs true, coloured by block.
#
# @param tab    Output of recovery_table() (possibly rbind of multiple).
# @param main   Plot title.
# ------------------------------------------------------------------
recovery_scatter = function(tab, main = "Parameter recovery") {
  block_cols = c(
    mux = "#E41A1C", muy = "#377EB8", Kxx = "#4DAF4A",
    Kxy = "#984EA3", Kyy = "#FF7F00"
  )
  rng = range(c(tab$true, tab$estimate), na.rm = TRUE)
  rng = rng + diff(rng) * c(-0.05, 0.05)

  plot(tab$true, tab$estimate,
    pch = 19, cex = 0.9,
    col = adjustcolor(block_cols[tab$block], 0.7),
    xlim = rng, ylim = rng, asp = 1,
    xlab = "True value", ylab = "Posterior mean",
    main = main
  )
  abline(0, 1, lty = 2, col = "grey40")

  r = cor(tab$true, tab$estimate)
  rmse = sqrt(mean(tab$diff^2))
  legend("topleft",
    legend = c(
      names(block_cols),
      sprintf("r = %.3f", r),
      sprintf("RMSE = %.3f", rmse)
    ),
    col = c(block_cols, NA, NA),
    pch = c(rep(19, 5), NA, NA),
    bty = "n", cex = 0.8
  )
}

# ------------------------------------------------------------------
# summarise_recovery
# ------------------------------------------------------------------
# Print a one-line summary of recovery quality.
#
# @param tab    Output of recovery_table().
# @param label  String label.
#
# Returns: invisible data.frame with summary stats.
# ------------------------------------------------------------------
summarise_recovery = function(tab, label = "") {
  r = cor(tab$true, tab$estimate)
  bias = mean(tab$diff)
  rmse = sqrt(mean(tab$diff^2))
  max_diff = max(abs(tab$diff))
  cat(sprintf(
    "[%s] r = %.4f | bias = %.4f | RMSE = %.4f | max|diff| = %.4f\n",
    label, r, bias, rmse, max_diff
  ))

  # Per-block summary
  blocks = unique(tab$block)
  block_summary = do.call(rbind, lapply(blocks, function(b) {
    sub = tab[tab$block == b, ]
    r_val = if(nrow(sub) >= 3 && sd(sub$true) > 0 && sd(sub$estimate) > 0) {
      cor(sub$true, sub$estimate)
    } else {
      NA_real_
    }
    data.frame(
      block = b,
      n = nrow(sub),
      r = r_val,
      bias = mean(sub$diff),
      rmse = sqrt(mean(sub$diff^2)),
      max_diff = max(abs(sub$diff)),
      stringsAsFactors = FALSE
    )
  }))

  cat("  Per-block breakdown:\n")
  print(block_summary, row.names = FALSE)
  invisible(block_summary)
}

# ------------------------------------------------------------------
# agreement_scatter
# ------------------------------------------------------------------
# Scatterplot of method A vs method B estimates, coloured by block.
#
# @param tab_a  recovery_table for method A.
# @param tab_b  recovery_table for method B.
# @param label_a, label_b  Axis labels.
# @param main   Plot title.
# ------------------------------------------------------------------
agreement_scatter = function(tab_a, tab_b, label_a = "Method A",
                             label_b = "Method B", main = "Agreement") {
  block_cols = c(
    mux = "#E41A1C", muy = "#377EB8", Kxx = "#4DAF4A",
    Kxy = "#984EA3", Kyy = "#FF7F00"
  )
  rng = range(c(tab_a$estimate, tab_b$estimate), na.rm = TRUE)
  rng = rng + diff(rng) * c(-0.05, 0.05)

  plot(tab_a$estimate, tab_b$estimate,
    pch = 19, cex = 0.9,
    col = adjustcolor(block_cols[tab_a$block], 0.7),
    xlim = rng, ylim = rng, asp = 1,
    xlab = label_a, ylab = label_b, main = main
  )
  abline(0, 1, lty = 2, col = "grey40")

  r = cor(tab_a$estimate, tab_b$estimate)
  rmse = sqrt(mean((tab_a$estimate - tab_b$estimate)^2))
  legend("topleft",
    legend = c(
      names(block_cols),
      sprintf("r = %.3f", r),
      sprintf("RMSE = %.3f", rmse)
    ),
    col = c(block_cols, NA, NA),
    pch = c(rep(19, 5), NA, NA),
    bty = "n", cex = 0.8
  )
}

# ------------------------------------------------------------------
# trace_panel
# ------------------------------------------------------------------
# Plot trace + density panels for selected parameters from raw samples.
#
# @param samples  Matrix of MCMC samples (iterations x parameters).
# @param names    Column names to plot.
# @param true_vals  Optional named vector of true values.
# @param main     Overall title.
# ------------------------------------------------------------------
trace_panel = function(samples, names = NULL, true_vals = NULL, main = "") {
  if(is.null(names)) names = colnames(samples)[seq_len(min(6, ncol(samples)))]
  n_par = length(names)
  par(mfrow = c(n_par, 2), mar = c(3, 3, 2, 1), mgp = c(2, 0.6, 0))
  for(nm in names) {
    vals = samples[, nm]
    # Trace
    plot(vals,
      type = "l", col = adjustcolor("grey30", 0.5),
      main = paste(nm, "trace"), xlab = "Iteration", ylab = nm, cex.main = 0.9
    )
    if(!is.null(true_vals) && nm %in% names(true_vals)) {
      abline(h = true_vals[nm], col = "red", lwd = 2)
    }
    abline(h = mean(vals), col = "steelblue", lwd = 1.5)
    # Density
    d = density(vals)
    plot(d, main = paste(nm, "density"), xlab = nm, cex.main = 0.9)
    if(!is.null(true_vals) && nm %in% names(true_vals)) {
      abline(v = true_vals[nm], col = "red", lwd = 2)
    }
    abline(v = mean(vals), col = "steelblue", lwd = 1.5)
  }
}
