# --------------------------------------------------------------------------- #
# Tests for the mixed MRF RATTLE projection methods.
#
# Validates:
#  1. logp_and_gradient_full against central finite differences
#  2. SHAKE position projection (idempotency, constraint enforcement)
#  3. RATTLE momentum projection (idempotency, orthogonality)
#  4. Consistency between active-space and full-space gradient
# --------------------------------------------------------------------------- #


# ---- Helpers ----------------------------------------------------------------

# Compute the full-space parameter dimension for a given setup
mixed_full_dim = function(num_cats, is_ord, p, q) {
  n_main = sum(ifelse(is_ord == 1L, num_cats, 2L))
  n_pairwise_xx = p * (p - 1) / 2
  n_cross = p * q
  n_chol = q * (q + 1) / 2
  n_main + n_pairwise_xx + q + n_cross + n_chol
}

# Finite-difference gradient for logp_and_gradient_full
mixed_fd_gradient_full = function(params, x, y, num_cats, is_ord, base_cat,
                                  edge_ind, pl_mode, scale, eps = 1e-5) {
  n_total = length(params)
  fd = numeric(n_total)
  for(k in seq_len(n_total)) {
    p_plus = params
    p_minus = params
    p_plus[k] = p_plus[k] + eps
    p_minus[k] = p_minus[k] - eps
    fp = mixed_test_logp_and_gradient_full(
      p_plus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, pl_mode, scale
    )$value
    fm = mixed_test_logp_and_gradient_full(
      p_minus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, pl_mode, scale
    )$value
    fd[k] = (fp - fm) / (2 * eps)
  }
  fd
}

mixed_check_gradient_full = function(params, x, y, num_cats, is_ord, base_cat,
                                     edge_ind, pl_mode, scale, eps = 1e-5) {
  res = mixed_test_logp_and_gradient_full(
    params, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, pl_mode, scale
  )
  fd = mixed_fd_gradient_full(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, pl_mode, scale, eps
  )
  ag = res$gradient
  denom = pmax(abs(ag), abs(fd), 1)
  rel_err = abs(ag - fd) / denom
  max(rel_err)
}

# Compute excluded indices in the full parameter vector for a given edge setup
# Returns list with kxx_indices, kxy_indices (1-based R indices)
mixed_excluded_indices = function(num_cats, is_ord, p, q, edge_ind) {
  n_main = sum(ifelse(is_ord == 1L, num_cats, 2L))
  n_pairwise_xx = p * (p - 1) / 2

  # Kxx excluded indices (1-based)
  kxx_offset = n_main
  kxx_excluded = integer(0)
  idx = 0L
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      if(edge_ind[i, j] == 0L) {
        kxx_excluded = c(kxx_excluded, kxx_offset + idx + 1L)
      }
      idx = idx + 1L
    }
  }

  # Kxy excluded indices (1-based)
  kxy_offset = n_main + n_pairwise_xx + q
  kxy_excluded = integer(0)
  idx = 0L
  for(i in seq_len(p)) {
    for(j in seq_len(q)) {
      if(edge_ind[i, p + j] == 0L) {
        kxy_excluded = c(kxy_excluded, kxy_offset + idx + 1L)
      }
      idx = idx + 1L
    }
  }

  list(kxx = kxx_excluded, kxy = kxy_excluded)
}


# ---- Full-space gradient FD tests ------------------------------------------ #

test_that("full gradient matches FD, dense edges (p=3, q=2, conditional)", {
  set.seed(42)
  n = 80
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(123)
  params = rnorm(dim, sd = 0.3)
  err = mixed_check_gradient_full(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, "conditional", 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("full gradient matches FD, dense edges (p=3, q=2, marginal)", {
  set.seed(42)
  n = 80
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(123)
  params = rnorm(dim, sd = 0.3)
  err = mixed_check_gradient_full(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, "marginal", 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("full gradient matches FD, sparse edges (p=3, q=2, conditional)", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  # Sparse: only edge (1,2) discrete, edges (1,4) and (2,5) cross, no Gyy edges
  edge_ind = matrix(0L, total, total)
  edge_ind[1, 2] = edge_ind[2, 1] = 1L
  edge_ind[1, 4] = edge_ind[4, 1] = 1L
  edge_ind[2, 5] = edge_ind[5, 2] = 1L
  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(4)
  params = rnorm(dim, sd = 0.3)
  err = mixed_check_gradient_full(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, "conditional", 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("full gradient matches FD, sparse with Gyy edges (p=3, q=3, conditional)", {
  set.seed(77)
  n = 80
  p = 3
  q = 3
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  # Sparse Gyy: only edge (4,5) = continuous (1,2)
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[5, 6] = edge_ind[6, 5] = 0L # remove continuous edge (2,3)
  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(5)
  params = rnorm(dim, sd = 0.2)
  err = mixed_check_gradient_full(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, "conditional", 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("full gradient matches FD, mixed ord+BC (p=3, q=2, conditional)", {
  set.seed(99)
  n = 70
  p = 3
  q = 2
  x = matrix(sample(0:3, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = c(3L, 3L, 3L)
  is_ord = c(1L, 0L, 1L)
  base_cat = c(0L, 1L, 0L)
  total = p + q
  # Sparse: remove some edges
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L # remove discrete edge (1,3)
  edge_ind[2, 4] = edge_ind[4, 2] = 0L # remove cross edge (2,1)
  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(2)
  params = rnorm(dim, sd = 0.3)
  err = mixed_check_gradient_full(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, "conditional", 2.5
  )
  expect_lt(err, 1e-5)
})


# ---- Full gradient value consistency ---------------------------------------- #

test_that("full gradient agrees with active gradient when all edges present", {
  set.seed(42)
  n = 80
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(123)
  params = rnorm(dim, sd = 0.3)

  res_full = mixed_test_logp_and_gradient_full(
    params, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  res_active = mixed_test_logp_and_gradient(
    params, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  # With all edges active, full and active dimensions are the same
  expect_equal(res_full$value, res_active$value, tolerance = 1e-10)
  expect_equal(res_full$gradient, res_active$gradient, tolerance = 1e-10)
})


# ---- SHAKE position projection tests --------------------------------------- #

test_that("SHAKE projection zeros excluded Kxx and Kxy entries", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  # Remove edges (1,2) discrete, (1,4) cross
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 2] = edge_ind[2, 1] = 0L
  edge_ind[1, 4] = edge_ind[4, 1] = 0L

  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(10)
  params = rnorm(dim, sd = 0.5)
  inv_mass = rep(1.0, dim)

  excl = mixed_excluded_indices(num_cats, is_ord, p, q, edge_ind)
  # Verify some entries are non-zero before projection
  expect_true(any(params[excl$kxx] != 0) || any(params[excl$kxy] != 0))

  res = mixed_test_project_position(
    params, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  proj = res$projected

  # Excluded entries must be zero
  if(length(excl$kxx) > 0) expect_true(all(proj[excl$kxx] == 0))
  if(length(excl$kxy) > 0) expect_true(all(proj[excl$kxy] == 0))
})

test_that("SHAKE projection is idempotent", {
  set.seed(42)
  n = 60
  p = 3
  q = 3
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  # Sparse: remove Gyy edge (4,6) and Gxx edge (1,3)
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[4, 6] = edge_ind[6, 4] = 0L

  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(10)
  params = rnorm(dim, sd = 0.5)
  inv_mass = rep(1.0, dim)

  # First projection
  res1 = mixed_test_project_position(
    params, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  # Second projection of the result
  res2 = mixed_test_project_position(
    res1$projected, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  expect_equal(res1$projected, res2$projected, tolerance = 1e-12)
})

test_that("SHAKE projection is identity when all edges present", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L

  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(10)
  params = rnorm(dim, sd = 0.5)
  inv_mass = rep(1.0, dim)

  res = mixed_test_project_position(
    params, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  expect_equal(as.vector(res$projected), params, tolerance = 1e-12)
})


# ---- RATTLE momentum projection tests -------------------------------------- #

test_that("RATTLE projection zeros excluded Kxx and Kxy momentum entries", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 2] = edge_ind[2, 1] = 0L
  edge_ind[2, 5] = edge_ind[5, 2] = 0L

  dim = mixed_full_dim(num_cats, is_ord, p, q)
  set.seed(10)
  # Position must be on the constraint manifold first
  pos = rnorm(dim, sd = 0.3)
  inv_mass = rep(1.0, dim)
  pos_proj = mixed_test_project_position(
    pos, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )$projected

  # Random momentum
  set.seed(20)
  mom = rnorm(dim, sd = 1.0)

  excl = mixed_excluded_indices(num_cats, is_ord, p, q, edge_ind)

  res = mixed_test_project_momentum(
    mom, pos_proj, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  proj_mom = res$projected

  # Excluded Kxx and Kxy momentum must be zero
  if(length(excl$kxx) > 0) expect_true(all(abs(proj_mom[excl$kxx]) < 1e-10))
  if(length(excl$kxy) > 0) expect_true(all(abs(proj_mom[excl$kxy]) < 1e-10))
})

test_that("RATTLE momentum projection is idempotent", {
  set.seed(42)
  n = 60
  p = 3
  q = 3
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[4, 6] = edge_ind[6, 4] = 0L

  dim = mixed_full_dim(num_cats, is_ord, p, q)
  inv_mass = rep(1.0, dim)

  set.seed(10)
  pos = rnorm(dim, sd = 0.3)
  pos_proj = mixed_test_project_position(
    pos, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )$projected

  set.seed(20)
  mom = rnorm(dim, sd = 1.0)

  res1 = mixed_test_project_momentum(
    mom, pos_proj, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  res2 = mixed_test_project_momentum(
    res1$projected, pos_proj, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  expect_equal(res1$projected, res2$projected, tolerance = 1e-10)
})

test_that("RATTLE projection is identity when all edges present", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L

  dim = mixed_full_dim(num_cats, is_ord, p, q)
  inv_mass = rep(1.0, dim)

  set.seed(10)
  pos = rnorm(dim, sd = 0.3)
  set.seed(20)
  mom = rnorm(dim, sd = 1.0)

  res = mixed_test_project_momentum(
    mom, pos, inv_mass, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", 2.5
  )
  expect_equal(as.vector(res$projected), mom, tolerance = 1e-12)
})


# ---- Integration smoke test ------------------------------------------------ #

test_that("bgm with NUTS + edge_selection runs on mixed data", {
  skip_on_cran()
  set.seed(42)
  n = 80
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  data = data.frame(x, y)
  colnames(data) = c(paste0("d", seq_len(p)), paste0("c", seq_len(q)))
  variable_type = c(rep("ordinal", p), rep("continuous", q))

  fit = bgm(
    data,
    variable_type = variable_type,
    iter = 100,
    warmup = 50,
    edge_selection = TRUE,
    update_method = "nuts"
  )

  expect_s3_class(fit, "bgms")
  expect_true(!is.null(fit$posterior_mean_associations))
})
