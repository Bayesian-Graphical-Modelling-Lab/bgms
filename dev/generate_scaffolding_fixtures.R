# ==============================================================================
# Generate Golden-Snapshot Fixtures for R Scaffolding Refactor
# ==============================================================================
#
# Phase A-0 of the scaffolding refactor (dev/scaffolding/plan.md).
#
# This script captures the INTERMEDIATE outputs of the current validation and
# preprocessing functions — check_model(), check_compare_model(),
# reformat_data(), and compare_reformat_data() — for representative inputs.
# These fixtures are mechanical oracles: every refactored validator must
# reproduce exactly the same outputs.
#
# Unlike the refactor fixtures (which capture full model fits), these fixtures
# are FAST — no sampling is involved. They test the code we're about to
# restructure.
#
# Output: dev/fixtures/scaffolding/
#   - One .rds per fixture case
#   - A manifest.rds listing all cases
#
# Usage:
#   Rscript dev/generate_scaffolding_fixtures.R
#
# ==============================================================================

library(bgms)

# These are internal functions — access via :::
check_model          <- bgms:::check_model
check_compare_model  <- bgms:::check_compare_model
reformat_data        <- bgms:::reformat_data
compare_reformat_data <- bgms:::compare_reformat_data

fixture_dir <- file.path("dev", "fixtures", "scaffolding")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# ==============================================================================
# Helper: generate small synthetic datasets
# ==============================================================================

# Small ordinal dataset (5 variables, 3 categories each: 0, 1, 2)
make_ordinal_data <- function(n = 50, p = 5, max_cat = 2) {
  x <- matrix(sample(0:max_cat, n * p, replace = TRUE), nrow = n, ncol = p)
  colnames(x) <- paste0("V", seq_len(p))
  x
}

# Small continuous dataset
make_continuous_data <- function(n = 50, p = 5) {
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("V", seq_len(p))
  x
}

# Inject NAs into a dataset
inject_nas <- function(x, prop = 0.05) {
  n_na <- max(1, floor(nrow(x) * ncol(x) * prop))
  idx <- sample(length(x), n_na)
  x[idx] <- NA
  x
}

# Make a dataset where one category is missing from one group
make_missing_category_data <- function() {
  # Group 1: categories 0, 1, 2 for all variables
  x1 <- matrix(sample(0:2, 30 * 4, replace = TRUE), nrow = 30, ncol = 4)
  # Group 2: variable 1 only has categories 0, 1 (no 2)
  x2 <- matrix(sample(0:2, 30 * 4, replace = TRUE), nrow = 30, ncol = 4)
  x2[, 1] <- sample(0:1, 30, replace = TRUE)

  x <- rbind(x1, x2)
  colnames(x) <- paste0("V", 1:4)
  group <- c(rep(1L, 30), rep(2L, 30))
  list(x = x, group = group)
}

# ==============================================================================
# Fixture definition
# ==============================================================================

# Each fixture stores:
#   $id             - unique name
#   $desc           - human-readable description
#   $type           - "bgm" or "compare"
#   $input          - the exact arguments passed to the functions
#   $check_model    - return value of check_model() or check_compare_model()
#   $reformat_data  - return value of reformat_data() or compare_reformat_data()

fixtures <- list()

# ---------------------------------------------------------------------------
# 1. bgm / GGM / continuous / Bernoulli / listwise
# ---------------------------------------------------------------------------
x_cont <- make_continuous_data(n = 50, p = 5)
cm <- check_model(
  x = x_cont,
  variable_type   = "continuous",
  baseline_category = 0L,
  edge_selection  = TRUE,
  edge_prior      = "Bernoulli",
  inclusion_probability = 0.5
)
fixtures[["bgm_ggm_bernoulli_listwise"]] <- list(
  id   = "bgm_ggm_bernoulli_listwise",
  desc = "bgm / GGM continuous / Bernoulli / listwise / no NAs",
  type = "bgm",
  input = list(
    x = x_cont,
    variable_type = "continuous",
    baseline_category = 0L,
    na_action = "listwise",
    edge_selection = TRUE,
    edge_prior = "Bernoulli",
    inclusion_probability = 0.5
  ),
  check_model    = cm,
  reformat_data  = NULL  # GGM path doesn't call reformat_data()
)

# ---------------------------------------------------------------------------
# 2. bgm / GGM / continuous / Beta-Bernoulli / listwise / with NAs
# ---------------------------------------------------------------------------
x_cont_na <- inject_nas(make_continuous_data(n = 60, p = 4))
cm2 <- check_model(
  x = x_cont_na,
  variable_type   = "continuous",
  baseline_category = 0L,
  edge_selection  = TRUE,
  edge_prior      = "Beta-Bernoulli",
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  inclusion_probability = 0.5
)
fixtures[["bgm_ggm_betabern_listwise_na"]] <- list(
  id   = "bgm_ggm_betabern_listwise_na",
  desc = "bgm / GGM continuous / Beta-Bernoulli / listwise / with NAs",
  type = "bgm",
  input = list(
    x = x_cont_na,
    variable_type = "continuous",
    baseline_category = 0L,
    na_action = "listwise",
    edge_selection = TRUE,
    edge_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1
  ),
  check_model    = cm2,
  reformat_data  = NULL  # GGM path doesn't call reformat_data()
)

# ---------------------------------------------------------------------------
# 3. bgm / OMRF / ordinal / Bernoulli / listwise
# ---------------------------------------------------------------------------
x_ord <- make_ordinal_data(n = 50, p = 5, max_cat = 3)
cm3 <- check_model(
  x = x_ord,
  variable_type   = "ordinal",
  baseline_category = 0L,
  edge_selection  = TRUE,
  edge_prior      = "Bernoulli",
  inclusion_probability = 0.5
)
rd3 <- reformat_data(
  x = x_ord,
  na_action     = "listwise",
  variable_bool = cm3$variable_bool,
  baseline_category = cm3$baseline_category
)
fixtures[["bgm_omrf_ordinal_bernoulli_listwise"]] <- list(
  id   = "bgm_omrf_ordinal_bernoulli_listwise",
  desc = "bgm / OMRF ordinal / Bernoulli / listwise",
  type = "bgm",
  input = list(
    x = x_ord,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    edge_selection = TRUE,
    edge_prior = "Bernoulli",
    inclusion_probability = 0.5
  ),
  check_model    = cm3,
  reformat_data  = rd3
)

# ---------------------------------------------------------------------------
# 4. bgm / OMRF / ordinal / Beta-Bernoulli / impute / with NAs
# ---------------------------------------------------------------------------
x_ord_na <- inject_nas(make_ordinal_data(n = 60, p = 5, max_cat = 3))
cm4 <- check_model(
  x = x_ord_na,
  variable_type   = "ordinal",
  baseline_category = 0L,
  edge_selection  = TRUE,
  edge_prior      = "Beta-Bernoulli",
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  inclusion_probability = 0.5
)
# Impute path: need a fresh copy since reformat_data mutates x
x_ord_na_copy <- x_ord_na
rd4 <- reformat_data(
  x = x_ord_na_copy,
  na_action     = "impute",
  variable_bool = cm4$variable_bool,
  baseline_category = cm4$baseline_category
)
fixtures[["bgm_omrf_ordinal_betabern_impute_na"]] <- list(
  id   = "bgm_omrf_ordinal_betabern_impute_na",
  desc = "bgm / OMRF ordinal / Beta-Bernoulli / impute / with NAs",
  type = "bgm",
  input = list(
    x = x_ord_na,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "impute",
    edge_selection = TRUE,
    edge_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1
  ),
  check_model    = cm4,
  reformat_data  = rd4
)

# ---------------------------------------------------------------------------
# 5. bgm / OMRF / ordinal / SBM / listwise
# ---------------------------------------------------------------------------
x_ord5 <- make_ordinal_data(n = 50, p = 6, max_cat = 3)
cm5 <- check_model(
  x = x_ord5,
  variable_type   = "ordinal",
  baseline_category = 0L,
  edge_selection  = TRUE,
  edge_prior      = "Stochastic-Block",
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  beta_bernoulli_alpha_between = 1,
  beta_bernoulli_beta_between = 1,
  dirichlet_alpha = 1,
  lambda = 1,
  inclusion_probability = 0.5
)
rd5 <- reformat_data(
  x = x_ord5,
  na_action     = "listwise",
  variable_bool = cm5$variable_bool,
  baseline_category = cm5$baseline_category
)
fixtures[["bgm_omrf_ordinal_sbm_listwise"]] <- list(
  id   = "bgm_omrf_ordinal_sbm_listwise",
  desc = "bgm / OMRF ordinal / SBM / listwise",
  type = "bgm",
  input = list(
    x = x_ord5,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    edge_selection = TRUE,
    edge_prior = "Stochastic-Block",
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1,
    beta_bernoulli_alpha_between = 1,
    beta_bernoulli_beta_between = 1,
    dirichlet_alpha = 1,
    lambda = 1
  ),
  check_model    = cm5,
  reformat_data  = rd5
)

# ---------------------------------------------------------------------------
# 6. bgm / OMRF / blume-capel / Bernoulli / listwise / custom baseline
# ---------------------------------------------------------------------------
# BC variables with scores starting at 1 (not 0) — triggers recoding
x_bc <- make_ordinal_data(n = 50, p = 4, max_cat = 4)
x_bc <- x_bc + 1L  # shift to 1-based scores
vtype6 <- rep("blume-capel", 4)
cm6 <- check_model(
  x = x_bc,
  variable_type   = vtype6,
  baseline_category = 3L,
  edge_selection  = TRUE,
  edge_prior      = "Bernoulli",
  inclusion_probability = 0.5
)
rd6 <- reformat_data(
  x = x_bc,
  na_action     = "listwise",
  variable_bool = cm6$variable_bool,
  baseline_category = cm6$baseline_category
)
fixtures[["bgm_omrf_blumecapel_bernoulli_listwise"]] <- list(
  id   = "bgm_omrf_blumecapel_bernoulli_listwise",
  desc = "bgm / OMRF blume-capel / Bernoulli / listwise / custom baseline",
  type = "bgm",
  input = list(
    x = x_bc,
    variable_type = vtype6,
    baseline_category = 3L,
    na_action = "listwise",
    edge_selection = TRUE,
    edge_prior = "Bernoulli",
    inclusion_probability = 0.5
  ),
  check_model    = cm6,
  reformat_data  = rd6
)

# ---------------------------------------------------------------------------
# 7. bgm / OMRF / blume-capel / Beta-Bernoulli / impute / with NAs
# ---------------------------------------------------------------------------
x_bc_na <- inject_nas(x_bc)
cm7 <- check_model(
  x = x_bc_na,
  variable_type   = vtype6,
  baseline_category = 3L,
  edge_selection  = TRUE,
  edge_prior      = "Beta-Bernoulli",
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  inclusion_probability = 0.5
)
x_bc_na_copy <- x_bc_na
rd7 <- reformat_data(
  x = x_bc_na_copy,
  na_action     = "impute",
  variable_bool = cm7$variable_bool,
  baseline_category = cm7$baseline_category
)
fixtures[["bgm_omrf_blumecapel_betabern_impute_na"]] <- list(
  id   = "bgm_omrf_blumecapel_betabern_impute_na",
  desc = "bgm / OMRF blume-capel / Beta-Bernoulli / impute / with NAs",
  type = "bgm",
  input = list(
    x = x_bc_na,
    variable_type = vtype6,
    baseline_category = 3L,
    na_action = "impute",
    edge_selection = TRUE,
    edge_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1
  ),
  check_model    = cm7,
  reformat_data  = rd7
)

# ---------------------------------------------------------------------------
# 8. bgm / OMRF / mixed ordinal+BC / Bernoulli / listwise
# ---------------------------------------------------------------------------
x_mixed <- make_ordinal_data(n = 50, p = 5, max_cat = 4)
vtype8 <- c("ordinal", "ordinal", "blume-capel", "ordinal", "blume-capel")
bcat8 <- c(0L, 0L, 2L, 0L, 2L)
cm8 <- check_model(
  x = x_mixed,
  variable_type   = vtype8,
  baseline_category = bcat8,
  edge_selection  = TRUE,
  edge_prior      = "Bernoulli",
  inclusion_probability = 0.5
)
rd8 <- reformat_data(
  x = x_mixed,
  na_action     = "listwise",
  variable_bool = cm8$variable_bool,
  baseline_category = cm8$baseline_category
)
fixtures[["bgm_omrf_mixed_bernoulli_listwise"]] <- list(
  id   = "bgm_omrf_mixed_bernoulli_listwise",
  desc = "bgm / OMRF mixed ordinal+BC / Bernoulli / listwise",
  type = "bgm",
  input = list(
    x = x_mixed,
    variable_type = vtype8,
    baseline_category = bcat8,
    na_action = "listwise",
    edge_selection = TRUE,
    edge_prior = "Bernoulli",
    inclusion_probability = 0.5
  ),
  check_model    = cm8,
  reformat_data  = rd8
)

# ---------------------------------------------------------------------------
# 9. bgmCompare / ordinal / listwise / 2 groups
# ---------------------------------------------------------------------------
x_comp9 <- make_ordinal_data(n = 60, p = 4, max_cat = 2)
group9 <- rep(1:2, each = 30)
cm9 <- check_compare_model(
  x = x_comp9,
  y = NULL,
  group_indicator  = group9,
  difference_selection = TRUE,
  variable_type    = "ordinal",
  baseline_category = 0L,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
rd9 <- compare_reformat_data(
  x = cm9$x,
  group = cm9$group_indicator,
  na_action     = "listwise",
  variable_bool = cm9$variable_bool,
  baseline_category = cm9$baseline_category
)
fixtures[["compare_ordinal_listwise_2groups"]] <- list(
  id   = "compare_ordinal_listwise_2groups",
  desc = "bgmCompare / ordinal / listwise / 2 groups",
  type = "compare",
  input = list(
    x = x_comp9,
    group_indicator = group9,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5
  ),
  check_model    = cm9,
  reformat_data  = rd9
)

# ---------------------------------------------------------------------------
# 10. bgmCompare / ordinal / impute / 2 groups + NAs
# ---------------------------------------------------------------------------
x_comp10 <- inject_nas(make_ordinal_data(n = 60, p = 4, max_cat = 2))
group10 <- rep(1:2, each = 30)
cm10 <- check_compare_model(
  x = x_comp10,
  y = NULL,
  group_indicator  = group10,
  difference_selection = TRUE,
  variable_type    = "ordinal",
  baseline_category = 0L,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
x_comp10_copy <- cm10$x
rd10 <- compare_reformat_data(
  x = x_comp10_copy,
  group = cm10$group_indicator,
  na_action     = "impute",
  variable_bool = cm10$variable_bool,
  baseline_category = cm10$baseline_category
)
fixtures[["compare_ordinal_impute_2groups_na"]] <- list(
  id   = "compare_ordinal_impute_2groups_na",
  desc = "bgmCompare / ordinal / impute / 2 groups + NAs",
  type = "compare",
  input = list(
    x = x_comp10,
    group_indicator = group10,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "impute",
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5
  ),
  check_model    = cm10,
  reformat_data  = rd10
)

# ---------------------------------------------------------------------------
# 11. bgmCompare / blume-capel / listwise / 2 groups
# ---------------------------------------------------------------------------
x_comp11 <- make_ordinal_data(n = 60, p = 4, max_cat = 4)
group11 <- rep(1:2, each = 30)
vtype11 <- rep("blume-capel", 4)
cm11 <- check_compare_model(
  x = x_comp11,
  y = NULL,
  group_indicator  = group11,
  difference_selection = TRUE,
  variable_type    = vtype11,
  baseline_category = 2L,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
rd11 <- compare_reformat_data(
  x = cm11$x,
  group = cm11$group_indicator,
  na_action     = "listwise",
  variable_bool = cm11$variable_bool,
  baseline_category = cm11$baseline_category
)
fixtures[["compare_blumecapel_listwise_2groups"]] <- list(
  id   = "compare_blumecapel_listwise_2groups",
  desc = "bgmCompare / blume-capel / listwise / 2 groups",
  type = "compare",
  input = list(
    x = x_comp11,
    group_indicator = group11,
    variable_type = vtype11,
    baseline_category = 2L,
    na_action = "listwise",
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5
  ),
  check_model    = cm11,
  reformat_data  = rd11
)

# ---------------------------------------------------------------------------
# 12. bgmCompare / ordinal / listwise / >2 groups
# ---------------------------------------------------------------------------
x_comp12 <- make_ordinal_data(n = 90, p = 4, max_cat = 2)
group12 <- rep(1:3, each = 30)
cm12 <- check_compare_model(
  x = x_comp12,
  y = NULL,
  group_indicator  = group12,
  difference_selection = TRUE,
  variable_type    = "ordinal",
  baseline_category = 0L,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
rd12 <- compare_reformat_data(
  x = cm12$x,
  group = cm12$group_indicator,
  na_action     = "listwise",
  variable_bool = cm12$variable_bool,
  baseline_category = cm12$baseline_category
)
fixtures[["compare_ordinal_listwise_3groups"]] <- list(
  id   = "compare_ordinal_listwise_3groups",
  desc = "bgmCompare / ordinal / listwise / 3 groups",
  type = "compare",
  input = list(
    x = x_comp12,
    group_indicator = group12,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5
  ),
  check_model    = cm12,
  reformat_data  = rd12
)

# ---------------------------------------------------------------------------
# 13. bgmCompare / ordinal / listwise / categories missing in 1 group
# ---------------------------------------------------------------------------
mc_data <- make_missing_category_data()
cm13 <- check_compare_model(
  x = mc_data$x,
  y = NULL,
  group_indicator  = mc_data$group,
  difference_selection = TRUE,
  variable_type    = "ordinal",
  baseline_category = 0L,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
rd13 <- compare_reformat_data(
  x = cm13$x,
  group = cm13$group_indicator,
  na_action     = "listwise",
  variable_bool = cm13$variable_bool,
  baseline_category = cm13$baseline_category
)
fixtures[["compare_ordinal_missing_category"]] <- list(
  id   = "compare_ordinal_missing_category",
  desc = "bgmCompare / ordinal / listwise / category missing in group 2",
  type = "compare",
  input = list(
    x = mc_data$x,
    group_indicator = mc_data$group,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5
  ),
  check_model    = cm13,
  reformat_data  = rd13
)

# ---------------------------------------------------------------------------
# 14. bgmCompare / mixed ordinal+BC / listwise / missing categories + BC
# ---------------------------------------------------------------------------
x_comp14 <- mc_data$x
# Extend to 5 vars with an extra BC variable
x_extra <- matrix(sample(0:4, 60, replace = TRUE), nrow = 60, ncol = 1)
x_comp14 <- cbind(x_comp14, x_extra)
colnames(x_comp14) <- paste0("V", 1:5)
group14 <- mc_data$group
vtype14 <- c("ordinal", "ordinal", "ordinal", "ordinal", "blume-capel")
bcat14 <- c(0L, 0L, 0L, 0L, 2L)
cm14 <- check_compare_model(
  x = x_comp14,
  y = NULL,
  group_indicator  = group14,
  difference_selection = TRUE,
  variable_type    = vtype14,
  baseline_category = bcat14,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
rd14 <- compare_reformat_data(
  x = cm14$x,
  group = cm14$group_indicator,
  na_action     = "listwise",
  variable_bool = cm14$variable_bool,
  baseline_category = cm14$baseline_category
)
fixtures[["compare_mixed_missing_category_bc"]] <- list(
  id   = "compare_mixed_missing_category_bc",
  desc = "bgmCompare / mixed ord+BC / listwise / missing categories + BC",
  type = "compare",
  input = list(
    x = x_comp14,
    group_indicator = group14,
    variable_type = vtype14,
    baseline_category = bcat14,
    na_action = "listwise",
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5
  ),
  check_model    = cm14,
  reformat_data  = rd14
)

# ---------------------------------------------------------------------------
# 15. bgm / OMRF / ordinal / no edge selection
# ---------------------------------------------------------------------------
x_ord15 <- make_ordinal_data(n = 50, p = 4, max_cat = 2)
cm15 <- check_model(
  x = x_ord15,
  variable_type   = "ordinal",
  baseline_category = 0L,
  edge_selection  = FALSE,
  edge_prior      = "Bernoulli",
  inclusion_probability = 0.5
)
rd15 <- reformat_data(
  x = x_ord15,
  na_action     = "listwise",
  variable_bool = cm15$variable_bool,
  baseline_category = cm15$baseline_category
)
fixtures[["bgm_omrf_ordinal_no_edgesel"]] <- list(
  id   = "bgm_omrf_ordinal_no_edgesel",
  desc = "bgm / OMRF ordinal / no edge selection",
  type = "bgm",
  input = list(
    x = x_ord15,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    edge_selection = FALSE
  ),
  check_model    = cm15,
  reformat_data  = rd15
)

# ---------------------------------------------------------------------------
# 16. bgmCompare / ordinal / Beta-Bernoulli difference prior
# ---------------------------------------------------------------------------
x_comp16 <- make_ordinal_data(n = 60, p = 4, max_cat = 2)
group16 <- rep(1:2, each = 30)
cm16 <- check_compare_model(
  x = x_comp16,
  y = NULL,
  group_indicator  = group16,
  difference_selection = TRUE,
  variable_type    = "ordinal",
  baseline_category = 0L,
  difference_scale = 1,
  difference_prior = "Beta-Bernoulli",
  difference_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1
)
rd16 <- compare_reformat_data(
  x = cm16$x,
  group = cm16$group_indicator,
  na_action     = "listwise",
  variable_bool = cm16$variable_bool,
  baseline_category = cm16$baseline_category
)
fixtures[["compare_ordinal_betabern_diff"]] <- list(
  id   = "compare_ordinal_betabern_diff",
  desc = "bgmCompare / ordinal / Beta-Bernoulli difference prior",
  type = "compare",
  input = list(
    x = x_comp16,
    group_indicator = group16,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    difference_selection = TRUE,
    difference_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1
  ),
  check_model    = cm16,
  reformat_data  = rd16
)

# ---------------------------------------------------------------------------
# 17. bgmCompare / no difference selection
# ---------------------------------------------------------------------------
x_comp17 <- make_ordinal_data(n = 60, p = 4, max_cat = 2)
group17 <- rep(1:2, each = 30)
cm17 <- check_compare_model(
  x = x_comp17,
  y = NULL,
  group_indicator  = group17,
  difference_selection = FALSE,
  variable_type    = "ordinal",
  baseline_category = 0L,
  difference_scale = 1,
  difference_prior = "Bernoulli",
  difference_probability = 0.5
)
rd17 <- compare_reformat_data(
  x = cm17$x,
  group = cm17$group_indicator,
  na_action     = "listwise",
  variable_bool = cm17$variable_bool,
  baseline_category = cm17$baseline_category
)
fixtures[["compare_ordinal_no_diffsel"]] <- list(
  id   = "compare_ordinal_no_diffsel",
  desc = "bgmCompare / ordinal / no difference selection",
  type = "compare",
  input = list(
    x = x_comp17,
    group_indicator = group17,
    variable_type = "ordinal",
    baseline_category = 0L,
    na_action = "listwise",
    difference_selection = FALSE
  ),
  check_model    = cm17,
  reformat_data  = rd17
)

# ==============================================================================
# Save all fixtures
# ==============================================================================

manifest <- data.frame(
  id   = vapply(fixtures, `[[`, character(1), "id"),
  desc = vapply(fixtures, `[[`, character(1), "desc"),
  type = vapply(fixtures, `[[`, character(1), "type"),
  stringsAsFactors = FALSE
)

cat("\nSaving", nrow(manifest), "scaffolding fixtures to", fixture_dir, "\n\n")

for (f in fixtures) {
  path <- file.path(fixture_dir, paste0(f$id, ".rds"))
  saveRDS(f, path)
  cat(sprintf("  [%s] %s\n", f$id, f$desc))
}

saveRDS(manifest, file.path(fixture_dir, "manifest.rds"))
cat("\nManifest saved. Done.\n")
