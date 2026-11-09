# ==============================================================================
# Number of clusters under the Stochastic Block Model prior
# ==============================================================================
#
# The MFM-SBM prior carries two counts. B is the number of mixture components
# the prior makes available, B - 1 ~ Poisson(lambda). T is the number of
# clusters the variables occupy, the number of distinct labels in a
# partition, so T <= B. The responses inform the partition and hence T; they
# inform B only through T, because P(B = b | T = t, X) = P(B = b | T = t)
# (Miller & Harrison 2018, Eq. 3.7).
#
# This file computes the prior and posterior distributions of both counts over
# their full support, and their Bayes factors, from the saved allocations. For
# older fits without allocations it recovers T from the component summary,
# but only for the summary conventions of known bgms versions.
# ==============================================================================


# log p_B(b) for the shifted Poisson B - 1 ~ Poisson(lambda), the prior the
# sampler's partition coefficients (compute_Vn_mfm_sbm) use.
#' @importFrom stats dpois ppois
sbm_log_component_prior = function(b, lambda) {
  dpois(b - 1, lambda, log = TRUE)
}

log_sum_exp = function(x) {
  m = max(x)
  if(!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}


# ------------------------------------------------------------------
# sbm_log_partition_counts
# ------------------------------------------------------------------
# log S_eps(p, t), t = 1..p: the sum over partitions of p items into t
# blocks of prod_c eps^(|c|) (rising factorials), by the recursion
# S(n + 1, t) = eps S(n, t - 1) + (n + t eps) S(n, t), S(0, 0) = 1.
# At eps = 1 these are the unsigned Lah numbers.
# ------------------------------------------------------------------
sbm_log_partition_counts = function(p, eps) {
  log_s = c(0, rep(-Inf, p)) # S(0, t), t = 0..p
  for(n in seq_len(p) - 1L) {
    nxt = rep(-Inf, p + 1L)
    for(t in seq_len(n + 1L)) {
      nxt[t + 1L] = log_sum_exp(c(
        log(eps) + log_s[t],
        log(n + t * eps) + log_s[t + 1L]
      ))
    }
    log_s = nxt
  }
  log_s[-1L]
}


# ------------------------------------------------------------------
# sbm_count_kernel
# ------------------------------------------------------------------
# The prior side of the count inference, which depends only on p, the
# Dirichlet concentration eps and lambda:
#
#   log_vp[t]          log V_p(t) = log sum_b b_(t) / (b eps)^(p) p_B(b)
#   log_q[t]           log P(T = t) = log V_p(t) + log S_eps(p, t)
#   log_b_given_t[b,t] log P(B = b | T = t), b = 1..b_max
#
# b_max grows until the weights beyond it decay geometrically and their
# total, relative to V_p(t), is below `tol` for every t; that bound is
# returned as `omitted`.
# ------------------------------------------------------------------
sbm_count_kernel = function(num_variables, dirichlet_alpha, lambda,
                            tol = 1e-12) {
  p = num_variables
  eps = dirichlet_alpha
  b_max = max(p + 10, ceiling(lambda + 20 * sqrt(lambda) + 20))
  repeat {
    b = seq_len(b_max)
    base = sbm_log_component_prior(b, lambda) +
      lgamma(b * eps) - lgamma(b * eps + p)
    log_w = vapply(seq_len(p), function(t) {
      ifelse(b >= t, base + lgamma(b + 1) - lgamma(pmax(b - t, 0) + 1), -Inf)
    }, numeric(b_max))
    log_vp = apply(log_w, 2, log_sum_exp)
    # Successive weights beyond b_max shrink by at most
    # u = lambda / b_max * (b_max + 1) / (b_max + 1 - p): the Poisson ratio
    # times the largest falling-factorial ratio; the rising-factorial ratio
    # is below one.
    u = lambda / b_max * (b_max + 1) / (b_max + 1 - p)
    omitted = if(u < 1) exp(log_w[b_max, ] - log_vp) * u / (1 - u) else Inf
    if(all(omitted < tol)) break
    if(b_max > 1e6) {
      stop(
        "The prior on the number of components has too heavy a tail to ",
        "evaluate (lambda = ", lambda, ")."
      )
    }
    b_max = 2 * b_max
  }
  list(
    num_variables = p,
    dirichlet_alpha = eps,
    lambda = lambda,
    b_max = b_max,
    log_vp = log_vp,
    log_q = log_vp + sbm_log_partition_counts(p, eps),
    log_b_given_t = sweep(log_w, 2, log_vp),
    omitted = max(omitted)
  )
}


# ------------------------------------------------------------------
# sbm_summary_convention
# ------------------------------------------------------------------
# The prior on B that a fit's stored posterior_num_blocks assumed:
#   bgms 0.1.6.x            zero-truncated Poisson ("zero_truncated"),
#                           while the sampler used the shifted Poisson;
#   bgms 0.2.0.0 and later  shifted Poisson ("shifted").
# Both are restricted to B <= p and renormalized for each T. Versions
# before 0.1.6.0 sampled under a different prior on B: NA.
# ------------------------------------------------------------------
sbm_summary_convention = function(version) {
  if(is.null(version)) {
    return(NA_character_)
  }
  version = package_version(version)
  if(version < "0.1.6.0") {
    return(NA_character_)
  }
  if(version < "0.2.0.0") "zero_truncated" else "shifted"
}

# Q[k, t]: the legacy summary's P(B = k | T = t) on k <= p, so that the
# stored posterior_num_blocks equals Q %*% P(T | X).
sbm_legacy_summary_matrix = function(kernel, convention) {
  p = kernel$num_variables
  lw = kernel$log_b_given_t[seq_len(p), , drop = FALSE]
  if(identical(convention, "zero_truncated")) {
    # lambda^k / k! against lambda^(k - 1) / (k - 1)!
    lw = lw + log(kernel$lambda) - log(seq_len(p))
  }
  exp(sweep(lw, 2, apply(lw, 2, log_sum_exp)))
}

sbm_occupied_from_summary = function(reported, kernel, convention, version) {
  p = kernel$num_variables
  if(length(reported) != p || anyNA(reported) || any(reported < 0)) {
    stop("The stored number-of-blocks summary does not match the fit.")
  }
  Q = sbm_legacy_summary_matrix(kernel, convention)
  if(rcond(Q) < 1e-12) {
    stop(
      "The stored number-of-blocks summary is too ill-conditioned to invert. ",
      "Use a fit that carries its saved allocations."
    )
  }
  w = as.numeric(forwardsolve(Q, reported))
  slack = max(
    -min(c(w, 0)), abs(sum(w) - 1),
    max(abs(as.numeric(Q %*% pmax(w, 0)) - reported))
  )
  if(!is.finite(slack) || slack > 1e-6) {
    stop(
      "The stored number-of-blocks summary does not match the summary ",
      "convention of bgms ", version, " (discrepancy ",
      format(slack, digits = 2), "). Development builds of bgms 0.2.0.0 ",
      "from before July 2026 used the 0.1.6 convention. Compatibility ",
      "cannot be established from the summary alone; use a fit that ",
      "carries its saved allocations."
    )
  }
  w = pmax(w, 0)
  w / sum(w)
}


# ------------------------------------------------------------------
# sbm_fit_provenance
# ------------------------------------------------------------------
# What a fit's recorded version says about its SBM draws and its stored
# summary. `sampler` is "fixed" for fits from bgms >= 0.2.1.0, whose plain
# SBM update draws a new cluster's block probabilities from their
# conditional posterior; "affected" for earlier releases, which drew them
# from the prior; "unknown" without a version. Fits under the joint
# precision-graph prior use a separate update and are "not assessed".
# ------------------------------------------------------------------
sbm_fit_provenance = function(arguments) {
  version = arguments$version
  joint = identical(arguments$precision_graph_prior, "joint")
  sampler = if(is.null(version)) {
    "unknown"
  } else if(joint) {
    "not assessed"
  } else if(isTRUE(arguments$sbm_sampler_fix) ||
    package_version(version) >= "0.2.1.0") {
    "fixed"
  } else {
    "affected"
  }
  list(
    version = if(is.null(version)) NA_character_ else as.character(version),
    build_commit = if(is.null(arguments$build_commit)) NA_character_ else arguments$build_commit,
    sampler = sampler,
    summary_convention = sbm_summary_convention(version)
  )
}


# ------------------------------------------------------------------
# sbm_provenance_message
# ------------------------------------------------------------------
# One message per kind and session, matched to what is wrong:
#   affected  draws from a sampler with the new-cluster defect: refit;
#   unknown   no recorded version: compatibility cannot be established;
#   legacy    valid draws, but the caller is handed the stored component
#             summary (legacy_summary = TRUE): resummarize with
#             extract_sbm_counts().
# Joint precision-graph fits ("not assessed") get no message.
# ------------------------------------------------------------------
.sbm_messages_shown = new.env(parent = emptyenv())

sbm_provenance_message = function(arguments, legacy_summary = FALSE) {
  prov = sbm_fit_provenance(arguments)
  kind = if(prov$sampler == "affected") {
    "affected"
  } else if(prov$sampler == "unknown") {
    "unknown"
  } else if(legacy_summary && prov$sampler == "fixed") {
    "legacy"
  } else {
    return(invisible(NULL))
  }
  if(isTRUE(.sbm_messages_shown[[kind]])) {
    return(invisible(NULL))
  }
  text = switch(kind,
    affected = paste0(
      "This fit comes from bgms ", prov$version, ", whose stochastic block ",
      "model sampler drew the block probabilities of a newly opened cluster ",
      "from the prior instead of their conditional posterior (fixed in bgms ",
      "0.2.1.0). Its cluster allocations, cluster-count probabilities and ",
      "clustering Bayes factors can be biased, and so, to a lesser degree, can ",
      "its edge inclusion results. Refit the model with bgms 0.2.1.0 or later; ",
      "resummarizing cannot repair the draws."
    ),
    unknown = paste0(
      "This fit does not record the bgms version that produced it, so ",
      "whether its stochastic block model draws and summaries are compatible ",
      "with this version cannot be established automatically."
    ),
    legacy = paste0(
      "posterior_num_blocks is the stored summary of the number of components ",
      "B, restricted to B <= p and renormalized; it is kept unchanged for ",
      "compatibility. Use extract_sbm_counts() for the full posterior of B ",
      "and for the number of occupied clusters T."
    )
  )
  # Classed, so callers can muffle one kind (e.g. a package that reads the
  # counts through extract_sbm_counts() and never shows posterior_num_blocks).
  # Signalled as message() does, but counted as shown only once delivered, so
  # a muffled message still reaches a later direct call. knitr writes a
  # message into the document and then muffles it, so under knitr a muffled
  # message counts as shown too; otherwise every chunk would repeat it.
  cond = structure(
    class = c(
      paste0("bgms_sbm_", kind, "_message"), "bgms_sbm_provenance_message",
      "message", "condition"
    ),
    list(
      message = paste0(paste(strwrap(text, width = 78), collapse = "\n"), "\n"),
      call = NULL
    )
  )
  delivered = FALSE
  withRestarts(
    {
      signalCondition(cond)
      cat(conditionMessage(cond), file = stderr())
      delivered = TRUE
    },
    muffleMessage = function() NULL
  )
  if(delivered || isTRUE(getOption("knitr.in.progress"))) {
    assign(kind, TRUE, envir = .sbm_messages_shown)
  }
  invisible(NULL)
}


#' @title Number of Clusters Under the Stochastic Block Model Prior
#'
#' @description
#' Prior and posterior distributions of the two cluster counts of the
#' Stochastic Block Model (SBM) prior, from a fit with
#' `edge_prior = sbm_prior()` (or `difference_prior = sbm_prior()` in
#' [bgmCompare()]). \eqn{B} is the number of mixture components the prior
#' makes available; \eqn{T} is the number of clusters the variables actually
#' occupy, so \eqn{T \le B}. The data inform \eqn{T} directly and \eqn{B}
#' only through \eqn{T}.
#'
#' `extract_sbm_counts()` works from a fitted model; `sbm_counts()` from saved
#' allocations or a stored summary. Use [sbm_count_bf()] for Bayes factors
#' and `print()` for a side-by-side table.
#'
#' @details
#' The posterior of \eqn{T} is the share of saved allocation draws with
#' \eqn{t} distinct labels. Its prior is \eqn{P(T = t) = V_p(t)
#' S_\epsilon(p, t)}, where \eqn{V_p} are the partition coefficients of
#' the mixture of finite mixtures and \eqn{S_\epsilon(p, t)} sums the
#' Dirichlet weights over partitions with \eqn{t} blocks. The prior of
#' \eqn{B} is the shifted Poisson \eqn{B - 1 \sim \textrm{Poisson}(\lambda)},
#' so \eqn{E[B] = 1 + \lambda}, while \eqn{E[T]} also depends on \eqn{p} and
#' the Dirichlet concentration. The posterior of \eqn{B} is
#' \eqn{\sum_t P(B = b \mid T = t) P(T = t \mid X)}, evaluated over its full
#' support; the components table stops where both the prior and the
#' posterior mass beyond it fall below `tol`, and that remaining mass is
#' returned as `tail`.
#'
#' The number of clusters in a representative partition, such as the
#' posterior mean or mode allocations of a fit, need not equal the posterior
#' mode of \eqn{T} or of \eqn{B}.
#'
#' A fit without saved allocations is summarized from its stored
#' `posterior_num_blocks`, which is inverted under the summary convention of
#' the bgms version that produced it: bgms 0.1.6.x and 0.2.0.0 and later are
#' recognized. The summary itself cannot reveal which convention produced it,
#' so this rests on the recorded version; prefer saved allocations whenever
#' they exist. Fits from versions before 0.1.6.0 used a different prior on
#' the number of components and are refused.
#'
#' @param fit A fitted `bgms` or `bgmCompare` object with a Stochastic Block
#'   prior on its indicators.
#' @param allocations Saved cluster allocations: a matrix (draws x variables)
#'   or a list of such matrices, one per chain.
#' @param num_blocks The stored `posterior_num_blocks` of a fit without
#'   saved allocations (a vector, or a one-column matrix or data frame).
#'   Ignored when `allocations` is supplied.
#' @param lambda,dirichlet_alpha The fit's [sbm_prior()] settings.
#' @param version The bgms version the fit was produced with.
#' @param tol Largest prior or posterior mass of \eqn{B} that may lie beyond
#'   the last row of the components table.
#'
#' @return An object of class `"bgms_sbm_counts"`, a list with
#'   \describe{
#'     \item{`occupied`}{Data frame with `count` (\eqn{t = 1, \dots, p}),
#'       `prior`, `posterior` and its Monte Carlo standard error `mcse`
#'       (`NA` without allocations).}
#'     \item{`components`}{Data frame with `count`
#'       (\eqn{b = 1, \dots, b_{\max}}), `prior` and `posterior`.}
#'     \item{`tail`}{Prior and posterior mass of \eqn{B > b_{\max}}.}
#'     \item{`chain_counts`}{Matrix (chains x \eqn{p}) of the number of
#'       draws with each \eqn{T}, or `NULL` without allocations.}
#'     \item{`prior`}{The prior settings.}
#'     \item{`source`}{`"allocations"`, or the summary convention used.}
#'     \item{`provenance`}{The fit's bgms version and build commit, and
#'       what they imply for its draws (`extract_sbm_counts()` only).}
#'   }
#'
#' @seealso [sbm_count_bf()], [extract_sbm()], [sbm_prior()]
#' @family extractors
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], edge_prior = sbm_prior(), iter = 1000)
#' counts = extract_sbm_counts(fit)
#' counts
#' sbm_count_bf(counts, target = "occupied")
#' }
#' @export
extract_sbm_counts = function(fit, tol = 1e-10) {
  UseMethod("extract_sbm_counts")
}

#' @exportS3Method
#' @noRd
extract_sbm_counts.bgms = function(fit, tol = 1e-10) {
  arguments = extract_arguments(fit)
  if(!isTRUE(arguments$edge_selection) ||
    !identical(arguments$edge_prior, "Stochastic-Block")) {
    stop("Cluster counts need a fit with edge_selection = TRUE and edge_prior = sbm_prior().")
  }
  sbm_provenance_message(arguments)
  sbm_counts_from_fit(
    fit, arguments, arguments$lambda,
    arguments$dirichlet_alpha, tol
  )
}

#' @exportS3Method
#' @noRd
extract_sbm_counts.bgmCompare = function(fit, tol = 1e-10) {
  arguments = extract_arguments(fit)
  if(!isTRUE(arguments$difference_selection) ||
    !identical(arguments$difference_prior, "Stochastic-Block")) {
    stop("Cluster counts need a fit with difference_selection = TRUE and difference_prior = sbm_prior().")
  }
  sbm_provenance_message(arguments)
  sbm_counts_from_fit(
    fit, arguments, arguments$difference_lambda,
    arguments$difference_dirichlet_alpha, tol
  )
}

sbm_counts_from_fit = function(fit, arguments, lambda, dirichlet_alpha, tol) {
  allocations = tryCatch(fit$raw_samples$allocations, error = function(e) NULL)
  out = sbm_counts(
    allocations = allocations,
    num_blocks = fit$posterior_num_blocks,
    lambda = lambda,
    dirichlet_alpha = dirichlet_alpha,
    version = arguments$version,
    tol = tol
  )
  out$provenance = sbm_fit_provenance(arguments)
  out
}

#' @rdname extract_sbm_counts
#' @export
sbm_counts = function(allocations = NULL, num_blocks = NULL, lambda,
                      dirichlet_alpha, version, tol = 1e-10) {
  if(missing(version) || is.null(version)) {
    stop(
      "The bgms version that produced the fit is needed to know which prior ",
      "on the number of components it used; compatibility cannot be ",
      "established automatically."
    )
  }
  if(is.na(sbm_summary_convention(version))) {
    stop(
      "Fits from bgms ", version, " used a different prior on the number of ",
      "components than bgms 0.1.6.0 and later; their cluster counts cannot ",
      "be summarized here. Refit with the current version."
    )
  }
  positive = function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  if(!positive(lambda) || !positive(dirichlet_alpha)) {
    stop("'lambda' and 'dirichlet_alpha' must be positive numbers.")
  }

  chain_counts = NULL
  mcse = NULL
  ess = NULL
  if(!is.null(allocations)) {
    if(is.matrix(allocations)) allocations = list(allocations)
    p = ncol(allocations[[1L]])
    n_iter = vapply(allocations, nrow, integer(1))
    if(!all(vapply(allocations, ncol, integer(1)) == p) || any(n_iter < 1L)) {
      stop("The saved allocations have inconsistent dimensions.")
    }
    clusters = lapply(allocations, function(z) {
      apply(z, 1, function(row) length(unique(row)))
    })
    chain_counts = t(vapply(clusters, tabulate, integer(p), nbins = p))
    dimnames(chain_counts) = list(paste0("chain", seq_along(clusters)), seq_len(p))
    posterior_t = unname(colSums(chain_counts)) / sum(chain_counts)
    if(length(unique(n_iter)) == 1L) {
      indicator = array(0, dim = c(n_iter[1L], length(clusters), p))
      for(c in seq_along(clusters)) {
        indicator[, c, ] = outer(clusters[[c]], seq_len(p), "==") * 1
      }
      stats = .compute_indicator_ess_cpp(indicator)
      mcse = stats[, "mcse"]
      ess = stats[, "n_eff_mixt"]
    }
    source = "allocations"
  } else if(!is.null(num_blocks)) {
    reported = as.numeric(as.matrix(num_blocks)[, 1L])
    p = length(reported)
    source = sbm_summary_convention(version)
  } else {
    stop("Supply saved allocations or the stored number-of-blocks summary.")
  }

  kernel = sbm_count_kernel(p, dirichlet_alpha, lambda, tol = min(tol, 1e-12) / 10)
  note = NULL
  if(source != "allocations") {
    posterior_t = sbm_occupied_from_summary(reported, kernel, source, version)
    # The inversion cannot tell the two conventions apart (every column of
    # the summary matrix sums to one), so it rests on the version alone.
    if(package_version(version) == "0.2.0.0") {
      note = paste(
        "Recovered under the convention of the released bgms 0.2.0.0;",
        "development builds labelled 0.2.0.0 from before July 2026 used the",
        "0.1.6 convention instead."
      )
    }
  }
  prior_t = exp(kernel$log_q)

  b = seq_len(kernel$b_max)
  prior_b = exp(sbm_log_component_prior(b, lambda))
  posterior_b = as.numeric(exp(kernel$log_b_given_t) %*% posterior_t)
  tail_after = function(x) c(rev(cumsum(rev(x)))[-1L], 0)
  prior_tail = ppois(b - 1, lambda, lower.tail = FALSE)
  posterior_tail = tail_after(posterior_b) + kernel$omitted
  b_rows = which(prior_tail < tol & posterior_tail < tol)[1L]
  b_rows = max(p, if(is.na(b_rows)) kernel$b_max else b_rows)

  structure(
    list(
      occupied = data.frame(
        count = seq_len(p), prior = prior_t,
        posterior = posterior_t,
        mcse = if(is.null(mcse)) NA_real_ else mcse
      ),
      components = data.frame(
        count = seq_len(b_rows),
        prior = prior_b[seq_len(b_rows)],
        posterior = posterior_b[seq_len(b_rows)]
      ),
      tail = c(prior = prior_tail[b_rows], posterior = posterior_tail[b_rows]),
      chain_counts = chain_counts,
      ess = ess,
      prior = list(
        family = "shifted Poisson", lambda = lambda,
        dirichlet_alpha = dirichlet_alpha, num_variables = p
      ),
      source = source,
      note = note,
      provenance = list(version = as.character(version)),
      kernel = kernel
    ),
    class = "bgms_sbm_counts"
  )
}


#' @title Bayes Factors for the Number of Clusters
#'
#' @description
#' Bayes factors for the number of occupied clusters \eqn{T} or the number of
#' mixture components \eqn{B} of a Stochastic Block Model prior. With
#' `type = "complement"` the factor weighs more than one against exactly one;
#' with `type = "point"`, `count1` against `count2`. Values above 1 favour
#' the first hypothesis.
#'
#' @details
#' Each Bayes factor is the posterior odds divided by the prior odds of the
#' same two hypotheses: from \eqn{P(T = t)} for occupied clusters, and from
#' the shifted Poisson for components, whose prior odds of \eqn{B > 1}
#' against \eqn{B = 1} are \eqn{e^\lambda - 1}. Because \eqn{B > 1} still
#' allows every variable to sit in one cluster, the component factor for
#' \eqn{B > 1} against \eqn{B = 1} is at least \eqn{P(T = 1 \mid B > 1)}:
#' evidence for a single component is capped, evidence for a single occupied
#' cluster is not.
#'
#' When the chains never visited one of the two hypotheses the Bayes factor
#' is not estimable from the draws and is returned as `NA` with a bound,
#' based on an approximate 95% upper bound of \eqn{3 / n_{\mathrm{eff}}} for
#' the posterior probability of each unvisited count. The bound is a
#' heuristic: \eqn{n_{\mathrm{eff}}} is the smallest effective sample size
#' among the visited counts (or the number of draws when none is available),
#' which need not be that of the unvisited event. Run the chains longer
#' before reporting such a result. Per-chain Bayes factors are returned to
#' judge Monte Carlo error.
#'
#' @param x A `bgms_sbm_counts` object from [extract_sbm_counts()], or a fit
#'   it accepts.
#' @param target `"occupied"` for occupied clusters \eqn{T}, or
#'   `"components"` for mixture components \eqn{B}. There is no default.
#' @param type `"complement"` (more than one against exactly one) or
#'   `"point"`.
#' @param count1,count2 Whole numbers, required for `type = "point"`.
#'
#' @return An object of class `"bgms_sbm_count_bf"`, a list with the
#'   `target`, the `hypotheses`, `prior_odds`, `posterior_odds`, `bf`,
#'   `log_bf`, `status` (`"estimate"`, `"bound"` or `"not estimable"`),
#'   `bound` with its `bound_direction`, an `explanation`, and `per_chain`
#'   Bayes factors.
#'
#' @seealso [extract_sbm_counts()]
#' @family extractors
#' @export
sbm_count_bf = function(x, target, type = c("complement", "point"),
                        count1 = NULL, count2 = NULL) {
  if(!inherits(x, "bgms_sbm_counts")) x = extract_sbm_counts(x)
  if(missing(target)) {
    stop(
      "Choose target = \"occupied\" (T, the clusters the variables occupy) ",
      "or target = \"components\" (B, the components the prior makes available)."
    )
  }
  target = match.arg(target, c("occupied", "components"))
  type = match.arg(type)
  p = x$prior$num_variables
  kernel = x$kernel
  whole = function(v, hi) {
    is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 1 &&
      v == floor(v) && v <= hi
  }
  b_given_t = exp(kernel$log_b_given_t)

  # A hypothesis is a weight vector c over T: its posterior mass is
  # sum_t P(T = t | X) c(t), with c(t) = 1{t in H} for occupied clusters
  # and c(t) = P(B in H | T = t) for components.
  if(target == "occupied") {
    hi = p
    weights = function(set) as.numeric(seq_len(p) %in% set)
    prior_mass = function(set) sum(x$occupied$prior[seq_len(p) %in% set])
    sym = "T"
  } else {
    hi = kernel$b_max
    weights = function(set) colSums(b_given_t[set, , drop = FALSE])
    prior_mass = function(set) sum(exp(sbm_log_component_prior(set, x$prior$lambda)))
    sym = "B"
  }
  if(type == "complement") {
    h1 = seq.int(2L, hi)
    h2 = 1L
    labels = c(paste(sym, "> 1"), paste(sym, "= 1"))
    prior_odds = if(target == "components") {
      expm1(x$prior$lambda)
    } else {
      (1 - x$occupied$prior[1L]) / x$occupied$prior[1L]
    }
  } else {
    if(!whole(count1, hi) || !whole(count2, hi)) {
      stop("For type = \"point\", count1 and count2 must be whole numbers between 1 and ", hi, ".")
    }
    h1 = count1
    h2 = count2
    labels = c(paste(sym, "=", count1), paste(sym, "=", count2))
    prior_odds = prior_mass(h1) / prior_mass(h2)
  }
  c1 = weights(h1)
  c2 = weights(h2)
  if(type == "complement" && target == "components") c1 = 1 - c2

  bf_from = function(r) {
    m1 = sum(r * c1)
    m2 = sum(r * c2)
    if(m1 > 0 && m2 > 0) (m1 / m2) / prior_odds else NA_real_
  }
  r = x$occupied$posterior
  m1 = sum(r * c1)
  m2 = sum(r * c2)

  bound = NA_real_
  direction = NA_character_
  status = "estimate"
  explanation = NULL
  if(m1 > 0 && m2 > 0) {
    bf = (m1 / m2) / prior_odds
  } else {
    bf = NA_real_
    unvisited = r == 0
    ess = x$ess
    u = if(!is.null(ess) && any(!unvisited & r < 1 & is.finite(ess))) {
      3 / min(ess[!unvisited & r < 1 & is.finite(ess)])
    } else if(!is.null(x$chain_counts)) {
      3 / sum(x$chain_counts)
    } else {
      NA_real_
    }
    upper1 = u * sum(c1[unvisited])
    upper2 = u * sum(c2[unvisited])
    if(m1 == 0 && m2 == 0) {
      status = "not estimable"
      explanation = "The chains visited neither hypothesis."
    } else if(is.na(u)) {
      status = "not estimable"
      explanation = paste0(
        "The posterior gives no mass to ", labels[if(m1 == 0) 1L else 2L],
        ", and without saved draws the Monte Carlo error cannot be bounded."
      )
    } else if(m2 == 0) {
      status = "bound"
      direction = ">"
      bound = (m1 / upper2) / prior_odds
      explanation = paste0(
        "No saved draw is consistent with ", labels[2L], ". The bound assumes ",
        "the posterior probability of each unvisited count of occupied ",
        "clusters is below 3 / n_eff = ", format(u, digits = 2),
        "; run the chains longer before reporting it."
      )
    } else {
      status = "bound"
      direction = "<"
      bound = (upper1 / m2) / prior_odds
      explanation = paste0(
        "No saved draw is consistent with ", labels[1L], ". The bound assumes ",
        "the posterior probability of each unvisited count of occupied ",
        "clusters is below 3 / n_eff = ", format(u, digits = 2),
        "; run the chains longer before reporting it."
      )
    }
  }

  per_chain = if(!is.null(x$chain_counts)) {
    apply(x$chain_counts, 1, function(n) bf_from(n / sum(n)))
  } else {
    NULL
  }

  structure(
    list(
      target = target,
      type = type,
      hypotheses = labels,
      prior_odds = prior_odds,
      posterior_odds = if(m2 > 0) m1 / m2 else NA_real_,
      bf = bf,
      log_bf = log(bf),
      status = status,
      bound = bound,
      bound_direction = direction,
      explanation = explanation,
      per_chain = per_chain
    ),
    class = "bgms_sbm_count_bf"
  )
}


format_bf = function(x) {
  num = function(v) prettyNum(signif(v, 3), big.mark = ",")
  if(x$status == "estimate") {
    num(x$bf)
  } else if(x$status == "bound") {
    paste(if(x$bound_direction == ">") "at least" else "at most", num(x$bound))
  } else {
    "not estimable"
  }
}

#' @export
print.bgms_sbm_count_bf = function(x, ...) {
  what = if(x$target == "occupied") "occupied clusters (T)" else "components (B)"
  cat("Bayes factor for ", what, ": ", x$hypotheses[1L], " against ",
    x$hypotheses[2L], "\n",
    sep = ""
  )
  cat("  BF:", format_bf(x), "\n")
  cat("  prior odds:", format(x$prior_odds, digits = 4))
  if(!is.na(x$posterior_odds)) cat(", posterior odds:", format(x$posterior_odds, digits = 4))
  cat("\n")
  if(!is.null(x$explanation)) cat(paste0("  ", strwrap(x$explanation, width = 76)), sep = "\n")
  invisible(x)
}

#' @export
print.bgms_sbm_counts = function(x, digits = 3, ...) {
  occ = x$occupied
  comp = x$components
  # T runs to p, B can run beyond it: pad T with NA (printed blank).
  pad = function(v) c(v, rep(NA_real_, nrow(comp) - length(v)))
  probs = cbind(pad(occ$prior), pad(occ$posterior), comp$prior, comp$posterior)

  # Show counts up to the last one with a probability visible at `digits`;
  # the mass beyond goes into one ">k" row when it is visible.
  visible = 0.5 * 10^-digits
  shown = seq_len(max(which(apply(probs, 1, max, na.rm = TRUE) >= visible)))
  rest = colSums(probs[-shown, , drop = FALSE], na.rm = TRUE) +
    c(0, 0, x$tail[["prior"]], x$tail[["posterior"]])
  if(max(shown) >= nrow(occ)) rest[1:2] = NA
  fmt = function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  tab = data.frame(
    count = as.character(shown),
    prior_T = fmt(probs[shown, 1L]),
    posterior_T = fmt(probs[shown, 2L]),
    mcse_T = fmt(pad(occ$mcse)[shown]),
    prior_B = fmt(probs[shown, 3L]),
    posterior_B = fmt(probs[shown, 4L])
  )
  if(max(rest, na.rm = TRUE) >= visible) {
    tab = rbind(tab, data.frame(
      count = paste0(">", max(shown)), prior_T = fmt(rest[1L]),
      posterior_T = fmt(rest[2L]), mcse_T = "",
      prior_B = fmt(rest[3L]), posterior_B = fmt(rest[4L])
    ))
  }
  names(tab) = c("count", "prior T", "posterior T", "MCSE T", "prior B", "posterior B")
  if(all(tab[["MCSE T"]] == "")) tab[["MCSE T"]] = NULL

  cat("Number of clusters:\n")
  cat("  T: clusters the variables occupy.\n")
  cat("  B: mixture components the prior makes available (T <= B).\n")
  cat(sprintf(
    "  Prior: B - 1 ~ Poisson(%s), Dirichlet concentration %s.\n\n",
    format(x$prior$lambda), format(x$prior$dirichlet_alpha)
  ))
  print(tab, row.names = FALSE, right = TRUE)

  bf_t = sbm_count_bf(x, "occupied")
  bf_b = sbm_count_bf(x, "components")
  cat("\nBayes factors for more than one cluster against one:\n")
  cat("  T > 1 vs T = 1: ", format_bf(bf_t), "\n", sep = "")
  cat("  B > 1 vs B = 1: ", format_bf(bf_b), "\n", sep = "")
  if("bound" %in% c(bf_t$status, bf_b$status)) {
    cat("  A bound means the chains never visited one of the two hypotheses.\n")
  }
  if(x$source != "allocations") {
    cat("Recovered from the stored number-of-blocks summary; no MCSE available.\n")
    if(!is.null(x$note)) cat(strwrap(x$note, width = 76), sep = "\n")
  }
  cat("Use sbm_count_bf() for other comparisons.\n")
  invisible(x)
}
