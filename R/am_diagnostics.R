# ==============================================================================
# Adaptive-Metropolis diagnostics
# ==============================================================================
#
# Post-sampling diagnostics for the componentwise adaptive-Metropolis
# sampler. Mirrors the role of nuts_diagnostics.R for NUTS, but
# necessarily lives in its own file because the underlying sampler
# organisation is different: NUTS does one accept/reject per iteration
# over the whole parameter vector, whereas adaptive-Metropolis does
# one accept/reject per parameter per iteration. The natural per-parameter
# acceptance summary therefore differs in shape (parameter x chain)
# from the per-trajectory NUTS summary (chain x iteration).
# ==============================================================================


# ------------------------------------------------------------------------------
# summarize_am_diagnostics
# ------------------------------------------------------------------------------
# Combine and summarize adaptive-Metropolis diagnostics across chains.
#
# The acceptance probability for a componentwise Metropolis update of a
# single parameter is estimated as the empirical *move rate*: the
# fraction of post-warmup iterations on which that parameter's value
# changed from the previous iteration. Since the sampler proposes from
# a continuous Normal proposal, "no change" is essentially equivalent
# to "rejected" (a probability-zero exact-match accept is ignored).
#
# @param out             List of chain outputs. Each element is a named list
#                        that must contain "main_samples" and
#                        "pairwise_samples" matrices (iterations x params).
# @param names_main      Character vector of main-effect parameter names.
# @param names_pairwise  Character vector of pairwise interaction names.
# @param target_accept   Target acceptance rate the sampler was tuned to
#                        (default: 0.44 for componentwise random-walk MH).
#
# Returns: An invisible named list with:
#   - accept_prob:   Numeric matrix (parameter x chain) of per-parameter
#       per-chain empirical Metropolis acceptance rates. Rows are
#       labelled by parameter (main effects first, then pairwise);
#       columns are labelled "chain 1", "chain 2", ...
#   - target_accept: Numeric scalar; the target acceptance rate.
#   - summary:       List with mean_accept_prob (mean across all
#       parameters and chains).
# ------------------------------------------------------------------------------
summarize_am_diagnostics = function(out, names_main, names_pairwise,
                                    target_accept = 0.44) {
  am_chains = Filter(function(chain) {
    all(c("main_samples", "pairwise_samples") %in% names(chain))
  }, out)

  if(length(am_chains) == 0) {
    return(NULL)
  }

  per_chain_accept = function(chain) {
    main_acc = apply(
      chain$main_samples, 2,
      function(col) mean(diff(col) != 0)
    )
    pair_acc = apply(
      chain$pairwise_samples, 2,
      function(col) mean(diff(col) != 0)
    )
    c(main_acc, pair_acc)
  }

  accept_prob_mat = vapply(
    am_chains, per_chain_accept,
    numeric(length(names_main) + length(names_pairwise))
  )

  rownames(accept_prob_mat) = c(names_main, names_pairwise)
  colnames(accept_prob_mat) = paste0("chain ", seq_along(am_chains))

  invisible(list(
    accept_prob = accept_prob_mat,
    target_accept = target_accept,
    summary = list(
      mean_accept_prob = mean(accept_prob_mat, na.rm = TRUE)
    )
  ))
}
