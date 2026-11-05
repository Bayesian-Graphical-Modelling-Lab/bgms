# an empty group-by-category cell marks its difference rows

    Code
      print(summary(fit))
    Output
      Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):
      
      groups: 1 = 1 (n = 80), 2 = 2 (n = 80)
      
      Category thresholds:
      parameter mean mcse sd n_eff Rhat
      A (1)
      A (2)
      A (3)
      B (1)
      B (2)
      C (1)
      ... (use `summary(fit)$main` to see full output)
      
      Pairwise interactions:
      parameter mean mcse sd n_eff Rhat
      A-B
      A-C
      B-C
      
      Inclusion probabilities:
      parameter mean mcse sd n_eff Rhat n0->1 n1->0
      A (main)
      A-B (pairwise)
      A-C (pairwise)
      B (main)
      B-C (pairwise)
      C (main)
      Note: blank mcse, n_eff and Rhat cells mark indicators whose inclusion
      probability was 0 or 1 in every draw; there is no variation to estimate them
      from.
      Main-effect difference rows are blank when main_difference_selection = FALSE
      left their indicators unsampled.
      
      Group differences (main effects):
      parameter mean mcse sd n_eff share_incl Rhat
      * A (diff1; 1)
      * A (diff1; 2)
      * A (diff1; 3)
      B (diff1; 1)
      B (diff1; 2)
      C (diff1; 1)
      ... (use `summary(fit)$main_diff` to see full output)
      * a group lacks observations in this category or in the reference category; the estimate reflects the prior, not the data
      Note: blank cells mark differences whose composite Monte Carlo error could
      not be formed because the difference was never selected or its inclusion
      probability was 0 or 1 in every draw.
      `summary(fit)$main_diff` retains the NA values.
      
      Group differences (pairwise effects):
      parameter mean mcse sd n_eff share_incl Rhat
      A-B (diff1)
      A-C (diff1)
      B-C (diff1)
      Note: blank cells mark differences whose composite Monte Carlo error could
      not be formed because the difference was never selected or its inclusion
      probability was 0 or 1 in every draw.
      `summary(fit)$pairwise_diff` retains the NA values.
      
      Use `summary(fit)$<component>` to access full results.
      See the `easybgm` package for other summary and plotting tools.

