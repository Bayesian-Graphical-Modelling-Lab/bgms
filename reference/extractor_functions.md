# Extractor Functions for bgms Objects

Extractor Functions for bgms Objects

## Usage

``` r
extract_arguments(bgms_object)

# S3 method for class 'bgms'
extract_arguments(bgms_object)

# S3 method for class 'bgmCompare'
extract_arguments(bgms_object)

extract_indicators(bgms_object)

# S3 method for class 'bgms'
extract_indicators(bgms_object)

# S3 method for class 'bgmCompare'
extract_indicators(bgms_object)

extract_posterior_inclusion_probabilities(bgms_object)

# S3 method for class 'bgms'
extract_posterior_inclusion_probabilities(bgms_object)

extract_sbm(bgms_object)

# S3 method for class 'bgms'
extract_sbm(bgms_object)

# S3 method for class 'bgmCompare'
extract_posterior_inclusion_probabilities(bgms_object)

extract_indicator_priors(bgms_object)

# S3 method for class 'bgms'
extract_indicator_priors(bgms_object)

# S3 method for class 'bgmCompare'
extract_indicator_priors(bgms_object)

extract_pairwise_interactions(bgms_object)

# S3 method for class 'bgms'
extract_pairwise_interactions(bgms_object)

# S3 method for class 'bgmCompare'
extract_pairwise_interactions(bgms_object)

extract_category_thresholds(bgms_object)

# S3 method for class 'bgms'
extract_category_thresholds(bgms_object)

# S3 method for class 'bgmCompare'
extract_category_thresholds(bgms_object)

extract_group_params(bgms_object)

# S3 method for class 'bgmCompare'
extract_group_params(bgms_object)

extract_edge_indicators(bgms_object)

extract_pairwise_thresholds(bgms_object)

extract_rhat(bgms_object)

# S3 method for class 'bgms'
extract_rhat(bgms_object)

# S3 method for class 'bgmCompare'
extract_rhat(bgms_object)

extract_ess(bgms_object)

# S3 method for class 'bgms'
extract_ess(bgms_object)

# S3 method for class 'bgmCompare'
extract_ess(bgms_object)
```

## Arguments

- bgms_object:

  An object of class \`bgms\` or \`bgmCompare\`.

## Details

These functions extract various components from objects returned by the
\`bgm()\` function, such as edge indicators, posterior inclusion
probabilities, and parameter summaries.

Internally, indicator samples were stored in \`\$gamma\` (pre-0.1.4, now
defunct) and \`\$indicator\` (0.1.4–0.1.5, deprecated). As of \*\*bgms
0.1.6.0\*\*, they are stored in \`\$raw_samples\$indicator\`.

For `bgmCompare` objects, indicator samples were stored in
`$pairwise_difference_indicator` and `$main_difference_indicator`
(0.1.4–0.1.5, deprecated). As of \*\*bgms 0.1.6.0\*\*, they are stored
in `$raw_samples$indicator`.

Posterior inclusion probabilities are computed from edge indicators.

Internally, indicator samples were stored in \`\$gamma\` (pre-0.1.4, now
defunct) and \`\$indicator\` (0.1.4–0.1.5, deprecated). As of \*\*bgms
0.1.6.0\*\*, they are stored in \`\$raw_samples\$indicator\`.

Pairwise interactions were previously stored in \`\$pairwise_effects\`
(pre-0.1.4, now defunct) and \`\$posterior_mean_pairwise\` (0.1.4–0.1.5,
deprecated). As of \*\*bgms 0.1.6.0\*\*, they are stored in
\`\$raw_samples\$pairwise\` (raw samples) and
\`\$posterior_summary_pairwise\` (summaries).

For `bgmCompare` objects, pairwise interactions were stored in
`$interactions` (0.1.4–0.1.5, deprecated). As of \*\*bgms 0.1.6.0\*\*,
they are stored in `$raw_samples$pairwise`.

Category thresholds were previously stored in \`\$main_effects\`
(pre-0.1.4, now defunct) and \`\$posterior_mean_main\` (0.1.4–0.1.5,
deprecated). As of \*\*bgms 0.1.6.0\*\*, they are stored in
\`\$posterior_summary_main\`.

For `bgmCompare` objects, category thresholds were stored in
`$thresholds` (0.1.4–0.1.5, deprecated). As of \*\*bgms 0.1.6.0\*\*,
they are stored in `$raw_samples$main`.

## Functions

\- \`extract_arguments()\` – Extract model arguments -
\`extract_indicators()\` – Get sampled edge indicators -
\`extract_posterior_inclusion_probabilities()\` – Posterior edge
inclusion probabilities - \`extract_pairwise_interactions()\` –
Posterior mean of pairwise interactions -
\`extract_category_thresholds()\` – Posterior mean of category
thresholds - \`extract_indicator_priors()\` – Prior structure used for
edge indicators - \`extract_sbm()\` – Extract stochastic block model
parameters (if applicable) - \`extract_rhat()\` – Extract R-hat
convergence diagnostics - \`extract_ess()\` – Extract effective sample
size estimates
