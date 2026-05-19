#pragma once

#include <RcppArmadillo.h>

// Laplace + NLO closed-form approximations to log Z(G) for the
// spikeslab GGM prior with Gamma diagonal and determinant tilt:
//
//   K_ij ~ N(0, sigma^2) on off-diagonals (slab),
//   K_ii ~ Gamma(alpha, beta) on diagonals (rate parameterisation),
//   tilt |K|^delta with delta >= 0,
//   restricted to the PD cone supported by edge indicator matrix G.
//
// The full-recompute call costs O(q + |E| + non_edges_with_common_pred);
// the alpha = 1 incremental form costs O(deg(i)^2 + deg(i) * q) per
// single-edge toggle. Both reduce to the pre-tilt formula bit-exact at
// delta = 0.
//
// Port of:
//   ~/SV/Z/R/src/branchB_chain.cpp:470-620     (full-recompute)
//   ~/SV/Z/R/src/incremental_log_Z_NLO_gamma.h (alpha=1 incremental)
// validated SBC-clean in the z-project program update of 2026-05-17.

double log_Z_NLO_gamma(
    const arma::imat& G,
    double alpha, double beta, double sigma,
    bool include_F = false,
    double delta = 0.0);

// Toggle-endpoint reordering ("DEGORD"): relabel (i, j) to (0, 1) and
// permute all other vertices in their original order, then evaluate
// log_Z_NLO_gamma on the permuted graph. The closed-form is not
// permutation-invariant, so the chain's MH ratio must always evaluate
// both endpoint graphs in the same reordering.
double log_Z_NLO_gamma_degord(
    const arma::imat& G, int i, int j,
    double alpha, double beta, double sigma,
    bool include_F = false,
    double delta = 0.0);

// log Z_NLO(G_after) - log Z_NLO(G_before) under a single-edge toggle
// (i, j), at alpha = 1. Equivalent to the difference of two full-recompute
// calls to log_Z_NLO_gamma at alpha = 1 to machine precision, but evaluated
// via the §4.6 locality decomposition: only vertex min(i,j) contributes a
// change.
double log_Z_NLO_gamma_delta_incr_alpha1(
    const arma::imat& G_before, int i, int j,
    double beta, double sigma, double delta,
    bool include_F = false);

// log Z_NLO(G_after) - log Z_NLO(G_before) under a single-edge toggle
// (i, j), at general alpha (including alpha > 1).
//
// At alpha > 1 the alpha = 1 locality breaks: H_e gains a -4 beta (alpha - 1)
// / M_v[v] piece that depends on the larger endpoint, so toggling (i, j)
// cascades changes to H_e on every edge incident to i (forward AND backward).
//
// Implementation evaluates the partial sum
//
//   partial_log_Z_NLO_gamma_on_V(G, V_a)
//     = sum over log_Z_NLO_gamma term instances with at least one index
//       in V_a, with V_a = {i, j} ∪ N_G_before(i) (i, the toggled endpoint
//       j, and i's graph neighbours - the vertices whose H_e to/from i
//       cascade-changes under the toggle).
//
// then returns partial(G_after, V_a) - partial(G_before, V_a). Instances
// outside V_a are invariant under the toggle and cancel out.
//
// Reduces to log_Z_NLO_gamma_delta_incr_alpha1 at alpha = 1 to machine
// precision.
double log_Z_NLO_gamma_delta_incr_alphaN(
    const arma::imat& G_before, int i, int j,
    double alpha, double beta, double sigma, double delta,
    bool include_F = false);
