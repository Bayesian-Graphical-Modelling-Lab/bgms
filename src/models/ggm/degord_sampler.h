#pragma once

#include <RcppArmadillo.h>
#include <vector>

#include "rng/rng_utils.h"

// DEGORD-permuted Bartlett-Cholesky importance sampler for log Zhat(G).
//
// Port of ~/SV/Z/R/src/degord_sampler.h (v4 layout, 2026-05-18). See the
// z header for the derivation; this header preserves the v4 architecture
// (per-nu transcendental tables, pre-transposed noise pool, trimmed
// PiAux) and naming.
//
// Convention: the DEGORD permutation in *this* file sends the toggle
// endpoints (i, j) to positions (q - 2, q - 1), so the toggled edge sits
// at the trailing diagonal-plus-off-diagonal pair. This is different
// from log_Z_NLO_gamma_degord (log_z_nlo.h), which places the toggle at
// (0, 1) for the closed-form formula's vertex-ordered structure. Both
// conventions are correct in their respective uses: the analytical c
// in the V estimator and the inner Zhat-from-pool importance sampler
// are computed under their own permutations of the same G.

namespace degord {

// ----------------------------------------------------------------------
// Per-chain auxiliary constants.
//
// Built once at the start of a chain (or when (alpha, beta, sigma,
// delta) changes). Holds all (q, alpha, beta, sigma, delta)-dependent
// constants plus per-nu transcendental tables indexed by forward degree
// nu in 0..q. Tables are read inside the inner kernel via
// c.nu_X[a.nu_pi[r]], avoiding per-row arma::vec allocations.
// ----------------------------------------------------------------------
struct ChainAux {
  int    q;
  double alpha;
  double beta;
  double sigma;
  double delta;
  double sigma2;
  double inv_sigma2;
  double two_beta;
  double inv_two_beta;
  double log_sigma;
  double log_2pi;
  double slab_kernel_const;
  double sigma_diag;
  double inv_sigma_diag2;
  double half_log_2pi_sigma_diag2;
  // Off-diagonal importance distribution selector:
  //   0 (default) - sample phi_rs around the unshifted saddle (tau-based).
  //   1           - sample phi_rs around the slab-tilt-shifted saddle.
  int    slab_tilt_mode;
  // Per-nu transcendental tables. Indexed by forward-degree nu in 0..q.
  std::vector<double> nu_chi_df;
  std::vector<double> nu_half_k;
  std::vector<double> nu_km1;
  std::vector<double> nu_lgamma_half_k;
  std::vector<double> nu_mu_l;
  std::vector<double> nu_diag_const;
  std::vector<double> nu_H_e_saddle;
  std::vector<double> nu_inv_sqrt_H_e_saddle;
  std::vector<double> nu_mu_coef_saddle;
  std::vector<double> nu_sigma_star_sq_saddle;
  std::vector<double> nu_inv_sigma_star_sq_saddle;
  std::vector<double> nu_slab_const_saddle;
  std::vector<double> nu_per_vertex;
};


// Build a ChainAux from (q, alpha, beta, sigma, delta). slab_tilt_mode is
// initialised to 0 (caller may overwrite after construction).
ChainAux make_chain_aux(int q, double alpha, double beta,
                        double sigma, double delta);


// ----------------------------------------------------------------------
// Per-permutation auxiliary state.
//
// Trimmed v4 layout: stores only the permuted graph, per-row forward
// degree, edge count, and log_C0. Per-row constants are read from
// ChainAux's nu_X[] tables.
// ----------------------------------------------------------------------
struct PiAux {
  int                 q;
  arma::imat          G_pi;
  std::vector<int>    nu_pi;
  int                 E_count;
  double              log_C0;
};


PiAux make_pi_aux(const arma::imat& G_pi, const ChainAux& c);


// ----------------------------------------------------------------------
// Graph permutation helpers.
// ----------------------------------------------------------------------

// Apply a vertex permutation pi to G. pi[u] gives the new index of u.
arma::imat permute_graph(const arma::imat& G, const arma::ivec& pi);

// DEGORD permutation that sends (i, j) -> (q-2, q-1), with all other
// vertices keeping their original order in positions 0..q-3.
arma::ivec degord_permutation(int q, int i, int j);


// ----------------------------------------------------------------------
// Inner kernel.
//
// Given a length-(q + q(q-1)/2) noise vector (Gaussian deviates), fill
// the upper-triangular Bartlett-Cholesky factor Phi_pi and return the
// total log-importance-weight for this sample. row_logw is filled with
// the per-row contributions for caching.
//
// Layout of noise[]:
//   noise[0 .. q-1]                                : diagonal innovations
//                                                    (Phi(r, r) for r = 0..q-1).
//   noise[q + edge_offset(r, s)] for r < s         : off-diagonal innovation
//                                                    for the (r, s) slot.
//     edge_offset(r, s) = r*(q-1) - r*(r-1)/2 + (s - r - 1)
//
// Off-diagonal slots are allocated for all (r, s) with r < s. Slots
// corresponding to non-edges are still consumed (filled with the slaving
// Phi(r, s) = -S_rs / Phi(r, r)) so the noise indexing stays stable.
// ----------------------------------------------------------------------
double phi_pi_sample_from_noise(
    arma::mat& Phi_pi,
    arma::vec& row_logw,
    const double* noise,
    const PiAux& a,
    const ChainAux& c);


// ----------------------------------------------------------------------
// log Zhat(G_pi) from a pre-transposed noise pool of shape (dim x M),
// column-major (so each sample's noise vector is a contiguous column).
//
//   dim = q + q*(q - 1) / 2
//   M   = number of inner samples
//
// Returns log Zhat (= log_C0 + log mean(exp(log_w))). Returns -Inf if
// every sample's log_w is non-finite.
// ----------------------------------------------------------------------
double log_Zhat_pi_from_pool(
    const arma::mat& noise_pool_t,
    const PiAux& a,
    const ChainAux& c);


// ----------------------------------------------------------------------
// PoolCache: per-sample state from a log_Zhat_pi_from_pool_cache call
// that delta_log_Zhat_pi_toggle reuses to avoid recomputing the head
// rows (0..q-3) when only the trailing edge (q-2, q-1) toggles.
// ----------------------------------------------------------------------
struct PoolCache {
  arma::vec log_w;      // per-sample total log-importance weight
  arma::vec rw_head;    // per-sample sum of row_logw[0..q-3]
  arma::vec S_trail;    // per-sample dot product of columns q-2 and q-1 over rows 0..q-3
};


double log_Zhat_pi_from_pool_cache(
    const arma::mat& noise_pool_t,
    const PiAux& a,
    const ChainAux& c,
    PoolCache& cache);


// Re-evaluate row (q-2)'s log_w (diag + slab/slaved at the trailing edge)
// for a *new* edge indicator G_pi(q-2, q-1) given the cached S_trail.
//
// z_qm2 is the diagonal innovation at row q-2 (noise[q-2]); z_trail is
// the off-diagonal innovation at the (q-2, q-1) slot (only used when
// slab_tilt_mode == 1).
double row_qm2_logw_from_S(
    double z_qm2,
    double z_trail,
    double S_trail,
    const PiAux& a,
    const ChainAux& c);


// ----------------------------------------------------------------------
// Efficient delta: log Zhat(Gamma_star) - log Zhat(Gamma_curr) under a
// single-edge toggle (i, j), with G_pi_star differing from G_pi_curr
// only at the trailing slot (q-2, q-1).
//
// Takes BOTH pool views:
//   noise_pool   : M x dim, column-major. z_qm2 = noise_pool(s, q-2) is
//                  contiguous along s.
//   noise_pool_t : dim x M, column-major. The kernel's per-sample
//                  contiguous noise extraction. Chain wrappers maintain
//                  both views together.
// ----------------------------------------------------------------------
double delta_log_Zhat_pi_toggle(
    const arma::mat& noise_pool,
    const arma::mat& noise_pool_t,
    const arma::imat& G_curr,
    int i, int j,
    const ChainAux& c);


// ----------------------------------------------------------------------
// Per-sample noise dimension for a chain of size q.
// ----------------------------------------------------------------------
inline int bartlett_pool_dim(int q) {
  return q + q * (q - 1) / 2;
}


// ----------------------------------------------------------------------
// Draw an independent standard-normal Bartlett pool of shape
// (dim x M_inner), column-major. Each column is one inner sample's
// noise vector. Uses SafeRNG so chain seeds are deterministic.
// ----------------------------------------------------------------------
arma::mat draw_bartlett_pool(SafeRNG& rng, int q, int M_inner);


}  // namespace degord
