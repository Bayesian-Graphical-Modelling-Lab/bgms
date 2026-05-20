#pragma once

#include <RcppArmadillo.h>
#include <utility>
#include <vector>

#include "models/ggm/degord_sampler.h"
#include "rng/rng_utils.h"

// Russian-Roulette V(Γ, U) estimator of 1/Z(Γ) for the hierarchical-spec
// GGM. Built on top of the Phase 2 DEGORD Bartlett-Cholesky sampler.
//
// Construction (Lyne 2015 / z-project Stage 3.2A):
//   V(Γ, U) = (1/c) · [1 + sum_{n=1..K_depth} (-1)^n · prod_{i=1..n} (Zhat_i - c)/c
//                                              / rho^n]
// where
//   Zhat_i  ~ log_Zhat_pi_from_pool(pools_t[i-1], pi_aux, chain_aux), expd
//   c       = kappa * exp(log_Z_NLO_degord(Γ))   -- analytic centring
//   K_depth ~ Geom(1 - rho)                       -- random truncation depth
//   rho     in (0, 1)                             -- truncation continuation prob
//
// E[V(Γ, U)] = 1/Z(Γ) exactly under (kappa, rho) chosen so the geometric
// series converges (radius + moment conditions; see route3a_V_helpers
// notes/exactness-proposition.md §2). The chain composes log|V_star/V_curr|
// into the between-edge MH ratio and tracks sign(V) separately for
// Lyne-style sign-corrected ergodic averaging.
//
// Port of:
//   ~/SV/Z/R/src/branchB_chain_route3a_degord.cpp:192-228
// validated on Z at q=10, delta=2 post the disable_log_r code-motion fix.

namespace degord {


// Per-Gamma V-estimator state, owned by the chain across between-step
// proposals. Updated lazily when the chain accepts a Γ-toggle.
struct ZRatioState {
    std::vector<arma::mat> pools_t;     // K_depth pools, each (dim x M_inner)
    int                    K_depth;     // Geom(1 - rho) draw
    double                 kappa;       // c = kappa * exp(log_Z_NLO)
    double                 rho;         // geometric truncation prob
    double                 log_Z_NLO_curr;   // analytic centring at Γ_curr
    int                    sign_curr;        // ±1, tracked from V(Γ_curr, U) sign
};


// V(Γ, U) at fixed (G_pi, pi_aux, c_val, rho). Returns the signed V.
//
// K_depth = 0 short-circuits to V = 1/c (the n=0 term only).
// Returns NaN if any Zhat_n is non-finite (caller treats as auto-reject).
double V_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const PiAux& pi_aux,
    const ChainAux& chain_aux,
    double c_val,
    double rho);


// Convenience overload: take G_pi instead of pre-built pi_aux (builds it
// internally). Used by R-callable test entry points; chain hot-path
// should pass the pre-built pi_aux to amortise.
double V_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const arma::imat& G_pi,
    const ChainAux& chain_aux,
    double c_val,
    double rho);


// Log-space variant of V(Γ, U). Factors V = (1/c) · S with
//
//   S = 1 + sum_{n=1..K_depth} (-1/rho)^n · prod_{m=1..n} (Zhat_m - c) / c
//     = 1 + sum_n term_n,
//
//   log|term_n| = -n · log(rho) + sum_{m=1..n} log|expm1(log_Zhat_m - log_c)|
//   sign(term_n) = (-1)^n · prod_m sign(expm1(log_Zhat_m - log_c))
//
// Computes log|S| via signed log-sum-exp over the K_depth + 1 terms, then
// returns (log|V|, sign(V)) = (-log_c + log|S|, sign(S)).
//
// The linear form V_at_Gamma_pi_degord exponentiates log_Z_NLO into c, which
// underflows to 0 at large p (e.g. log_Z_NLO ~ -3500 at p = 100 makes c = 0),
// silently breaking the MH ratio. The log-space form never materialises c.
//
// Returns {NaN, 0} on auto-reject (non-finite Zhat, S evaluates to 0, or
// any other non-finite intermediate). Caller treats as ln_alpha = -Inf.
//
// In the MH ratio, the log_kappa term cancels:
//   log|V_curr| - log|V_star|
//     = (log_Z_NLO_star - log_Z_NLO_curr) + (log|S_curr| - log|S_star|).
// Callers pass log_c = log(kappa) + log_Z_NLO at the relevant Gamma.
std::pair<double, int> V_log_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const PiAux& pi_aux,
    const ChainAux& chain_aux,
    double log_c,
    double rho);


// Convenience overload mirroring the linear form: take G_pi instead of a
// pre-built pi_aux. Used by R-callable test entry points.
std::pair<double, int> V_log_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const arma::imat& G_pi,
    const ChainAux& chain_aux,
    double log_c,
    double rho);


// Paired (log|V|, sign) for Gamma_curr / Gamma_star, computed with
// within-pool cache reuse. For each of K_depth pools, the inner loop
// builds Phi-state under a_curr ONCE (per log_Zhat_pi_from_pool_cache),
// then re-evaluates only row q-2 under a_star via row_qm2_logw_from_S
// using the cached (rw_head, S_trail). This halves the per-pool Phi
// rebuild cost vs two separate V_log_at_Gamma_pi_degord calls.
//
// G_pi_curr and G_pi_star must share q and may differ only at the
// trailing slot (q - 2, q - 1) — the DEGORD convention enforced by
// callers via degord_permutation(q, i, j). Returns {NaN, 0} on
// either side if any Zhat_n is non-finite, the signed sum collapses to
// zero, or log_c is non-finite.
struct LogSignedVPair {
    std::pair<double, int> curr;
    std::pair<double, int> star;
};


LogSignedVPair V_log_pair_at_Gamma_curr_star_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const PiAux& a_curr,
    const PiAux& a_star,
    const ChainAux& chain_aux,
    double log_c_curr,
    double log_c_star,
    double rho);


// Convenience overload mirroring the single-graph variant: take G_pi_curr
// and G_pi_star instead of pre-built PiAux'es. Used by R-callable test
// entry points; chain hot-path should pass pre-built pi_aux to amortise.
LogSignedVPair V_log_pair_at_Gamma_curr_star_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const arma::imat& G_pi_curr,
    const arma::imat& G_pi_star,
    const ChainAux& chain_aux,
    double log_c_curr,
    double log_c_star,
    double rho);


// Fresh U-pool draw: K_depth ~ Geom(1 - rho); pools_t[n] is (dim x M_inner)
// with iid N(0, 1) entries. Uses SafeRNG so chain seeds remain
// deterministic across platforms.
//
// Each pool is in pre-transposed (dim x M_inner) layout so the inner
// DEGORD kernel can access each sample's noise as a contiguous column.
void draw_U_degord_rr(
    SafeRNG& rng,
    int& K_depth,
    std::vector<arma::mat>& pools_t,
    int M_inner,
    int q,
    double rho);


}  // namespace degord
