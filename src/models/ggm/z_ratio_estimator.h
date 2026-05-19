#pragma once

#include <RcppArmadillo.h>
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
