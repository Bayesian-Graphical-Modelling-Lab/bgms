#pragma once

#include <RcppArmadillo.h>

// Gauss-Hermite quadrature variant of the L-space SD primitive. Reliable
// across all chain configurations (no Laplace failure modes / bimodality
// issues), at the cost of ~64 log-kernel evaluations per call.
//
// Mathematically: with f(x) = -A x² + B x + (alpha-1) log(s + x²), the
// substitution y = sqrt(A) (x - B/(2A)) gives
//   ∫ exp(f(x)) dx = exp(B²/(4A)) / sqrt(A)
//                  · ∫ (s + (y/sqrt(A) + B/(2A))²)^(alpha-1) e^{-y²} dy
// and the inner integral is approximated by N-point Gauss-Hermite
// quadrature (physicists' convention, weight e^{-y²}).
//
// log_density returned is log pi(x_eval) = f(x_eval) - log_Z_quadrature.
//
// Use as a drop-in replacement (or fallback) for the Laplace+NLO primitive
// in cells where Newton/Laplace is unreliable.

namespace ggm_sd {

struct LSDQuadResult {
    double log_density;
    double log_Z;
    int    status;        // 0 ok; 1 invalid A (must be > 0)
};

LSDQuadResult density_at_l_ji_gh(double x_eval,
                                  double A, double B,
                                  double s_jj,
                                  double alpha,
                                  int num_nodes = 64);

// Adaptive Gauss-Hermite quadrature variant. Locates the actual mode of f
// via the closed-form cubic solver (sd_density_cubic.h), centres the
// nodes there with curvature kappa = -f''(phi*), and runs a three-level
// escalation cascade through the nested Genz-Keister sequence
// GK_3 -> GK_9 -> GK_35 until the |I_high - I_low| difference between
// adjacent levels falls below tol_strict. The GK nodes are nested
// (GK_3 nodes are a subset of GK_9, which is a subset of GK_35) and all
// weights are strictly positive, so the cascade's worst-case kernel-eval
// cost is one GK_35 evaluation (35 evals), not 3 + 9 + 35 = 47.
//
// At alpha = 1 the cubic returns phi* = B/(2A), kappa = 2A, and every
// level returns the same value to machine roundoff (the AGHQ integrand
// collapses to a constant in y_k).
//
// Polynomial exactness per level (verified against the analytic moments
// of exp(-y^2)):
//   GK_3:  degree 5
//   GK_9:  degree 15
//   GK_35: degree >= 50
//
// Escalation cascade (see density_at_l_ji_aghq impl):
//   1.  Compute I_GK3.
//   2.  Compute I_GK9.  If |I_GK9  - I_GK3| < tol_strict: return I_GK9.
//   3.  Compute I_GK35. If |I_GK35 - I_GK9| < tol_strict: return I_GK35.
//   4.  Otherwise return I_GK35 with status = 4 (not converged); the
//       caller should PD-revert.
//
// Bimodal cells (Delta >= 0 from the cubic) use the global mode as a single
// Laplace reference in this Phase-3a implementation; Phase 3b adds the
// mixture branch over both cubic modes, with each component running through
// the same per-mode cascade.
//
// Status codes:
//   0  ok (converged at some level <= GK_35).
//   1  A <= 0 (PD-revert condition).
//   2  s_jj <= 0 (PD-revert condition).
//   3  numerical fallback to the alpha=1 reference Gaussian (phi* = B/(2A),
//      kappa = 2A); the cubic solver returned no usable mode (e.g. triple-
//      root degeneracy) or curvature was below kappa_floor.
//   4  cascade exhausted without converging; |I_GK35 - I_GK9| > tol_strict
//      at the final step. Note: tol_strict is calibrated to the cheaper
//      rule's truncation error (GK_9 degree 15), so status = 4 means
//      "GK_9 was off by more than tol; the integrand has non-Gaussian
//      structure that the single-Laplace AGHQ at the global mode cannot
//      capture cleanly." Typically: bimodal cells where Phase 3b's
//      mixture branch is needed.

struct LSDAGHQResult {
    double log_density;    ///< log pi(x_eval | rest, Y)
    double log_Z;          ///< log of the normaliser at the final level
    double log_Z_err_est;  ///< |I_final - I_previous|: the cascade's
                           ///< empirical truncation estimate
    double x_mode;         ///< mode used as the Laplace centre
    double curvature;      ///< -f''(x_mode); positive
    int    n_nodes_used;   ///< 3, 9, or 35 (Genz-Keister level reached)
    int    status;
};

LSDAGHQResult density_at_l_ji_aghq(double x_eval,
                                    double A, double B,
                                    double s_jj,
                                    double alpha);

}  // namespace ggm_sd
