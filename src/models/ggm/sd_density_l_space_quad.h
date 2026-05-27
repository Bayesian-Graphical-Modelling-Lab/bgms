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
// via the closed-form cubic solver (sd_density_cubic.h) and centres the GH
// nodes there with curvature kappa = -f''(phi*). At alpha = 1 the cubic
// returns phi* = B/(2A) and kappa = 2A; the AGHQ integrand becomes constant
// in y_k and the formula recovers the closed-form Gaussian normaliser at
// any N >= 1.
//
// At alpha > 1 with Delta >= 0 (bimodal cubic) this Phase-2 implementation
// uses only the global mode as the single Laplace reference; Phase 3 adds
// the mixture-AGHQ branch over both modes.
//
// Status codes:
//   0  ok.
//   1  A <= 0 (PD-revert condition).
//   2  s_jj <= 0 (PD-revert condition).
//   3  numerical fallback to the alpha=1 reference Gaussian (phi* = B/(2A),
//      kappa = 2A); the cubic solver returned no usable mode (e.g. triple-
//      root degeneracy) or curvature was non-positive at the returned mode.
//      log_Z is still computed via GH against the fallback Gaussian.

struct LSDAGHQResult {
    double log_density;    ///< log pi(x_eval | rest, Y)
    double log_Z;          ///< log of the normaliser
    double x_mode;         ///< mode used as the Laplace centre
    double curvature;      ///< -f''(x_mode); positive
    int    status;
};

LSDAGHQResult density_at_l_ji_aghq(double x_eval,
                                    double A, double B,
                                    double s_jj,
                                    double alpha,
                                    int num_nodes = 32);

}  // namespace ggm_sd
