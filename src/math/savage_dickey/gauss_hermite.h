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

namespace savage_dickey {

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

}  // namespace savage_dickey
