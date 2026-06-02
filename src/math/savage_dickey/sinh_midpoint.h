#pragma once

#include <RcppArmadillo.h>

/**
 * @file sinh_midpoint.h
 * @brief Sinh-substitution + midpoint-rule quadrature of the Savage-Dickey
 *        L-space normaliser.
 *
 * Computes the normaliser of the conditional density
 *   I(A, B, s_jj, alpha) = integral over R of
 *                          (s_jj + phi²)^(alpha-1) exp(-A phi² + B phi) dphi
 * via the substitution phi = sqrt(s_jj) * sinh(t). The substitution pins
 * the branch points of (s_jj + phi²)^(alpha-1) at fixed imaginary distance
 * pi/2 from the real axis in t, so midpoint rule on a truncated interval
 * centred at the Gaussian-envelope mode converges uniformly across the
 * (A, B, s_jj, alpha) plane.
 *
 * See the bgms manuscript / Savage-Dickey article for the derivation.
 *
 * Status codes (LSDSinhResult.status):
 *   0  ok
 *   1  A <= 0 (PD-revert condition)
 *   2  s_jj <= 0 (PD-revert condition)
 */

namespace savage_dickey {

struct LSDSinhResult {
    double log_density;  ///< log pi(x_eval | rest, Y)
    double log_Z;        ///< log of the normaliser
    int    status;
};

LSDSinhResult density_at_l_ji_sinh(double x_eval,
                                    double A, double B,
                                    double s_jj,
                                    double alpha,
                                    int num_nodes = 128);

}  // namespace savage_dickey
