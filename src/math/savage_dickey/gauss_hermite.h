#pragma once

#include <RcppArmadillo.h>

/**
 * @file gauss_hermite.h
 * @brief Gauss-Hermite quadrature of the Savage-Dickey L-space density.
 *
 * Drop-in alternative to the Laplace primitive (laplace.h) for cells where
 * Newton iteration is unreliable. Substitutes y = sqrt(A) (x - B/(2A)) so
 * the Gaussian weight e^{-y²} is exposed, then approximates the remaining
 * (s_jj + (y/sqrt(A) + B/(2A))²)^(alpha-1) factor by N-point Gauss-Hermite
 * quadrature (physicists' convention).
 *
 * See the bgms manuscript / Savage-Dickey article for the derivation.
 *
 * Status codes (LSDQuadResult.status):
 *   0  ok
 *   1  A <= 0 (must be positive for the substitution to be valid)
 */

namespace savage_dickey {

struct LSDQuadResult {
    double log_density;
    double log_Z;
    int    status;
};

LSDQuadResult density_at_l_ji_gh(double x_eval,
                                  double A, double B,
                                  double s_jj,
                                  double alpha,
                                  int num_nodes = 64);

}  // namespace savage_dickey
