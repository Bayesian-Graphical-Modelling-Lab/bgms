#pragma once

#include <RcppArmadillo.h>

/**
 * @file laplace.h
 * @brief Laplace approximation (with optional 1/n NLO correction) of the
 *        Savage-Dickey L-space conditional density at x_eval.
 *
 * Evaluates pi(x_eval | rest) from the log-kernel
 *   f(x) = -A x² + B x + (alpha - 1) log(s_jj + x²)
 * by locating the mode via Newton iteration (closed form at alpha = 1) and
 * applying Laplace's method; the NLO correction adds the Tierney-Kadane
 * f''''/f''' terms unless `nlo` is disabled. At alpha = 1 Laplace is exact
 * and the NLO contribution is identically zero.
 *
 * See the bgms manuscript / Savage-Dickey article for the derivation.
 *
 * Status codes (LSDResult.status):
 *   0  ok
 *   2  curvature non-positive at mode (Laplace invalid; chain should revert)
 *   3  Newton did not converge within newton_max_iter
 */

namespace savage_dickey {

struct LSDResult {
    double log_density;
    double x_mode;
    double curvature;     ///< -f''(x_mode)
    int    status;
};

LSDResult density_at_l_ji_one(
    double x_eval,
    double A,
    double B,
    double s_jj,
    double alpha,
    bool   nlo            = true,
    int    newton_max_iter = 50,
    double newton_tol      = 1e-10);

}  // namespace savage_dickey
