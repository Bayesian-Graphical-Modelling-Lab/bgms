#pragma once

#include <array>
#include <RcppArmadillo.h>

/**
 * @file cubic_mode.h
 * @brief Critical points of the Savage-Dickey L-space kernel.
 *
 * Solves f'(phi) = 0 for f(phi) = -A phi² + B phi + (alpha - 1) log(s_jj + phi²)
 * by reducing to a depressed cubic and dispatching to the trigonometric
 * branch (three real roots) or the hyperbolic branch (one real root) of
 * Cardano. Each root is returned with kernel value, second derivative, and
 * a mode-or-saddle label (mode iff f''(phi) < 0).
 *
 * See the bgms manuscript / Savage-Dickey article for the derivation.
 *
 * Status codes (CubicResult.status):
 *   0  ok
 *   1  A <= 0 (PD-revert condition)
 *   2  s_jj <= 0 (PD-revert condition)
 *   3  numerical fallback: caller should revert to the alpha = 1 Gaussian
 *      reference (mode B/(2A), curvature 2A)
 */

namespace savage_dickey {

struct CubicRoot {
    double phi;             ///< root value
    double ell;             ///< f(phi)
    double ell_pp;          ///< f''(phi); negative at modes, positive at saddle
    double curvature;       ///< -f''(phi); positive at modes, used by Laplace
    bool   is_mode;         ///< ell_pp < 0
};

struct CubicResult {
    int n_real_roots;                   ///< 1 or 3
    int n_modes;                        ///< 1 or 2 (3 only on numerical degeneracy)
    std::array<CubicRoot, 3> roots;     ///< valid entries 0 .. n_real_roots - 1
    int global_mode_index;              ///< largest ell among modes; -1 if n_modes == 0
    int local_mode_index;               ///< secondary mode; -1 if n_modes < 2
    int saddle_index;                   ///< -1 if n_modes < 2
    int status;
};

/// Solve the critical-point cubic; see file header for kernel and status codes.
CubicResult solve_sd_cubic(double A, double B, double s_jj, double alpha);

/// Evaluate the kernel at phi. Exposed for the quadrature primitives and unit tests.
double sd_log_kernel(double phi, double A, double B, double s_jj, double alpha);

/// Evaluate the second derivative of the kernel at phi.
double sd_log_kernel_pp(double phi, double A, double s_jj, double alpha);

}  // namespace savage_dickey
