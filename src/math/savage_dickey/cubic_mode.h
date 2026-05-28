#pragma once

#include <array>
#include <RcppArmadillo.h>

// Critical-point cubic for the L-space SD between-step kernel
//
//   f(phi) = -A phi^2 + B phi + (alpha - 1) log(s_jj + phi^2).
//
// Setting f'(phi) = 0 and clearing (s_jj + phi^2) gives the cubic
//
//   phi^3 + r phi^2 + (c_3 - s) phi + r c_3 = 0,
//   r = -B/(2A),  s = (alpha - 1)/A,  c_3 = s_jj.
//
// The cubic is solved by Cardano in depressed form (phi = t - r/3):
//
//   t^3 + p t + q = 0,
//   p = (c_3 - s) - r^2/3,
//   q = (r/3) * (2 r^2/9 + s + 2 c_3),
//   Delta = -4 p^3 - 27 q^2.
//
// Branches:
//   Delta >  0  ->  three distinct real roots; trigonometric form
//                   (requires p < 0, which Delta > 0 enforces).
//   Delta == 0  ->  repeated root; trigonometric form with clamp
//                   collapses two roots within float epsilon.
//   Delta <  0  ->  one real root; hyperbolic form (cosh/sinh by sign of p).
//
// Mode/saddle classification by sign of f''(phi):
//   f''(phi) = -2 A + 2 (alpha - 1) (c_3 - phi^2) / (c_3 + phi^2)^2.
//   f''(phi) <  0  ->  mode      (curvature = -f''(phi) > 0).
//   f''(phi) >  0  ->  saddle.
//
// At alpha = 1, s = 0 and the cubic factors as (phi^2 + c_3)(phi + r) = 0,
// returning a single real root phi = -r = B/(2A) exactly. The solver returns
// this without iteration.
//
// Status codes:
//   0  ok.
//   1  A <= 0 (PD-revert condition).
//   2  s_jj <= 0 (PD-revert condition).
//   3  numerical fallback (acos arg clamp beyond eps, or polynomial residual
//      check failed): caller should fall back to the alpha=1 Gaussian
//      reference (B/(2A), curvature 2A).

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
    int global_mode_index;              ///< roots[idx] is the global mode (largest ell among modes); -1 if n_modes == 0
    int local_mode_index;               ///< roots[idx] is the secondary mode; -1 if n_modes < 2
    int saddle_index;                   ///< roots[idx] is the saddle; -1 if n_modes < 2
    int status;                         ///< see header
};

CubicResult solve_sd_cubic(double A, double B, double s_jj, double alpha);

// Helper: evaluate the kernel and its second derivative at phi. Exposed for
// the AGHQ primitive and for unit tests.
double sd_log_kernel(double phi, double A, double B, double s_jj, double alpha);
double sd_log_kernel_pp(double phi, double A, double s_jj, double alpha);

}  // namespace savage_dickey
