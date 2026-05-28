#include "cubic_mode.h"

#include <algorithm>
#include <cmath>

namespace savage_dickey {

double sd_log_kernel(double phi, double A, double B,
                     double s_jj, double alpha) {
    const double a1 = alpha - 1.0;
    double v = -A * phi * phi + B * phi;
    if (a1 != 0.0) v += a1 * std::log(s_jj + phi * phi);
    return v;
}

double sd_log_kernel_pp(double phi, double A, double s_jj, double alpha) {
    const double a1 = alpha - 1.0;
    double v = -2.0 * A;
    if (a1 != 0.0) {
        const double q  = s_jj + phi * phi;
        const double q2 = q * q;
        v += 2.0 * a1 * (s_jj - phi * phi) / q2;
    }
    return v;
}

namespace {

// Polynomial residual r(phi) = phi^3 + b phi^2 + c phi + d.
// Used as a self-check at each returned root.
inline double cubic_residual(double phi, double b, double c, double d) {
    return ((phi + b) * phi + c) * phi + d;
}

// Fill ell, ell_pp, curvature, is_mode at a root in-place.
inline void fill_root_metadata(CubicRoot& root, double A, double B,
                               double s_jj, double alpha) {
    root.ell       = sd_log_kernel    (root.phi, A, B, s_jj, alpha);
    root.ell_pp    = sd_log_kernel_pp (root.phi, A,    s_jj, alpha);
    root.curvature = -root.ell_pp;
    root.is_mode   = (root.ell_pp < 0.0);
}

// Classify mode/saddle indices on a 3-root result. Sets global_mode_index,
// local_mode_index, saddle_index, and n_modes. Caller has already filled
// per-root metadata.
inline void classify_modes_three(CubicResult& out) {
    int modes[3] = {-1, -1, -1};
    int n_modes = 0;
    int saddle_idx = -1;
    for (int k = 0; k < 3; ++k) {
        if (out.roots[k].is_mode) {
            modes[n_modes++] = k;
        } else {
            saddle_idx = k;
        }
    }
    out.n_modes = n_modes;
    if (n_modes == 2) {
        const int g = (out.roots[modes[0]].ell >= out.roots[modes[1]].ell)
                    ? modes[0] : modes[1];
        const int l = (g == modes[0]) ? modes[1] : modes[0];
        out.global_mode_index = g;
        out.local_mode_index  = l;
        out.saddle_index      = saddle_idx;
    } else if (n_modes == 1) {
        // Numerical degeneracy: the saddle and the other mode collided.
        // Return the single classifiable mode as global.
        out.global_mode_index = modes[0];
        out.local_mode_index  = -1;
        out.saddle_index      = saddle_idx;
    } else if (n_modes == 3) {
        // Mathematically impossible; pick top two by ell defensively.
        int order[3] = {0, 1, 2};
        std::sort(order, order + 3, [&](int a, int b){
            return out.roots[a].ell > out.roots[b].ell;
        });
        out.global_mode_index = order[0];
        out.local_mode_index  = order[1];
        out.saddle_index      = order[2];
        out.n_modes           = 2;
    } else {
        // n_modes == 0: all roots classified as saddle. Numerical pathology.
        out.global_mode_index = -1;
        out.local_mode_index  = -1;
        out.saddle_index      = -1;
    }
}

}  // namespace

CubicResult solve_sd_cubic(double A, double B, double s_jj, double alpha) {
    CubicResult out;
    out.n_real_roots      = 0;
    out.n_modes           = 0;
    out.global_mode_index = -1;
    out.local_mode_index  = -1;
    out.saddle_index      = -1;
    out.status            = 0;
    for (auto& root : out.roots) {
        root.phi        = 0.0;
        root.ell        = 0.0;
        root.ell_pp     = 0.0;
        root.curvature  = 0.0;
        root.is_mode    = false;
    }

    if (!(A > 0.0))    { out.status = 1; return out; }
    if (!(s_jj > 0.0)) { out.status = 2; return out; }

    const double r   = -B / (2.0 * A);
    const double s   = (alpha - 1.0) / A;
    const double c_3 = s_jj;

    // Coefficients of phi^3 + b phi^2 + c phi + d = 0 (for residual checks).
    const double b_coef = r;
    const double c_coef = c_3 - s;
    const double d_coef = r * c_3;

    // alpha = 1 fast path: cubic factors as (phi^2 + c_3)(phi + r) = 0;
    // single real root phi = -r = B/(2A), no iteration, no roundoff.
    if (alpha == 1.0) {
        out.n_real_roots = 1;
        out.roots[0].phi = -r;
        fill_root_metadata(out.roots[0], A, B, s_jj, alpha);
        out.n_modes           = out.roots[0].is_mode ? 1 : 0;
        out.global_mode_index = out.roots[0].is_mode ? 0  : -1;
        return out;
    }

    // Depressed cubic phi = t - r/3, then t^3 + p t + q = 0.
    const double p = c_coef - (r * r) / 3.0;
    const double q = (r / 3.0) * ((2.0 * r * r) / 9.0 + s + 2.0 * c_3);

    // Discriminant of the depressed cubic.
    //   Delta = -4 p^3 - 27 q^2.
    // Delta > 0 -> three real roots; Delta = 0 -> repeated; Delta < 0 -> one.
    const double Delta = -4.0 * p * p * p - 27.0 * q * q;

    if (Delta >= 0.0) {
        // p == q == 0 (numerically): triple root at t = 0 of the depressed
        // cubic, i.e. phi = -r/3 with multiplicity three. The kernel has a
        // third-order critical point (ell' = ell'' = 0 simultaneously);
        // there is no isolated mode in the Laplace sense and the caller
        // must fall back. We still return the triple root location for
        // diagnostics; n_modes = 0 signals "no Laplace mode available."
        const double pq_eps = 1e-14;
        if (std::abs(p) < pq_eps && std::abs(q) < pq_eps) {
            const double phi_triple = -r / 3.0;
            for (int k = 0; k < 3; ++k) {
                out.roots[k].phi = phi_triple;
                fill_root_metadata(out.roots[k], A, B, s_jj, alpha);
            }
            out.n_real_roots = 3;
            // is_mode requires ell_pp < 0 strictly; at the triple root
            // ell_pp = 0, so n_modes = 0 and indices stay -1.
            out.n_modes = 0;
            return out;
        }
        // Three real roots (or repeated at Delta == 0). Trigonometric form.
        // p must be <= 0 here; Delta >= 0 with p > 0 is mathematically
        // impossible (-4 p^3 < 0 and -27 q^2 <= 0). Guard for safety.
        if (!(p < 0.0)) {
            // Treat as single real root via cbrt fallback.
            out.n_real_roots = 1;
            out.roots[0].phi = -std::cbrt(q) - r / 3.0;
            fill_root_metadata(out.roots[0], A, B, s_jj, alpha);
            out.n_modes = out.roots[0].is_mode ? 1 : 0;
            out.global_mode_index = out.roots[0].is_mode ? 0 : -1;
            out.status = 3;
            return out;
        }
        const double sqrt_neg_p_third = std::sqrt(-p / 3.0);
        // arg = (3 q / (2 p)) * sqrt(-3 / p). Clamp to [-1, 1] for stability
        // at Delta ~ 0 (where the cubic has a repeated root and arg = +/-1).
        const double raw_arg = (3.0 * q / (2.0 * p)) * std::sqrt(-3.0 / p);
        const double clamp_eps = 1e-9;
        if (std::abs(raw_arg) > 1.0 + clamp_eps) {
            // Beyond the clamp tolerance -> numerical inconsistency
            // between Delta sign and arg magnitude. Fall back.
            out.status = 3;
        }
        const double arg = std::clamp(raw_arg, -1.0, 1.0);
        const double theta_over_3 = std::acos(arg) / 3.0;
        constexpr double two_pi_over_3 = 2.0943951023931953;  // 2 pi / 3
        for (int k = 0; k < 3; ++k) {
            const double t_k = 2.0 * sqrt_neg_p_third
                             * std::cos(theta_over_3 - k * two_pi_over_3);
            out.roots[k].phi = t_k - r / 3.0;
            fill_root_metadata(out.roots[k], A, B, s_jj, alpha);
        }
        out.n_real_roots = 3;
        classify_modes_three(out);
    } else {
        // Delta < 0: one real root. Hyperbolic form, branched by sign of p.
        double t_0;
        if (p < 0.0) {
            const double sqrt_neg_p_third = std::sqrt(-p / 3.0);
            // arg for arccosh: 3 |q| / (-2 p) * sqrt(-3/p) >= 1 when Delta<0.
            const double inner = (3.0 * std::abs(q) / (-2.0 * p))
                                 * std::sqrt(-3.0 / p);
            const double inner_safe = (inner < 1.0) ? 1.0 : inner;
            if (inner < 1.0 - 1e-9) {
                // Numerical inconsistency (Delta < 0 should imply inner >= 1).
                out.status = 3;
            }
            const double sgnq = (q >= 0.0) ? 1.0 : -1.0;
            t_0 = -sgnq * 2.0 * sqrt_neg_p_third
                       * std::cosh(std::acosh(inner_safe) / 3.0);
        } else if (p > 0.0) {
            const double sqrt_p_third = std::sqrt(p / 3.0);
            const double inner = (3.0 * q / (2.0 * p)) * std::sqrt(3.0 / p);
            t_0 = -2.0 * sqrt_p_third * std::sinh(std::asinh(inner) / 3.0);
        } else {
            // p == 0: t^3 + q = 0  =>  t = -cbrt(q).
            t_0 = -std::cbrt(q);
        }
        out.roots[0].phi = t_0 - r / 3.0;
        fill_root_metadata(out.roots[0], A, B, s_jj, alpha);
        out.n_real_roots = 1;
        out.n_modes = out.roots[0].is_mode ? 1 : 0;
        out.global_mode_index = out.roots[0].is_mode ? 0 : -1;
    }

    // Self-check: polynomial residual at each returned root must be small.
    // Scale tolerance by the magnitude of the coefficients.
    const double scale = 1.0 + std::abs(b_coef) + std::abs(c_coef)
                              + std::abs(d_coef);
    const double tol = 1e-7 * scale * (1.0 + std::abs(out.roots[0].phi));
    for (int k = 0; k < out.n_real_roots; ++k) {
        const double rk = cubic_residual(out.roots[k].phi,
                                         b_coef, c_coef, d_coef);
        if (!(std::abs(rk) < tol)) {
            out.status = 3;
        }
    }

    return out;
}

}  // namespace savage_dickey
