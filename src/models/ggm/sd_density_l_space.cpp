#include "sd_density_l_space.h"

#include <cmath>
#include <limits>

namespace ggm_sd {

namespace {

// log-kernel derivatives at x. With g(x) = log(s + x²):
//   f'(x)   = -2 A x + B + (alpha - 1) g'(x)
//   f''(x)  = -2 A      + (alpha - 1) g''(x)
//   f'''(x) =             (alpha - 1) g'''(x)
//   f''''(x) =            (alpha - 1) g''''(x)
//
// g'(x)    = 2 x / (s + x²)
// g''(x)   = 2 (s - x²) / (s + x²)²
// g'''(x)  = 4 x (x² - 3 s) / (s + x²)³
// g''''(x) = -12 (x⁴ - 6 s x² + s²) / (s + x²)⁴
struct Derivs {
    double f1, f2, f3, f4;
};

inline Derivs log_kernel_derivs(double x, double A, double B,
                                double s, double alpha) {
    const double a1 = alpha - 1.0;
    Derivs d;
    if (a1 == 0.0) {
        d.f1 = -2.0 * A * x + B;
        d.f2 = -2.0 * A;
        d.f3 = 0.0;
        d.f4 = 0.0;
        return d;
    }
    const double q     = s + x * x;
    const double q2    = q * q;
    const double q3    = q2 * q;
    const double q4    = q2 * q2;
    const double smx2  = s - x * x;
    const double g1    = 2.0 * x / q;
    const double g2    = 2.0 * smx2 / q2;
    const double g3    = 4.0 * x * (x * x - 3.0 * s) / q3;
    const double xsq   = x * x;
    const double g4    = -12.0 * (xsq * xsq - 6.0 * s * xsq + s * s) / q4;
    d.f1 = -2.0 * A * x + B + a1 * g1;
    d.f2 = -2.0 * A          + a1 * g2;
    d.f3 =                     a1 * g3;
    d.f4 =                     a1 * g4;
    return d;
}

inline double log_kernel(double x, double A, double B,
                         double s, double alpha) {
    const double a1 = alpha - 1.0;
    double val = -A * x * x + B * x;
    if (a1 != 0.0) val += a1 * std::log(s + x * x);
    return val;
}

}  // namespace

LSDResult density_at_l_ji_one(double x_eval,
                              double A, double B,
                              double s_jj,
                              double alpha,
                              bool   nlo,
                              int    newton_max_iter,
                              double newton_tol) {
    LSDResult out;
    out.log_density = arma::datum::nan;
    out.x_mode      = arma::datum::nan;
    out.curvature   = arma::datum::nan;
    out.status      = 0;

    // Newton iteration. Start at the Gaussian mode (closed-form when alpha=1).
    double x = (A > 0.0) ? B / (2.0 * A) : 0.0;
    bool   converged = (alpha == 1.0);  // exact at alpha=1, no iteration needed
    if (!converged) {
        for (int it = 0; it < newton_max_iter; ++it) {
            const Derivs d = log_kernel_derivs(x, A, B, s_jj, alpha);
            if (!(std::abs(d.f2) > 1e-14)) break;
            const double step = d.f1 / d.f2;
            x -= step;
            if (std::abs(step) < newton_tol) { converged = true; break; }
        }
    }

    const Derivs d_mode = log_kernel_derivs(x, A, B, s_jj, alpha);
    const double curvature = -d_mode.f2;
    out.x_mode    = x;
    out.curvature = curvature;
    if (!std::isfinite(curvature) || !(curvature > 0.0)) {
        out.status = 2;
        return out;
    }

    // Laplace log-Z ≈ f(x_mode) + ½ log(2π) − ½ log(curvature).
    const double f_mode   = log_kernel(x, A, B, s_jj, alpha);
    const double log_Z_LP = f_mode + 0.5 * std::log(2.0 * arma::datum::pi)
                                   - 0.5 * std::log(curvature);

    // Tierney-Kadane 1/n correction for log ∫ exp(f(x)) dx. Derivation: let
    // κ = −f''(x*), substitute u = (x − x*) √κ, expand exp(f) around the
    // mode, take E[·] under N(0, 1). To leading 1/n:
    //   log_Z ≈ log_Z_LP + f''''(x*)/(8 κ²) + (5/24) (f'''(x*))² / κ³
    // Both terms positive at a max (κ > 0, signs absorbed by squares).
    // At α=1 the (α−1) factor in f''' and f'''' is zero ⇒ NLO ≡ 0 ⇒
    // Laplace is exact, matching the closed-form Gaussian.
    double nlo_correction = 0.0;
    if (nlo && (alpha != 1.0)) {
        const double L3      = d_mode.f3;
        const double L4      = d_mode.f4;
        const double abs_Lpp = curvature;          // = -f''(x_mode), > 0
        nlo_correction = L4 / (8.0 * abs_Lpp * abs_Lpp)
                       + 5.0 * L3 * L3 / (24.0 * abs_Lpp * abs_Lpp * abs_Lpp);
    }

    // log π(x_eval | rest) = f(x_eval) − log_Z.
    const double f_eval = log_kernel(x_eval, A, B, s_jj, alpha);
    out.log_density = f_eval - (log_Z_LP + nlo_correction);
    out.status      = converged ? 0 : 3;
    return out;
}

}  // namespace ggm_sd
