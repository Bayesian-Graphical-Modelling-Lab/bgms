#include "sinh_midpoint.h"

#include <cmath>

namespace savage_dickey {

namespace {

// log of the original integrand factor: -A phi^2 + B phi + (alpha-1) log(s + phi^2).
inline double log_kernel(double phi, double A, double B,
                         double s, double alpha) {
    const double a1 = alpha - 1.0;
    double v = -A * phi * phi + B * phi;
    if (a1 != 0.0) v += a1 * std::log(s + phi * phi);
    return v;
}

// Number of "standard deviations" (in the sinh coordinate, where the
// Gaussian envelope -A s sinh^2(t) + B sqrt(s) sinh(t) is exactly
// Gaussian with variance 1/(2 A s)) to include in the integration
// interval. n=10 gives Gaussian-envelope tail mass of roughly exp(-50)
// ≈ 2e-22. The arcsinh mapping back to t naturally compresses for
// small A * s_jj (envelope wide in sinh, bounded in t).
//
// Note: an alternative is to use t-sigma via the integrand's Laplace
// curvature kappa_t = 2 A s + B^2/(2 A). That captures the local
// concavity of the integrand at the mode but underestimates the
// integrand's actual width when the integrand isn't well-approximated
// by a Gaussian in t (e.g., for moderate alpha > 1 with sinh(t*)
// significantly nonzero). The sinh-sigma criterion is conservatively
// wider and tracks the true Gaussian-in-sinh decay of the envelope.
constexpr double kNSigmaSinh = 10.0;

}  // namespace

LSDSinhResult density_at_l_ji_sinh(double x_eval, double A, double B,
                                    double s_jj, double alpha,
                                    int num_nodes) {
    LSDSinhResult out;
    out.log_density = arma::datum::nan;
    out.log_Z       = arma::datum::nan;
    out.status      = 0;
    if (!(A > 0.0))    { out.status = 1; return out; }
    if (!(s_jj > 0.0)) { out.status = 2; return out; }

    const double sqrt_s = std::sqrt(s_jj);
    const double As     = A * s_jj;

    // Gaussian-envelope mode in t: sinh(t*) = B / (2 A sqrt(s)).
    // Variance in the sinh coordinate is 1/(2 A s).
    const double sinh_t_star = B / (2.0 * A * sqrt_s);
    const double sigma_sinh  = 1.0 / std::sqrt(2.0 * As);

    // Truncation interval: +- kNSigmaSinh standard deviations of the
    // envelope in the sinh coordinate, mapped back to t via arcsinh.
    const double sinh_left  = sinh_t_star - kNSigmaSinh * sigma_sinh;
    const double sinh_right = sinh_t_star + kNSigmaSinh * sigma_sinh;
    const double T_left     = std::asinh(sinh_left);
    const double T_right    = std::asinh(sinh_right);
    const double h          = (T_right - T_left) / num_nodes;

    // Midpoint rule: evaluate at t_k = T_left + (k + 1/2) * h.
    // log g(t_k) = (2 alpha - 1) log cosh(t_k)
    //            - A s sinh^2(t_k)
    //            + B sqrt(s) sinh(t_k).
    const double two_alpha_minus_1 = 2.0 * alpha - 1.0;
    arma::vec logterms(num_nodes);
    for (int k = 0; k < num_nodes; ++k) {
        const double t_k     = T_left + (static_cast<double>(k) + 0.5) * h;
        const double sinh_tk = std::sinh(t_k);
        const double cosh_tk = std::cosh(t_k);
        logterms(k) = two_alpha_minus_1 * std::log(cosh_tk)
                    - As * sinh_tk * sinh_tk
                    + B * sqrt_s * sinh_tk;
    }
    const double max_lt = logterms.max();
    const double lse =
        max_lt + std::log(arma::accu(arma::exp(logterms - max_lt)));

    // log Z = log(s^(alpha - 1/2)) + log(h) + log( sum_k g(t_k) ).
    const double log_prefactor = (alpha - 0.5) * std::log(s_jj);
    out.log_Z       = log_prefactor + std::log(h) + lse;
    out.log_density = log_kernel(x_eval, A, B, s_jj, alpha) - out.log_Z;
    out.status      = 0;
    return out;
}

}  // namespace savage_dickey
