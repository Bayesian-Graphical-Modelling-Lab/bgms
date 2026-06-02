#include "cauchy_omega.h"

#include <cmath>
#include <limits>

namespace savage_dickey {

namespace {

constexpr double kSliceStepWidth = 1.0;
constexpr int    kSliceMaxExpand = 50;
constexpr int    kSliceMaxShrink = 100;

// Log-target on u = log(omega).  Up to a constant in u.  See the file
// header of cauchy_omega.h for the derivation.
//
//   f(u) = -K^2 / (2 sigma^2 omega) - 1/(2 omega) - 0.5 u
//          + 0.5 log A_K(omega) - B_K^2 / (4 A_K(omega))
//
// where A_K(omega) = 0.5 / (sigma^2 omega) + A_diag_K. The IG(1/2, 1/2)
// prior on omega contributes -1.5 log omega - 1/(2 omega); combined with
// the +u Jacobian for the u = log(omega) reparameterisation we get the
// -0.5 u term.
inline double log_target_u(double u,
                           double sigma2,
                           double A_diag_K,
                           double half_K2,
                           double half_B2) {
    const double inv_w = std::exp(-u);
    const double A_K   = 0.5 * inv_w / sigma2 + A_diag_K;
    if (!(A_K > 0.0) || !std::isfinite(A_K)) {
        return -std::numeric_limits<double>::infinity();
    }
    return -half_K2 * inv_w / sigma2
         - 0.5 * inv_w
         - 0.5 * u
         + 0.5 * std::log(A_K)
         - half_B2 / A_K;
}

}  // namespace

double slice_sample_cauchy_omega_active(
    double K,
    double sigma2,
    double A_diag_K,
    double B_K,
    double omega_curr,
    SafeRNG& rng
) {
    const double half_K2 = 0.5 * K * K;
    const double half_B2 = 0.25 * B_K * B_K;

    double u_curr = std::log(omega_curr);
    if (!std::isfinite(u_curr)) u_curr = 0.0;

    const double log_y =
        log_target_u(u_curr, sigma2, A_diag_K, half_K2, half_B2)
        - rexp(rng, 1.0);

    double L = u_curr - kSliceStepWidth * runif(rng);
    double R = L + kSliceStepWidth;
    for (int k = 0; k < kSliceMaxExpand; ++k) {
        if (log_target_u(L, sigma2, A_diag_K, half_K2, half_B2) <= log_y) break;
        L -= kSliceStepWidth;
    }
    for (int k = 0; k < kSliceMaxExpand; ++k) {
        if (log_target_u(R, sigma2, A_diag_K, half_K2, half_B2) <= log_y) break;
        R += kSliceStepWidth;
    }

    double u_new = u_curr;
    bool   accepted = false;
    for (int k = 0; k < kSliceMaxShrink; ++k) {
        u_new = L + (R - L) * runif(rng);
        if (log_target_u(u_new, sigma2, A_diag_K, half_K2, half_B2) > log_y) {
            accepted = true;
            break;
        }
        if (u_new < u_curr) L = u_new; else R = u_new;
    }
    if (!accepted) return omega_curr;

    const double w_new = std::exp(u_new);
    return (std::isfinite(w_new) && w_new > 0.0) ? w_new : omega_curr;
}

double sample_cauchy_omega_prior(SafeRNG& rng) {
    const double v = rgamma(rng, 0.5, 0.5);
    return 1.0 / std::max(v, 1e-300);
}

}  // namespace savage_dickey
