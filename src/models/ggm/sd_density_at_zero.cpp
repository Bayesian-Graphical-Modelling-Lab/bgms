#include "sd_density_at_zero.h"

#include <Rcpp.h>
#include <cmath>
#include <limits>

namespace ggm_sd {

SDResult density_at_zero_one(
    const arma::mat& K, int i, int j,
    const arma::mat& S, int n_obs,
    double delta, double sigma,
    bool nlo,
    bool apply_pd_truncation,
    int newton_max_iter,
    double newton_tol)
{
    SDResult out;
    out.status      = 0;
    out.log_density = arma::datum::nan;
    out.x_mode      = arma::datum::nan;
    out.curvature   = arma::datum::nan;
    out.x_minus     = arma::datum::nan;
    out.x_plus      = arma::datum::nan;
    out.log_Z_trunc = arma::datum::nan;

    const int q = K.n_rows;

    arma::mat K0 = K;
    K0(i, j) = 0.0;
    K0(j, i) = 0.0;

    arma::mat L_chol;
    if (!arma::chol(L_chol, K0, "lower")) {
        out.log_density = -arma::datum::inf;
        out.status      = 1;
        return out;
    }

    arma::vec e_i(q, arma::fill::zeros); e_i(i) = 1.0;
    arma::vec e_j(q, arma::fill::zeros); e_j(j) = 1.0;
    arma::vec tmp_i        = arma::solve(arma::trimatl(L_chol), e_i);
    arma::vec sigma_col_i  = arma::solve(arma::trimatu(L_chol.t()), tmp_i);
    arma::vec tmp_j        = arma::solve(arma::trimatl(L_chol), e_j);
    arma::vec sigma_col_j  = arma::solve(arma::trimatu(L_chol.t()), tmp_j);

    const double sig_ii = sigma_col_i(i);
    const double sig_jj = sigma_col_j(j);
    const double sig_ij = sigma_col_i(j);

    const double c1 = sig_ij;
    const double c2 = sig_ii * sig_jj - sig_ij * sig_ij;

    // PD-feasible interval for K_ij given K_{-ij}. Roots of D(x) = 0:
    //   x_pm = (c1 +/- sqrt(c1^2 + c2)) / c2,    c1^2 + c2 = sig_ii * sig_jj > 0.
    const double disc      = c1 * c1 + c2;
    const double sqrt_disc = std::sqrt(disc);
    out.x_minus = (c1 - sqrt_disc) / c2;
    out.x_plus  = (c1 + sqrt_disc) / c2;

    const double log_det_factor = delta + n_obs / 2.0;
    const double S_ij           = S(i, j);
    const double inv_sigma2     = 1.0 / (sigma * sigma);

    double x = 0.0;
    bool converged = false;
    for (int iter = 0; iter < newton_max_iter; ++iter) {
        const double denom = 1.0 + 2.0 * c1 * x - c2 * x * x;
        if (denom <= 0.0) break;
        const double gp  = 2.0 * c1 - 2.0 * c2 * x;
        const double Lp  = log_det_factor * gp / denom - S_ij - x * inv_sigma2;
        const double Lpp = log_det_factor
                           * (-2.0 * c2 / denom - (gp * gp) / (denom * denom))
                           - inv_sigma2;
        if (std::abs(Lpp) < 1e-14) break;
        const double step = Lp / Lpp;
        x -= step;
        if (std::abs(step) < newton_tol) {
            converged = true;
            break;
        }
    }

    const double denom_m = 1.0 + 2.0 * c1 * x - c2 * x * x;
    if (denom_m <= 0.0) {
        out.x_mode = x;
        out.status = 2;
        return out;
    }

    const double gp_m = 2.0 * c1 - 2.0 * c2 * x;
    const double u    = gp_m / denom_m;
    const double v    = -2.0 * c2 / denom_m;
    const double Lpp  = log_det_factor * (v - u * u) - inv_sigma2;
    const double curvature = -Lpp;

    if (!std::isfinite(curvature) || curvature <= 0.0) {
        out.x_mode    = x;
        out.curvature = curvature;
        out.status    = 2;
        return out;
    }

    const double L_at_mode = log_det_factor * std::log(denom_m) - S_ij * x
                             - 0.5 * x * x * inv_sigma2;
    const double log_norm  = L_at_mode + 0.5 * std::log(2.0 * M_PI / curvature);

    double nlo_correction = 0.0;
    if (nlo) {
        const double phi3 = -3.0 * u * v + 2.0 * u * u * u;
        const double phi4 = -3.0 * v * v + 12.0 * u * u * v
                            - 6.0 * std::pow(u, 4);
        const double L3       = log_det_factor * phi3;
        const double L4       = log_det_factor * phi4;
        const double abs_Lpp  = -Lpp;
        nlo_correction = (1.0 / 8.0)  * L4 / (Lpp * Lpp)
                       + (5.0 / 24.0) * (L3 * L3)
                                      / (abs_Lpp * abs_Lpp * abs_Lpp);
    }

    // log Z_trunc = log P(N(x_mode, 1/sqrt(curvature)) in (x_minus, x_plus)).
    // Closed-form via the normal CDF. Use R::pnorm5(q, mu, sigma, lower=1, log_p=0)
    // for numerical stability in the tails; combine via pmax of the difference
    // to guard against negative arguments from cancellation.
    const double sd_prop = 1.0 / std::sqrt(curvature);
    const double F_lo    = R::pnorm5(out.x_minus, x, sd_prop, 1, 0);
    const double F_hi    = R::pnorm5(out.x_plus,  x, sd_prop, 1, 0);
    const double mass    = F_hi - F_lo;
    out.log_Z_trunc = (mass > 0.0 && std::isfinite(mass)) ? std::log(mass) : 0.0;

    double base = log_norm + nlo_correction;
    if (apply_pd_truncation) base += out.log_Z_trunc;
    out.log_density = -base;
    out.x_mode      = x;
    out.curvature   = curvature;
    out.status      = converged ? 0 : 3;
    return out;
}

}  // namespace ggm_sd
