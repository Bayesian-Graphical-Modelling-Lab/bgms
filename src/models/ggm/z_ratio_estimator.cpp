// Russian-Roulette V(Γ, U) estimator of 1/Z(Γ). Built on the Phase 2 DEGORD
// Bartlett-Cholesky sampler. See z_ratio_estimator.h for the construction.
//
// Port of ~/SV/Z/R/src/branchB_chain_route3a_degord.cpp:192-228 (V function)
// and :215-228 (pool draw), with R::rnorm / R::rgeom replaced by SafeRNG
// for chain-seed portability.

#include "models/ggm/z_ratio_estimator.h"
#include <cmath>
#include <limits>

#include "rng/rng_utils.h"

namespace degord {


double V_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const PiAux& pi_aux,
    const ChainAux& chain_aux,
    double c_val,
    double rho
) {
    if (K_depth == 0) return 1.0 / c_val;
    double acc = 1.0 / c_val;
    double running_prod = 1.0;
    for (int n = 1; n <= K_depth; ++n) {
        double log_Zhat_n = log_Zhat_pi_from_pool(
            pools_t[n - 1], pi_aux, chain_aux);
        if (!std::isfinite(log_Zhat_n))
            return std::numeric_limits<double>::quiet_NaN();
        double Zhat_n = std::exp(log_Zhat_n);
        running_prod *= (Zhat_n - c_val) / c_val;
        double sgn = (n % 2 == 0) ? 1.0 : -1.0;
        acc += sgn * running_prod / (c_val * std::pow(rho, static_cast<double>(n)));
    }
    return acc;
}


double V_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const arma::imat& G_pi,
    const ChainAux& chain_aux,
    double c_val,
    double rho
) {
    PiAux pi_aux = make_pi_aux(G_pi, chain_aux);
    return V_at_Gamma_pi_degord(K_depth, pools_t, pi_aux, chain_aux, c_val, rho);
}


void draw_U_degord_rr(
    SafeRNG& rng,
    int& K_depth,
    std::vector<arma::mat>& pools_t,
    int M_inner,
    int q,
    double rho
) {
    // K_depth ~ Geom(1 - rho). boost::random doesn't ship a geometric directly
    // here; we draw via inverse-CDF on a uniform: K = floor(log(U) / log(rho)).
    // This matches R::rgeom(1 - rho) in distribution (number of failures
    // before the first success when success prob = 1 - rho).
    double u = runif(rng);
    if (u <= 0.0) u = std::numeric_limits<double>::min();  // guard log(0)
    K_depth = static_cast<int>(std::floor(std::log(u) / std::log(rho)));

    pools_t.clear();
    pools_t.reserve(static_cast<size_t>(K_depth));
    for (int n = 0; n < K_depth; ++n) {
        pools_t.push_back(draw_bartlett_pool(rng, q, M_inner));
    }
}


}  // namespace degord
