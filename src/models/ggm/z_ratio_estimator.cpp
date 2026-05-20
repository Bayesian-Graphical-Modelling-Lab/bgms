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


// log|expm1(x)| with explicit sign(expm1(x)).
//   x > 0: expm1(x) > 0, log|.| = x + log1p(-exp(-x))
//   x < 0: expm1(x) < 0, log|.| = log1p(-exp(x))
//   x == 0: expm1(0) = 0; caller treats as sign=0, contribution=0.
// Both numerically stable near x = 0 (log1p(-exp(-|x|)) ~ log|x|).
static inline std::pair<double, int> log_abs_expm1_signed(double x) {
    if (!std::isfinite(x)) {
        return {std::numeric_limits<double>::quiet_NaN(), 0};
    }
    if (x == 0.0) {
        return {-std::numeric_limits<double>::infinity(), 0};
    }
    if (x > 0.0) {
        return {x + std::log1p(-std::exp(-x)), +1};
    }
    return {std::log1p(-std::exp(x)), -1};
}


// Build (log|V|, sign(V)) from a sequence of log_Zhat_n values (n = 1..K)
// at fixed (log_c, rho). Computes log|S|, sign(S) over the K+1 truncated
// series terms via signed log-sum-exp, then returns (-log_c + log|S|,
// sign(S)). K = 0 short-circuits to {-log_c, +1}.
//
// Returns {NaN, 0} on any non-finite intermediate or S = 0 collapse.
static std::pair<double, int> V_log_from_log_Zhats(
    const std::vector<double>& log_Zhats,
    double log_c,
    double rho
) {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double neg_inf = -std::numeric_limits<double>::infinity();

    if (!std::isfinite(log_c) || !(rho > 0.0 && rho < 1.0)) {
        return {NaN, 0};
    }

    const int K_depth = static_cast<int>(log_Zhats.size());

    // Accumulate (log|term_n|, sign(term_n)) for n = 0..K_depth, then resolve
    // via signed log-sum-exp at the end.
    std::vector<double> log_abs;
    std::vector<int>    sgn;
    log_abs.reserve(static_cast<size_t>(K_depth) + 1u);
    sgn.reserve(static_cast<size_t>(K_depth) + 1u);

    // Term n = 0: +1.
    log_abs.push_back(0.0);
    sgn.push_back(+1);

    if (K_depth > 0) {
        const double log_rho = std::log(rho);
        double log_abs_prod = 0.0;
        int    sign_prod = +1;
        for (int n = 1; n <= K_depth; ++n) {
            double log_Zhat_n = log_Zhats[static_cast<size_t>(n - 1)];
            if (!std::isfinite(log_Zhat_n)) return {NaN, 0};

            double x = log_Zhat_n - log_c;
            auto em1 = log_abs_expm1_signed(x);
            double log_abs_em1 = em1.first;
            int    sgn_em1     = em1.second;
            if (!std::isfinite(log_abs_em1) && sgn_em1 == 0 && x == 0.0) {
                // expm1(x) = 0 exactly → all further (and this) term_n vanish.
                // Subsequent terms inherit the same zero factor, so we are
                // done extending the truncated series.
                break;
            }
            if (sgn_em1 == 0) {
                // Non-finite x; bail out.
                return {NaN, 0};
            }
            log_abs_prod += log_abs_em1;
            sign_prod    *= sgn_em1;

            double log_abs_term = -static_cast<double>(n) * log_rho + log_abs_prod;
            int    sign_term    = ((n & 1) ? -1 : +1) * sign_prod;

            if (!std::isfinite(log_abs_term)) {
                // Overflow in magnitude (e.g. tail of the series ran away).
                // log_Zhat overflow already screened above; this guards the
                // log1p(-exp) edge cases.
                if (log_abs_term == neg_inf) continue;
                return {NaN, 0};
            }

            log_abs.push_back(log_abs_term);
            sgn.push_back(sign_term);
        }
    }

    // Signed log-sum-exp.
    double M = neg_inf;
    for (double la : log_abs) {
        if (la > M) M = la;
    }
    if (M == neg_inf) {
        // All terms exactly zero; S = 0, sign undefined.
        return {NaN, 0};
    }
    double s = 0.0;
    for (size_t k = 0; k < log_abs.size(); ++k) {
        if (log_abs[k] == neg_inf) continue;
        s += static_cast<double>(sgn[k]) * std::exp(log_abs[k] - M);
    }
    if (s == 0.0 || !std::isfinite(s)) {
        return {NaN, 0};
    }
    double log_abs_S = M + std::log(std::abs(s));
    int    sign_S    = (s > 0.0) ? +1 : -1;

    double log_abs_V = -log_c + log_abs_S;
    if (!std::isfinite(log_abs_V)) {
        return {NaN, 0};
    }
    return {log_abs_V, sign_S};
}


std::pair<double, int> V_log_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const PiAux& pi_aux,
    const ChainAux& chain_aux,
    double log_c,
    double rho
) {
    std::vector<double> log_Zhats;
    log_Zhats.reserve(static_cast<size_t>(K_depth));
    for (int n = 0; n < K_depth; ++n) {
        double log_Zhat_n = log_Zhat_pi_from_pool(
            pools_t[n], pi_aux, chain_aux);
        if (!std::isfinite(log_Zhat_n)) {
            return {std::numeric_limits<double>::quiet_NaN(), 0};
        }
        log_Zhats.push_back(log_Zhat_n);
    }
    return V_log_from_log_Zhats(log_Zhats, log_c, rho);
}


std::pair<double, int> V_log_at_Gamma_pi_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const arma::imat& G_pi,
    const ChainAux& chain_aux,
    double log_c,
    double rho
) {
    PiAux pi_aux = make_pi_aux(G_pi, chain_aux);
    return V_log_at_Gamma_pi_degord(
        K_depth, pools_t, pi_aux, chain_aux, log_c, rho);
}


LogSignedVPair V_log_pair_at_Gamma_curr_star_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const PiAux& a_curr,
    const PiAux& a_star,
    const ChainAux& chain_aux,
    double log_c_curr,
    double log_c_star,
    double rho
) {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    LogSignedVPair out;
    out.curr = {NaN, 0};
    out.star = {NaN, 0};

    // Loop once over the K_depth pools, building cache_curr from a_curr and
    // reusing it to evaluate log_Zhat_n at a_star without rebuilding Phi.
    std::vector<double> log_Zhats_curr, log_Zhats_star;
    log_Zhats_curr.reserve(static_cast<size_t>(K_depth));
    log_Zhats_star.reserve(static_cast<size_t>(K_depth));
    for (int n = 0; n < K_depth; ++n) {
        PoolCache cache_curr;
        double log_Zhat_n_curr = log_Zhat_pi_from_pool_cache(
            pools_t[n], a_curr, chain_aux, cache_curr);
        if (!std::isfinite(log_Zhat_n_curr)) return out;
        double log_Zhat_n_star = log_Zhat_star_from_cache(
            pools_t[n], a_star, chain_aux, cache_curr);
        if (!std::isfinite(log_Zhat_n_star)) return out;
        log_Zhats_curr.push_back(log_Zhat_n_curr);
        log_Zhats_star.push_back(log_Zhat_n_star);
    }

    out.curr = V_log_from_log_Zhats(log_Zhats_curr, log_c_curr, rho);
    out.star = V_log_from_log_Zhats(log_Zhats_star, log_c_star, rho);
    return out;
}


LogSignedVPair V_log_pair_at_Gamma_curr_star_degord(
    int K_depth,
    const std::vector<arma::mat>& pools_t,
    const arma::imat& G_pi_curr,
    const arma::imat& G_pi_star,
    const ChainAux& chain_aux,
    double log_c_curr,
    double log_c_star,
    double rho
) {
    PiAux a_curr = make_pi_aux(G_pi_curr, chain_aux);
    PiAux a_star = make_pi_aux(G_pi_star, chain_aux);
    return V_log_pair_at_Gamma_curr_star_degord(
        K_depth, pools_t, a_curr, a_star, chain_aux,
        log_c_curr, log_c_star, rho);
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
