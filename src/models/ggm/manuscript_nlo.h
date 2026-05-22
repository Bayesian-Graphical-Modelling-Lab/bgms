// Closed-form NLO under the manuscript App C decomposition.
//
// eq:NLO-decomp, eq:Be / eq:Cel / eq:Di / eq:Elm. A_l is omitted because the
// exact log_C0 (Bartlett base, lgamma terms) already handles the diagonal
// class-(i) integrals; including A_l would double-count and shift log Z(empty)
// away from 0.
//
// Caveat: accurate in the sparse regime (n_edges <= ~q*q/4). Over-corrects at
// dense / high-max-degree centres (asymptotic-series breakdown). At p=50
// operational density (~5%) the sparse regime applies.
//
// Companion-AI delivery: ~/SV/Z/R/src/manuscript_NLO.h (2026-05-21).

#pragma once

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

namespace ggm_nlo {

// Manuscript NLO at alpha = 1 (no Na/Nb/Nc terms; those vanish at alpha=1).
// Returns log_C0 + log_LO + delta_NLO_manuscript.
inline double log_Z_manuscript_NLO_alpha1(
    const arma::imat& G,
    double beta, double sigma, double delta
) {
    int q = G.n_rows;
    std::vector<int> nu(q, 0);
    int E_count = 0;
    for (int l = 0; l < q; ++l)
        for (int m = l + 1; m < q; ++m)
            if (G(l, m) == 1) { ++nu[l]; ++E_count; }

    const double sigma2 = sigma * sigma;
    const double sigma4 = sigma2 * sigma2;
    const double lambda = beta;
    const double lambda2 = lambda * lambda;

    // log C0 (Bartlett base, exact; alpha = 1)
    double log_C0 = 0.5 * static_cast<double>(E_count) * std::log(M_PI)
                    - 0.5 * static_cast<double>(E_count) * std::log(2.0 * M_PI * sigma2)
                    - static_cast<double>(E_count) * std::log(beta)
                    - static_cast<double>(q) * std::lgamma(1.0)      // alpha = 1 -> 0
                    - static_cast<double>(q) * delta * std::log(beta);
    for (int l = 0; l < q; ++l)
        log_C0 += std::lgamma((static_cast<double>(nu[l]) + 2.0 * (1.0 + delta)) / 2.0);
    if (E_count == 0) return log_C0;

    // M_tilde_l = nu_l + 1 + 2 delta;  H_l^off = 2 lambda + M_tilde_l/(2 sigma^2 lambda)
    std::vector<double> M_tilde(q), H_off(q);
    for (int l = 0; l < q; ++l) {
        M_tilde[l] = static_cast<double>(nu[l]) + 1.0 + 2.0 * delta;
        H_off[l]   = 2.0 * lambda + M_tilde[l] / (2.0 * sigma2 * lambda);
    }

    // log_LO = sum_e 0.5 log(2 lambda / H_e), with H_e = H_{i_e}^off
    double log_LO = 0.0;
    for (int i = 0; i < q; ++i)
        for (int j = i + 1; j < q; ++j)
            if (G(i, j) == 1) {
                double H_e = H_off[i];
                if (H_e <= 0.0) return R_NegInf;
                log_LO += 0.5 * std::log(2.0 * lambda / H_e);
            }

    // delta_NLO = sum_e B_e + sum_(e, l in C_e^off) C_{e,l}
    //           + sum_{i: nu_i >= 1} D_i + sum_{(l,m) not in E, l<m} E_{l,m}
    double dNLO = 0.0;

    // B_e
    for (int i = 0; i < q; ++i)
        for (int j = i + 1; j < q; ++j)
            if (G(i, j) == 1) {
                double H_e = H_off[i];
                dNLO += -1.0 / (4.0 * lambda * sigma2 * H_e);
            }

    // C_{e,l}: edge e = (i, j), common predecessor l < i with G(l, i) = G(l, j) = 1
    for (int i = 1; i < q; ++i)
        for (int j = i + 1; j < q; ++j) {
            if (G(i, j) != 1) continue;
            double H_e = H_off[i];
            for (int l = 0; l < i; ++l)
                if (G(l, i) == 1 && G(l, j) == 1) {
                    double Hl = H_off[l];
                    double Hl2 = Hl * Hl;
                    dNLO += -1.0 / (2.0 * sigma2 * Hl2)
                            + M_tilde[i] / (4.0 * lambda * sigma4 * H_e * Hl2);
                }
        }

    // D_i (only vertices with nu_i >= 1)
    for (int i = 0; i < q; ++i)
        if (nu[i] >= 1) {
            double Hi = H_off[i];
            dNLO += static_cast<double>(nu[i]) * (static_cast<double>(nu[i]) + 2.0)
                    * M_tilde[i]
                    / (16.0 * lambda2 * sigma4 * Hi * Hi);
        }

    // E_{l,m}: non-edges (l, m), l < m, common predecessors k < l with G(k,l)=G(k,m)=1
    for (int l = 0; l < q - 1; ++l)
        for (int m = l + 1; m < q; ++m) {
            if (G(l, m) == 1) continue;
            if (l < 1) continue;  // need k < l
            double sum_inv_H2 = 0.0;
            for (int k = 0; k < l; ++k)
                if (G(k, l) == 1 && G(k, m) == 1) {
                    double Hk = H_off[k];
                    sum_inv_H2 += 1.0 / (Hk * Hk);
                }
            if (sum_inv_H2 > 0.0) {
                dNLO += -2.0 * lambda2 / M_tilde[l] * sum_inv_H2;
            }
        }

    return log_C0 + log_LO + dNLO;
}


// Toggle-endpoint reordering ("DEGORD"): same convention as
// log_Z_NLO_gamma_degord — relabel (i, j) to (0, 1) and permute remaining
// vertices in their original order before applying the closed form. Required
// for the chain MH ratio because the closed form is not permutation-invariant.
inline double log_Z_manuscript_NLO_alpha1_degord(
    const arma::imat& G, int i, int j,
    double beta, double sigma, double delta
) {
    int q = G.n_rows;
    std::vector<int> perm;
    perm.reserve(q);
    perm.push_back(i);
    perm.push_back(j);
    for (int v = 0; v < q; ++v) if (v != i && v != j) perm.push_back(v);
    arma::imat G_perm(q, q, arma::fill::zeros);
    for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b)
            G_perm(a, b) = G(perm[a], perm[b]);
    return log_Z_manuscript_NLO_alpha1(G_perm, beta, sigma, delta);
}


}  // namespace ggm_nlo
