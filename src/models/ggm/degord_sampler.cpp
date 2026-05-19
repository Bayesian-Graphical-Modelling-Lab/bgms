// DEGORD-permuted Bartlett-Cholesky importance sampler for log Zhat(G).
// Direct port of ~/SV/Z/R/src/degord_sampler.h (v4, 2026-05-18) with the
// header-only static-inline bodies moved into one translation unit, and
// a SafeRNG-based Bartlett pool draw added at the end.
//
// The kernel is bit-identical to the z reference up to floating-point
// reordering when given the same noise pool and ChainAux/PiAux state.

#include "models/ggm/degord_sampler.h"
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <algorithm>

namespace degord {

ChainAux make_chain_aux(int q, double alpha, double beta,
                        double sigma, double delta) {
  ChainAux c;
  c.q = q;
  c.alpha = alpha; c.beta = beta; c.sigma = sigma; c.delta = delta;
  c.sigma2 = sigma * sigma;
  c.inv_sigma2 = 1.0 / c.sigma2;
  c.two_beta = 2.0 * beta;
  c.inv_two_beta = 1.0 / c.two_beta;
  c.log_sigma = std::log(sigma);
  c.log_2pi = std::log(2.0 * M_PI);
  c.slab_kernel_const = -0.5 * std::log(2.0 * M_PI * c.sigma2);
  c.sigma_diag = 1.0 / std::sqrt(4.0 * beta);
  c.inv_sigma_diag2 = 4.0 * beta;
  c.half_log_2pi_sigma_diag2 = 0.5 * std::log(2.0 * M_PI / (4.0 * beta));
  c.slab_tilt_mode = 0;

  int nmax = q + 1;
  c.nu_chi_df.assign(nmax, 0.0);
  c.nu_half_k.assign(nmax, 0.0);
  c.nu_km1.assign(nmax, 0.0);
  c.nu_lgamma_half_k.assign(nmax, 0.0);
  c.nu_mu_l.assign(nmax, 0.0);
  c.nu_diag_const.assign(nmax, 0.0);
  c.nu_H_e_saddle.assign(nmax, 0.0);
  c.nu_inv_sqrt_H_e_saddle.assign(nmax, 0.0);
  c.nu_mu_coef_saddle.assign(nmax, 0.0);
  c.nu_sigma_star_sq_saddle.assign(nmax, 0.0);
  c.nu_inv_sigma_star_sq_saddle.assign(nmax, 0.0);
  c.nu_slab_const_saddle.assign(nmax, 0.0);
  c.nu_per_vertex.assign(nmax, 0.0);

  double log_beta = std::log(beta);
  for (int nu = 0; nu < nmax; ++nu) {
    double k    = static_cast<double>(nu) + 2.0 + 2.0 * delta;
    double hk   = 0.5 * k;
    double km1  = k - 1.0;
    double lgk  = std::lgamma(hk);
    double mu_l = std::sqrt(km1 * c.inv_two_beta);
    double mu_l2 = km1 * c.inv_two_beta;
    double H_e  = c.two_beta + mu_l2 * c.inv_sigma2;
    double inv_sqrt_H_e = 1.0 / std::sqrt(H_e);
    double mu_coef = -mu_l * c.inv_sigma2 / H_e;
    double sig_st2 = mu_l2 * c.inv_two_beta + c.sigma2;
    double inv_sig_st2 = 1.0 / sig_st2;
    double slab_const_saddle = c.log_sigma - 0.5 * std::log(sig_st2);
    double diag_const = std::log(2.0) + hk * log_beta - lgk
                        + c.half_log_2pi_sigma_diag2;
    double per_vertex = lgk - hk * log_beta;
    c.nu_chi_df[nu]                   = k;
    c.nu_half_k[nu]                   = hk;
    c.nu_km1[nu]                      = km1;
    c.nu_lgamma_half_k[nu]            = lgk;
    c.nu_mu_l[nu]                     = mu_l;
    c.nu_diag_const[nu]               = diag_const;
    c.nu_H_e_saddle[nu]               = H_e;
    c.nu_inv_sqrt_H_e_saddle[nu]      = inv_sqrt_H_e;
    c.nu_mu_coef_saddle[nu]           = mu_coef;
    c.nu_sigma_star_sq_saddle[nu]     = sig_st2;
    c.nu_inv_sigma_star_sq_saddle[nu] = inv_sig_st2;
    c.nu_slab_const_saddle[nu]        = slab_const_saddle;
    c.nu_per_vertex[nu]               = per_vertex;
  }
  return c;
}


PiAux make_pi_aux(const arma::imat& G_pi, const ChainAux& c) {
  int q = G_pi.n_rows;
  PiAux a;
  a.q = q;
  a.G_pi = G_pi;
  a.nu_pi.assign(q, 0);
  int E = 0;
  for (int r = 0; r < q - 1; ++r) {
    for (int s = r + 1; s < q; ++s) {
      if (G_pi(r, s) == 1) { a.nu_pi[r] += 1; ++E; }
    }
  }
  a.E_count = E;
  double log_beta = std::log(c.beta);
  double log_pi_over_beta = std::log(M_PI) - log_beta;
  double per_vertex = 0.0;
  for (int l = 0; l < q; ++l) per_vertex += c.nu_per_vertex[a.nu_pi[l]];
  double log_C0 = static_cast<double>(q) * log_beta
                + (-0.5 * static_cast<double>(E)) * std::log(2.0 * M_PI * c.sigma2)
                + per_vertex
                + 0.5 * static_cast<double>(E) * log_pi_over_beta;
  a.log_C0 = log_C0;
  return a;
}


arma::imat permute_graph(const arma::imat& G, const arma::ivec& pi) {
  int q = G.n_rows;
  arma::imat G_pi(q, q, arma::fill::zeros);
  for (int u = 0; u < q; ++u) {
    int pu = pi[u];
    for (int v = u + 1; v < q; ++v) {
      int pv = pi[v];
      if (G(u, v) == 1) {
        int rmin = std::min(pu, pv), rmax = std::max(pu, pv);
        G_pi(rmin, rmax) = 1;
        G_pi(rmax, rmin) = 1;
      }
    }
  }
  return G_pi;
}


arma::ivec degord_permutation(int q, int i, int j) {
  int lo = std::min(i, j), hi = std::max(i, j);
  arma::ivec pi(q, arma::fill::zeros);
  int next = 0;
  for (int v = 0; v < q; ++v) {
    if (v == lo || v == hi) continue;
    pi[v] = next++;
  }
  pi[lo] = q - 2;
  pi[hi] = q - 1;
  return pi;
}


double phi_pi_sample_from_noise(
    arma::mat& Phi_pi,
    arma::vec& row_logw,
    const double* noise,
    const PiAux& a,
    const ChainAux& c
) {
  int q = a.q;
  auto edge_offset = [q](int r, int s) -> int {
    return r * (q - 1) - r * (r - 1) / 2 + (s - r - 1);
  };
  double total_logw = 0.0;
  double half_inv_s2  = 0.5 * c.inv_sigma2;

  for (int r = 0; r < q - 1; ++r) {
    int nu_r = a.nu_pi[r];
    double mu_l_r        = c.nu_mu_l[nu_r];
    double km1_r         = c.nu_km1[nu_r];
    double diag_const_r  = c.nu_diag_const[nu_r];
    double slab_const_saddle_r        = c.nu_slab_const_saddle[nu_r];
    double inv_sqrt_H_e_r             = c.nu_inv_sqrt_H_e_saddle[nu_r];
    double mu_coef_saddle_r           = c.nu_mu_coef_saddle[nu_r];
    double inv_sigma_star_sq_saddle_r = c.nu_inv_sigma_star_sq_saddle[nu_r];

    double phi_rr = mu_l_r + c.sigma_diag * noise[r];
    if (phi_rr < 1e-12) phi_rr = 1e-12;
    Phi_pi(r, r) = phi_rr;

    double phi_rr2       = phi_rr * phi_rr;
    double sigma_star2   = phi_rr2 * c.inv_two_beta + c.sigma2;
    double inv_sigma_st2 = 1.0 / sigma_star2;
    double tau           = c.two_beta + phi_rr2 * c.inv_sigma2;
    double inv_sqrt_tau  = 1.0 / std::sqrt(tau);
    double mu_coef       = -phi_rr * c.inv_sigma2 / tau;
    double inv_phi_rr    = 1.0 / phi_rr;
    double slab_per_edge = c.log_sigma - 0.5 * std::log(sigma_star2);

    double dphi = phi_rr - mu_l_r;
    double diag_logw = diag_const_r
                     + km1_r * std::log(phi_rr)
                     - c.beta * phi_rr2
                     + 0.5 * dphi * dphi * c.inv_sigma_diag2;

    double phi_minus_mu = phi_rr - mu_l_r;
    double phi_plus_mu  = phi_rr + mu_l_r;

    double slab_const_for_row = (c.slab_tilt_mode == 1)
        ? slab_const_saddle_r
        : slab_per_edge;
    double rw = static_cast<double>(nu_r) * slab_const_for_row + diag_logw;

    const double* col_r = Phi_pi.colptr(r);
    for (int s = r + 1; s < q; ++s) {
      const double* col_s = Phi_pi.colptr(s);
      double S_rs = 0.0;
      for (int k = 0; k < r; ++k) S_rs += col_r[k] * col_s[k];
      if (a.G_pi(r, s) == 1) {
        int idx = q + edge_offset(r, s);
        double phi_rs;
        if (c.slab_tilt_mode == 1) {
          phi_rs = mu_coef_saddle_r * S_rs
                 + noise[idx] * inv_sqrt_H_e_r;
          Phi_pi(r, s) = phi_rs;
          rw -= 0.5 * S_rs * S_rs * inv_sigma_star_sq_saddle_r;
          rw -= phi_minus_mu * phi_rs
              * (phi_plus_mu * phi_rs + 2.0 * S_rs) * half_inv_s2;
        } else {
          phi_rs = mu_coef * S_rs + noise[idx] * inv_sqrt_tau;
          Phi_pi(r, s) = phi_rs;
          rw -= 0.5 * S_rs * S_rs * inv_sigma_st2;
        }
      } else {
        double phi_rs = -S_rs * inv_phi_rr;
        Phi_pi(r, s) = phi_rs;
        rw -= c.beta * phi_rs * phi_rs;
      }
    }
    row_logw[r] = rw;
    total_logw += rw;
  }

  // Trailing row q-1.
  {
    int r = q - 1;
    int nu_r = a.nu_pi[r];
    double mu_l_r       = c.nu_mu_l[nu_r];
    double km1_r        = c.nu_km1[nu_r];
    double diag_const_r = c.nu_diag_const[nu_r];
    double phi_rr = mu_l_r + c.sigma_diag * noise[r];
    if (phi_rr < 1e-12) phi_rr = 1e-12;
    Phi_pi(r, r) = phi_rr;
    double dphi = phi_rr - mu_l_r;
    double diag_logw = diag_const_r
                     + km1_r * std::log(phi_rr)
                     - c.beta * phi_rr * phi_rr
                     + 0.5 * dphi * dphi * c.inv_sigma_diag2;
    row_logw[r] = diag_logw;
    total_logw += diag_logw;
  }
  return total_logw;
}


double log_Zhat_pi_from_pool(
    const arma::mat& noise_pool_t,
    const PiAux& a,
    const ChainAux& c
) {
  int q = a.q;
  int M = static_cast<int>(noise_pool_t.n_cols);
  double neg_inf = -std::numeric_limits<double>::infinity();
  arma::vec log_w(M);
  log_w.fill(neg_inf);
  arma::mat Phi(q, q, arma::fill::zeros);
  arma::vec row_logw(q, arma::fill::zeros);
  double m = neg_inf;
  int n_finite = 0;
  for (int s = 0; s < M; ++s) {
    double lw = phi_pi_sample_from_noise(
        Phi, row_logw, noise_pool_t.colptr(s), a, c);
    if (std::isfinite(lw)) {
      log_w[s] = lw;
      ++n_finite;
      if (lw > m) m = lw;
    }
  }
  if (n_finite == 0) return neg_inf;
  double acc = 0.0;
  for (int s = 0; s < M; ++s)
    if (std::isfinite(log_w[s])) acc += std::exp(log_w[s] - m);
  return a.log_C0 + m + std::log(acc) - std::log(static_cast<double>(M));
}


double log_Zhat_pi_from_pool_cache(
    const arma::mat& noise_pool_t,
    const PiAux& a,
    const ChainAux& c,
    PoolCache& cache
) {
  int q = a.q;
  int M = static_cast<int>(noise_pool_t.n_cols);
  double neg_inf = -std::numeric_limits<double>::infinity();
  cache.log_w.set_size(M);    cache.log_w.fill(neg_inf);
  cache.rw_head.set_size(M);
  cache.S_trail.set_size(M);
  arma::mat Phi(q, q, arma::fill::zeros);
  arma::vec row_logw(q, arma::fill::zeros);
  double m = neg_inf;
  int n_finite = 0;
  for (int s = 0; s < M; ++s) {
    double lw = phi_pi_sample_from_noise(
        Phi, row_logw, noise_pool_t.colptr(s), a, c);
    if (std::isfinite(lw)) {
      cache.log_w[s] = lw;
      ++n_finite;
      if (lw > m) m = lw;
    }
    double rh = 0.0;
    for (int r = 0; r < q - 2; ++r) rh += row_logw[r];
    rh += row_logw[q - 1];   // r_qm1 is invariant across curr/star (nu_pi[q-1] = 0 always);
                             // include it in the head so delta_log_Zhat_pi_toggle's star
                             // aggregation matches direct log_Zhat(star) - log_Zhat(curr).
    cache.rw_head[s] = rh;
    double s_trail = 0.0;
    for (int k = 0; k < q - 2; ++k) s_trail += Phi(k, q - 2) * Phi(k, q - 1);
    cache.S_trail[s] = s_trail;
  }
  if (n_finite == 0) return neg_inf;
  double acc = 0.0;
  for (int s = 0; s < M; ++s)
    if (std::isfinite(cache.log_w[s])) acc += std::exp(cache.log_w[s] - m);
  return a.log_C0 + m + std::log(acc) - std::log(static_cast<double>(M));
}


double row_qm2_logw_from_S(
    double z_qm2,
    double z_trail,
    double S_trail,
    const PiAux& a,
    const ChainAux& c
) {
  int q = a.q;
  int r = q - 2;
  int nu_r = a.nu_pi[r];
  double mu_l_r       = c.nu_mu_l[nu_r];
  double km1_r        = c.nu_km1[nu_r];
  double diag_const_r = c.nu_diag_const[nu_r];
  double phi_rr = mu_l_r + c.sigma_diag * z_qm2;
  if (phi_rr < 1e-12) phi_rr = 1e-12;
  double phi_rr2 = phi_rr * phi_rr;
  double dphi = phi_rr - mu_l_r;
  double diag_logw = diag_const_r
                   + km1_r * std::log(phi_rr)
                   - c.beta * phi_rr2
                   + 0.5 * dphi * dphi * c.inv_sigma_diag2;
  if (a.G_pi(r, q - 1) == 1) {
    if (c.slab_tilt_mode == 1) {
      double mu_coef_saddle_r           = c.nu_mu_coef_saddle[nu_r];
      double inv_sqrt_H_e_r             = c.nu_inv_sqrt_H_e_saddle[nu_r];
      double slab_const_saddle_r        = c.nu_slab_const_saddle[nu_r];
      double inv_sigma_star_sq_saddle_r = c.nu_inv_sigma_star_sq_saddle[nu_r];
      double mu_rs   = mu_coef_saddle_r * S_trail;
      double phi_rs  = mu_rs + z_trail * inv_sqrt_H_e_r;
      double slab_logw = slab_const_saddle_r
                       - 0.5 * S_trail * S_trail * inv_sigma_star_sq_saddle_r;
      double phi_minus_mu = phi_rr - mu_l_r;
      double phi_plus_mu  = phi_rr + mu_l_r;
      double mismatch     = phi_minus_mu * phi_rs
                          * (phi_plus_mu * phi_rs + 2.0 * S_trail)
                          * 0.5 * c.inv_sigma2;
      return diag_logw + slab_logw - mismatch;
    } else {
      (void) z_trail;
      double sigma_star2 = phi_rr2 * c.inv_two_beta + c.sigma2;
      return diag_logw
           + c.log_sigma - 0.5 * std::log(sigma_star2)
           - 0.5 * S_trail * S_trail / sigma_star2;
    }
  } else {
    (void) z_trail;
    double phi_rs = -S_trail / phi_rr;
    return diag_logw - c.beta * phi_rs * phi_rs;
  }
}


double delta_log_Zhat_pi_toggle(
    const arma::mat& noise_pool,
    const arma::mat& noise_pool_t,
    const arma::imat& G_curr,
    int i, int j,
    const ChainAux& c
) {
  int q = c.q;
  double neg_inf = -std::numeric_limits<double>::infinity();
  arma::ivec pi = degord_permutation(q, i, j);
  arma::imat G_pi_curr = permute_graph(G_curr, pi);
  arma::imat G_pi_star = G_pi_curr;
  int toggled = 1 - G_pi_curr(q - 2, q - 1);
  G_pi_star(q - 2, q - 1) = toggled;
  G_pi_star(q - 1, q - 2) = toggled;

  PiAux a_curr = make_pi_aux(G_pi_curr, c);
  PiAux a_star = make_pi_aux(G_pi_star, c);

  PoolCache cache_curr;
  double log_Zhat_curr = log_Zhat_pi_from_pool_cache(
      noise_pool_t, a_curr, c, cache_curr);

  int M = static_cast<int>(noise_pool.n_rows);
  // z_qm2 access is contiguous along s for fixed d = q-2 in the M x dim layout.
  const double* col_qm2 = noise_pool.colptr(q - 2);
  // z_trail: slab noise slot for the trailing off-diagonal (q-2, q-1). The
  // off-diagonal section begins at offset q in the noise vector and uses
  // edge_offset(r, s) = r*(q-1) - r*(r-1)/2 + (s - r - 1). At
  // (r, s) = (q-2, q-1) this collapses to (q-2)*(q+1)/2; total slab index
  // is q + (q-2)*(q+1)/2. Used only when slab_tilt_mode == 1 (saddle-shifted
  // IS path inside row_qm2_logw_from_S); slab_tilt_mode = 0 path ignores
  // z_trail via (void) z_trail. Hardcoding z_trail = 0.0 here was the v4
  // mode-1 sibling of the rw_head-misses-row-q-1 bug.
  int slab_idx = q + (q - 2) * (q + 1) / 2;
  const double* col_slab = noise_pool.colptr(slab_idx);
  arma::vec log_w_star(M);
  log_w_star.fill(neg_inf);
  double m = neg_inf;
  int n_finite = 0;
  for (int s = 0; s < M; ++s) {
    double z_qm2 = col_qm2[s];
    double z_trail = col_slab[s];
    double rw_qm2_star = row_qm2_logw_from_S(
        z_qm2, z_trail, cache_curr.S_trail[s], a_star, c);
    double total = cache_curr.rw_head[s] + rw_qm2_star;
    if (std::isfinite(total)) {
      log_w_star[s] = total;
      ++n_finite;
      if (total > m) m = total;
    }
  }
  double log_Zhat_star;
  if (n_finite == 0) {
    log_Zhat_star = neg_inf;
  } else {
    double acc = 0.0;
    for (int s = 0; s < M; ++s)
      if (std::isfinite(log_w_star[s])) acc += std::exp(log_w_star[s] - m);
    log_Zhat_star = a_star.log_C0 + m + std::log(acc)
                  - std::log(static_cast<double>(M));
  }
  return log_Zhat_star - log_Zhat_curr;
}


arma::mat draw_bartlett_pool(SafeRNG& rng, int q, int M_inner) {
  int dim = bartlett_pool_dim(q);
  return arma_rnorm_mat(rng, static_cast<arma::uword>(dim),
                        static_cast<arma::uword>(M_inner));
}


}  // namespace degord
