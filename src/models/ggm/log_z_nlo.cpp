// Closed-form log Z(G) approximations for the spikeslab GGM prior with
// Gamma diagonal + determinant tilt. See log_z_nlo.h for derivation
// references; algorithm is a direct port of branchB_chain.cpp:470-620 and
// incremental_log_Z_NLO_gamma.h in ~/SV/Z/R/src/.

#include "models/ggm/log_z_nlo.h"
#include <vector>
#include <cmath>
#include <limits>

// ----------------------------------------------------------------------
// Local helpers (alpha = 1 specialisations for the incremental form).
// ----------------------------------------------------------------------
static inline int nu_at(int l, const arma::imat& G) {
  int q = G.n_rows;
  int n = 0;
  for (int m = l + 1; m < q; ++m) if (G(l, m) == 1) n += 1;
  return n;
}

static inline double Mv_at_alpha1(int l, const arma::imat& G, double delta) {
  return static_cast<double>(nu_at(l, G)) + 2.0 * (1.0 + delta) - 1.0;
}

static inline double He_alpha1(int k, const arma::imat& G,
                                double beta, double sigma, double delta) {
  double Mv_k = Mv_at_alpha1(k, G, delta);
  return 2.0 * beta + Mv_k / (2.0 * sigma * sigma * beta);
}

// ----------------------------------------------------------------------
// Full-recompute log Z_LO + log Z_NLO at general alpha, with optional
// NNLO_H F-piece (off by default; on alpha > 1 the Exp-style F over-
// corrects). At delta = 0 collapses to the pre-tilt formula bit-exact.
// ----------------------------------------------------------------------
double log_Z_NLO_gamma(
    const arma::imat& G,
    double alpha, double beta, double sigma,
    bool include_F,
    double delta
) {
  int q = G.n_rows;
  std::vector<int> nu(q, 0);
  int E_count = 0;
  for (int l = 0; l < q; ++l)
    for (int m = l + 1; m < q; ++m)
      if (G(l, m) == 1) { ++nu[l]; ++E_count; }

  double sigma2 = sigma * sigma;

  // Diagonal Gamma integral per vertex + tilt-induced constants.
  // Tilt shifts the lgamma argument by 2 delta and adds -q delta log beta
  // to the diagonal piece summed over the q vertices.
  double log_C0 = 0.5 * static_cast<double>(E_count) * std::log(M_PI)
                  - 0.5 * static_cast<double>(E_count) * std::log(2.0 * M_PI * sigma2)
                  - static_cast<double>(E_count) * std::log(beta)
                  - static_cast<double>(q) * std::lgamma(alpha)
                  - static_cast<double>(q) * delta * std::log(beta);
  for (int l = 0; l < q; ++l)
    log_C0 += std::lgamma((static_cast<double>(nu[l]) + 2.0 * (alpha + delta)) / 2.0);

  if (E_count == 0) return log_C0;

  // Tilt-shifted saddle exponent M_v = nu_v + 2 (alpha + delta) - 1.
  std::vector<double> M_v(q);
  for (int l = 0; l < q; ++l)
    M_v[l] = static_cast<double>(nu[l]) + 2.0 * (alpha + delta) - 1.0;

  // H_e per edge (k, v), k < v:
  //   H_e = 2 beta + M_k / (2 sigma^2 beta) - 4 beta (alpha - 1) / M_v.
  arma::mat H_e_lookup(q, q, arma::fill::zeros);
  for (int k = 0; k < q; ++k)
    for (int v = k + 1; v < q; ++v)
      if (G(k, v) == 1) {
        double h = 2.0 * beta
                   + M_v[k] / (2.0 * sigma2 * beta)
                   - 4.0 * beta * (alpha - 1.0) / M_v[v];
        H_e_lookup(k, v) = h;
        H_e_lookup(v, k) = h;
      }

  // LO slab: 0.5 sum_e log(2 beta / H_e). Guard against H_e <= 0.
  double log_LO = 0.0;
  for (int k = 0; k < q; ++k)
    for (int v = k + 1; v < q; ++v)
      if (G(k, v) == 1) {
        double h = H_e_lookup(k, v);
        if (h <= 0.0) return -std::numeric_limits<double>::infinity();
        log_LO += 0.5 * std::log(2.0 * beta / h);
      }

  // T (outgoing edges, v as smaller endpoint) and S (incoming).
  std::vector<double> T1(q, 0.0), T2(q, 0.0), S1(q, 0.0), S2(q, 0.0);
  for (int v = 0; v < q; ++v) {
    for (int w = v + 1; w < q; ++w)
      if (G(v, w) == 1) {
        double inv = 1.0 / H_e_lookup(v, w);
        T1[v] += inv; T2[v] += inv * inv;
      }
    for (int k = 0; k < v; ++k)
      if (G(k, v) == 1) {
        double inv = 1.0 / H_e_lookup(k, v);
        S1[v] += inv; S2[v] += inv * inv;
      }
  }

  double delta_NLO = 0.0;

  // B_e per edge: -(nu_i + alpha) / (4 beta sigma^2 M_i H_e), e = (i, j).
  for (int i = 0; i < q; ++i)
    for (int j = i + 1; j < q; ++j)
      if (G(i, j) == 1)
        delta_NLO -= (static_cast<double>(nu[i]) + alpha)
                     / (4.0 * beta * sigma2 * M_v[i] * H_e_lookup(i, j));

  // C_{e, k} per (edge e=(i,j), common predecessor k < i with G(k,i)=G(k,j)=1).
  for (int i = 0; i < q; ++i) {
    for (int j = i + 1; j < q; ++j) {
      if (G(i, j) != 1) continue;
      double H_e = H_e_lookup(i, j);
      for (int k = 0; k < i; ++k)
        if (G(k, i) == 1 && G(k, j) == 1) {
          double H_ki = H_e_lookup(k, i), H_kj = H_e_lookup(k, j);
          delta_NLO += -1.0 / (2.0 * sigma2 * H_ki * H_kj)
                       + M_v[i] / (4.0 * beta * sigma2 * sigma2 * H_e * H_ki * H_kj);
        }
    }
  }

  // D_v per vertex with at least one outgoing edge.
  for (int v = 0; v < q; ++v)
    if (T1[v] > 0.0)
      delta_NLO += M_v[v] / (16.0 * beta * beta * sigma2 * sigma2)
                   * (2.0 * T2[v] + T1[v] * T1[v]);

  // E_{lm} (+ optional F_{lm}) per non-edge with common predecessors.
  for (int l = 1; l < q - 1; ++l) {
    for (int m = l + 1; m < q; ++m) {
      if (G(l, m) != 0) continue;
      double s1 = 0.0, s2 = 0.0;
      bool has_common = false;
      for (int k = 0; k < l; ++k)
        if (G(k, l) == 1 && G(k, m) == 1) {
          double H_kl = H_e_lookup(k, l), H_km = H_e_lookup(k, m);
          double inv_prod = 1.0 / (H_kl * H_km);
          s1 += inv_prod; s2 += inv_prod * inv_prod;
          has_common = true;
        }
      if (!has_common) continue;
      double Ml = M_v[l], b2 = beta * beta, b4 = b2 * b2;
      delta_NLO -= 2.0 * b2 / Ml * s1;
      if (include_F) {
        delta_NLO -= 2.0 * b2 / (Ml * Ml) * s1;
        delta_NLO += 12.0 * b4 / (Ml * Ml) * s2;
        delta_NLO += 4.0 * b4 / (Ml * Ml) * s1 * s1;
      }
    }
  }

  // N_v per vertex (alpha > 1 Gamma piece; zero at alpha = 1).
  if (alpha != 1.0) {
    double am1 = alpha - 1.0;
    for (int v = 0; v < q; ++v) {
      double Mv = M_v[v];
      double bracket_4 = 2.0 * S2[v] + S1[v] * S1[v];
      double Mv2 = Mv * Mv, Mv3 = Mv2 * Mv;
      double Na = -am1 * ( 3.0 / (8.0 * Mv2)
                          - 3.0 * beta * S1[v] / Mv2
                          + 2.0 * beta * beta * bracket_4 / Mv2 );
      double Nb = am1 * am1 * ( 5.0 / (12.0 * Mv3)
                                - 2.0 * beta * S1[v] / Mv3
                                + 4.0 * beta * beta * bracket_4 / Mv3 );
      double Nc = am1 * (static_cast<double>(nu[v]) + 1.0)
                  * ( 5.0 / (12.0 * Mv3) - beta * S1[v] / Mv3 );
      delta_NLO += Na + Nb + Nc;
    }
    // M_e per edge.
    for (int i = 0; i < q; ++i)
      for (int j = i + 1; j < q; ++j)
        if (G(i, j) == 1)
          delta_NLO += am1 / (sigma2 * M_v[i] * H_e_lookup(i, j))
                       * (-1.0 / (4.0 * beta) + S1[i]);
  }

  return log_C0 + log_LO + delta_NLO;
}

// ----------------------------------------------------------------------
// DEGORD reordering: permute toggle endpoints to positions (0, 1) and
// keep all other vertices in their original order. The closed-form Z
// is not permutation-invariant, so MH ratios computing log Z(G+)-log Z(G-)
// must use the same reordering in both terms.
// ----------------------------------------------------------------------
double log_Z_NLO_gamma_degord(
    const arma::imat& G, int i, int j,
    double alpha, double beta, double sigma,
    bool include_F, double delta
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
  return log_Z_NLO_gamma(G_perm, alpha, beta, sigma, include_F, delta);
}

// ----------------------------------------------------------------------
// alpha = 1 incremental form. Sum of all log Z_NLO contributions whose
// value depends on (nu_i, M_v[i], H_e for edges adjacent to i), evaluated
// at vertex i in graph G. At alpha = 1 a single-edge toggle (i, j) with
// i < j changes only these terms (the §4.6 locality result).
// ----------------------------------------------------------------------
static double vertex_i_contribs_alpha1(
    const arma::imat& G, int i,
    double beta, double sigma, double delta,
    bool include_F
) {
  int q = G.n_rows;
  double sigma2 = sigma * sigma;
  double sigma4 = sigma2 * sigma2;
  double beta2  = beta * beta;
  double beta4  = beta2 * beta2;

  std::vector<int> fwd_i;
  for (int v = i + 1; v < q; ++v) if (G(i, v) == 1) fwd_i.push_back(v);
  int nu_i = static_cast<int>(fwd_i.size());
  double Mv_i = static_cast<double>(nu_i) + 2.0 * (1.0 + delta) - 1.0;
  double H_i  = 2.0 * beta + Mv_i / (2.0 * sigma2 * beta);

  double total = 0.0;

  // LO + C0 contributions.
  double C0_per_edge = 0.5 * std::log(M_PI)
                       - 0.5 * std::log(2.0 * M_PI * sigma2)
                       - std::log(beta);
  total += static_cast<double>(nu_i) * C0_per_edge;
  total += std::lgamma((static_cast<double>(nu_i) + 2.0 * (1.0 + delta)) / 2.0);
  if (nu_i > 0)
    total += 0.5 * static_cast<double>(nu_i) * std::log(2.0 * beta / H_i);

  // Block 1 slab term per forward edge of i:
  // -(nu_i + 1) / (4 beta sigma^2 M_v[i] H_e(i, v)). At alpha = 1, H_e(i, v) = H_i.
  double block1_per_edge = -(static_cast<double>(nu_i) + 1.0)
                            / (4.0 * beta * sigma2 * Mv_i * H_i);
  total += static_cast<double>(nu_i) * block1_per_edge;

  // Block 2 cross-predecessor per (forward edge (i, v), common pred k < i).
  for (int v_idx = 0; v_idx < nu_i; ++v_idx) {
    int v = fwd_i[v_idx];
    for (int k = 0; k < i; ++k) {
      if (G(k, i) == 1 && G(k, v) == 1) {
        double H_ki = He_alpha1(k, G, beta, sigma, delta);
        double H_kv = He_alpha1(k, G, beta, sigma, delta);
        total += -1.0 / (2.0 * sigma2 * H_ki * H_kv)
                 + Mv_i / (4.0 * beta * sigma4 * H_i * H_ki * H_kv);
      }
    }
  }

  // Block 3 D_v at v = i. At alpha = 1, H_e(i, w) = H_i for all forward
  // edges of i, so T1[i] = nu_i / H_i, T2[i] = nu_i / H_i^2.
  {
    double T1_i = static_cast<double>(nu_i) / H_i;
    double T2_i = static_cast<double>(nu_i) / (H_i * H_i);
    total += Mv_i / (16.0 * beta2 * sigma4) * (2.0 * T2_i + T1_i * T1_i);
  }

  // Block 4 non-edges (i, m) where i is the smaller endpoint.
  if (i > 0) {
    for (int m = i + 1; m < q; ++m) {
      if (G(i, m) != 0) continue;
      double s1 = 0.0, s2 = 0.0;
      bool any = false;
      for (int k = 0; k < i; ++k) {
        if (G(k, i) == 1 && G(k, m) == 1) {
          double H_ki = He_alpha1(k, G, beta, sigma, delta);
          double H_km = He_alpha1(k, G, beta, sigma, delta);
          double inv = 1.0 / (H_ki * H_km);
          s1 += inv;
          s2 += inv * inv;
          any = true;
        }
      }
      if (!any) continue;
      total += -2.0 * beta2 / Mv_i * s1;
      if (include_F) {
        total += -2.0 * beta2 / (Mv_i * Mv_i) * s1;
        total += 12.0 * beta4 / (Mv_i * Mv_i) * s2;
        total += 4.0 * beta4 / (Mv_i * Mv_i) * s1 * s1;
      }
    }
  }

  // Block 2 (k = i case) AND Block 4 (k = i case): pairs of forward
  // neighbours (a, b) of i with i < a < b.
  if (nu_i >= 2) {
    for (int idx_a = 0; idx_a < nu_i - 1; ++idx_a) {
      int a = fwd_i[idx_a];
      double Mv_a = Mv_at_alpha1(a, G, delta);
      double H_ia = He_alpha1(i, G, beta, sigma, delta);
      for (int idx_b = idx_a + 1; idx_b < nu_i; ++idx_b) {
        int b = fwd_i[idx_b];
        double H_ib = He_alpha1(i, G, beta, sigma, delta);
        if (G(a, b) == 1) {
          double H_ab = He_alpha1(a, G, beta, sigma, delta);
          total += -1.0 / (2.0 * sigma2 * H_ia * H_ib)
                   + Mv_a / (4.0 * beta * sigma4 * H_ab * H_ia * H_ib);
        } else {
          double v_ab = 1.0 / (H_ia * H_ib);
          double sum_other = 0.0;
          for (int k = 0; k < a; ++k) {
            if (k == i) continue;
            if (G(k, a) == 1 && G(k, b) == 1) {
              double H_ka = He_alpha1(k, G, beta, sigma, delta);
              double H_kb = He_alpha1(k, G, beta, sigma, delta);
              sum_other += 1.0 / (H_ka * H_kb);
            }
          }
          total += -2.0 * beta2 / Mv_a * v_ab;
          if (include_F) {
            total += -2.0 * beta2 / (Mv_a * Mv_a) * v_ab;
            total += 12.0 * beta4 / (Mv_a * Mv_a) * (v_ab * v_ab);
            total += 4.0 * beta4 / (Mv_a * Mv_a) * (v_ab * v_ab);
            total += 4.0 * beta4 / (Mv_a * Mv_a) * 2.0 * v_ab * sum_other;
          }
        }
      }
    }
  }

  return total;
}

// Public alpha = 1 incremental log Z_NLO ratio under single-edge toggle (i, j).
// Computes the before-side first so that any accidental aliasing between the
// two graphs (Armadillo copy-on-write) does not corrupt the before evaluation.
double log_Z_NLO_gamma_delta_incr_alpha1(
    const arma::imat& G_before, int i, int j,
    double beta, double sigma, double delta,
    bool include_F
) {
  int i_min = std::min(i, j);
  double before_i = vertex_i_contribs_alpha1(
      G_before, i_min, beta, sigma, delta, include_F);
  arma::imat G_after(G_before);
  G_after(i, j) = 1 - G_before(i, j);
  G_after(j, i) = G_after(i, j);
  double after_i = vertex_i_contribs_alpha1(
      G_after, i_min, beta, sigma, delta, include_F);
  return after_i - before_i;
}

// ----------------------------------------------------------------------
// Partial log Z_NLO on an affected-vertex set V_a.
//
// Returns the sum of all log_Z_NLO_gamma term instances whose index set
// intersects V_a. Instances disjoint from V_a are invariant under any
// toggle that only changes edges incident to V_a, so they cancel in the
// difference partial(G_after) - partial(G_before).
//
// Loop bounds are restricted to V_a's neighbourhood, using a canonical-
// representative rule (the lowest-index V_a member of an instance owns
// its enumeration) to count each qualifying instance exactly once. M_v,
// H_e, and per-vertex T/S sums are precomputed on the full graph in
// O(|E|) once the adjacency lists are built; the per-term loops are then
// O(deg^2 + deg * q) instead of O(q^2 + q^3).
//
// Inputs: in_V (size q membership bitmap) and V_a_idx (V_a indices in
// ascending order) — caller supplies both.
// ----------------------------------------------------------------------
static double partial_log_Z_NLO_gamma_on_V(
    const arma::imat& G,
    const std::vector<bool>& in_V,
    const std::vector<int>& V_a_idx,
    double alpha, double beta, double sigma,
    bool include_F, double delta
) {
  int q = G.n_rows;
  double sigma2 = sigma * sigma;
  double sigma4 = sigma2 * sigma2;
  double beta2  = beta * beta;
  double beta4  = beta2 * beta2;

  // Adjacency lists for fast neighbour iteration. O(|E|) to build.
  std::vector<std::vector<int>> fwd(q), bwd(q);
  for (int l = 0; l < q; ++l)
    for (int m = l + 1; m < q; ++m)
      if (G(l, m) == 1) {
        fwd[l].push_back(m);
        bwd[m].push_back(l);
      }

  std::vector<int> nu(q);
  for (int l = 0; l < q; ++l) nu[l] = static_cast<int>(fwd[l].size());

  // Count edges incident to V_a via the canonical-rep rule (each edge owned
  // by its lowest-index V_a endpoint).
  int E_count_in_V = 0;
  for (int p : V_a_idx) {
    E_count_in_V += static_cast<int>(fwd[p].size());        // (p, m), p smaller
    for (int l : bwd[p]) if (!in_V[l]) ++E_count_in_V;       // (l, p), l smaller, l not in V_a
  }
  int q_in_V = static_cast<int>(V_a_idx.size());

  // log_C0 partial.
  double log_C0 = 0.5 * static_cast<double>(E_count_in_V) * std::log(M_PI)
                  - 0.5 * static_cast<double>(E_count_in_V) * std::log(2.0 * M_PI * sigma2)
                  - static_cast<double>(E_count_in_V) * std::log(beta)
                  - static_cast<double>(q_in_V) * std::lgamma(alpha)
                  - static_cast<double>(q_in_V) * delta * std::log(beta);
  for (int p : V_a_idx)
    log_C0 += std::lgamma((static_cast<double>(nu[p]) + 2.0 * (alpha + delta)) / 2.0);

  // Total edge count: needed to short-circuit when graph is empty.
  int E_count_total = 0;
  for (int l = 0; l < q; ++l) E_count_total += nu[l];
  if (E_count_total == 0) return log_C0;

  // M_v[l] for all l (any v adjacent to V_a can appear as an index of an
  // included term).
  std::vector<double> M_v(q);
  for (int l = 0; l < q; ++l)
    M_v[l] = static_cast<double>(nu[l]) + 2.0 * (alpha + delta) - 1.0;

  // Full H_e lookup (some non-V_a-incident edges are still cross-referenced
  // in C_{e=(a, b), k} triples where neither k nor a is in V_a but b is).
  arma::mat H_e_lookup(q, q, arma::fill::zeros);
  for (int k = 0; k < q; ++k)
    for (int v : fwd[k]) {
      double h = 2.0 * beta
                 + M_v[k] / (2.0 * sigma2 * beta)
                 - 4.0 * beta * (alpha - 1.0) / M_v[v];
      H_e_lookup(k, v) = h;
      H_e_lookup(v, k) = h;
    }

  // T (forward, v as smaller) and S (backward, v as larger). For M_e in
  // the alpha > 1 branch we need S1[a] at edges (a, b) where a is a
  // backward neighbour of some V_a vertex (so a may be outside V_a). Keep
  // T/S for all q vertices; cost O(|E|).
  std::vector<double> T1(q, 0.0), T2(q, 0.0), S1(q, 0.0), S2(q, 0.0);
  for (int v = 0; v < q; ++v) {
    for (int w : fwd[v]) {
      double inv = 1.0 / H_e_lookup(v, w);
      T1[v] += inv; T2[v] += inv * inv;
    }
    for (int k : bwd[v]) {
      double inv = 1.0 / H_e_lookup(k, v);
      S1[v] += inv; S2[v] += inv * inv;
    }
  }

  double log_LO = 0.0;
  double delta_NLO = 0.0;

  // ---- Edge-indexed terms (log_LO, B_e, M_e): enumerate V_a-incident
  // edges via canonical rep (smaller V_a endpoint owns).
  bool alpha_nontrivial = (alpha != 1.0);
  double am1 = alpha - 1.0;

  auto process_edge = [&](int a, int b) {
    double h = H_e_lookup(a, b);
    if (h <= 0.0) {
      log_LO = -std::numeric_limits<double>::infinity();
      return;
    }
    log_LO    += 0.5 * std::log(2.0 * beta / h);
    delta_NLO -= (static_cast<double>(nu[a]) + alpha)
                 / (4.0 * beta * sigma2 * M_v[a] * h);
    if (alpha_nontrivial)
      delta_NLO += am1 / (sigma2 * M_v[a] * h)
                   * (-1.0 / (4.0 * beta) + S1[a]);
  };

  for (int p : V_a_idx) {
    // (p, m) with p smaller: p owns.
    for (int m : fwd[p]) process_edge(p, m);
    if (!std::isfinite(log_LO)) return log_LO;
    // (l, p) with l < p smaller, l not in V_a: p owns.
    for (int l : bwd[p])
      if (!in_V[l]) process_edge(l, p);
    if (!std::isfinite(log_LO)) return log_LO;
  }

  // ---- C_{(a, b), k}: triples (k, a, b), k < a < b, all three edges
  // present. Canonical rep = lowest-index V_a member of {k, a, b}.
  auto process_C_triple = [&](int k, int a, int b) {
    double H_ka = H_e_lookup(k, a);
    double H_kb = H_e_lookup(k, b);
    double H_ab = H_e_lookup(a, b);
    delta_NLO += -1.0 / (2.0 * sigma2 * H_ka * H_kb)
                 + M_v[a] / (4.0 * beta * sigma4 * H_ab * H_ka * H_kb);
  };

  for (int p : V_a_idx) {
    // Case p = k: triples (p, a, b), p < a < b, all edges present. p is the
    // smallest of the triple so p is canon (lowest V_a member).
    int nf = static_cast<int>(fwd[p].size());
    for (int ai = 0; ai < nf; ++ai) {
      int a = fwd[p][ai];
      for (int bi = ai + 1; bi < nf; ++bi) {
        int b = fwd[p][bi];
        if (G(a, b) == 1) process_C_triple(p, a, b);
      }
    }
    // Case p = a: triples (k, p, b), k < p < b. Canon = p iff k not in V_a.
    for (int k : bwd[p]) {
      if (in_V[k]) continue;
      for (int b : fwd[p])
        if (G(k, b) == 1) process_C_triple(k, p, b);
    }
    // Case p = b: triples (k, a, p), k < a < p. Canon = p iff neither k nor
    // a is in V_a.
    int nb = static_cast<int>(bwd[p].size());
    for (int ai = 0; ai < nb; ++ai) {
      int a = bwd[p][ai];
      if (in_V[a]) continue;
      for (int ki = 0; ki < ai; ++ki) {
        int k = bwd[p][ki];
        if (in_V[k]) continue;
        if (G(k, a) == 1) process_C_triple(k, a, p);
      }
    }
  }

  // ---- D_v: vertex term. v in V_a, T1[v] > 0.
  for (int p : V_a_idx)
    if (T1[p] > 0.0)
      delta_NLO += M_v[p] / (16.0 * beta2 * sigma4)
                   * (2.0 * T2[p] + T1[p] * T1[p]);

  // ---- E_{lm} (+ F): non-edges (l, m), l < m, with at least one common
  // predecessor. Canonical rep = lowest-index V_a member of {l, m}.
  auto process_E_pair = [&](int l, int m) {
    double s1 = 0.0, s2 = 0.0;
    bool has_common = false;
    // Common predecessors: k < l with G(k, l) = G(k, m) = 1. Use bwd[l]
    // (predecessors of l) and check G(k, m).
    for (int k : bwd[l]) {
      if (G(k, m) == 1) {
        double H_kl = H_e_lookup(k, l), H_km = H_e_lookup(k, m);
        double inv_prod = 1.0 / (H_kl * H_km);
        s1 += inv_prod; s2 += inv_prod * inv_prod;
        has_common = true;
      }
    }
    if (!has_common) return;
    double Ml = M_v[l];
    delta_NLO -= 2.0 * beta2 / Ml * s1;
    if (include_F) {
      delta_NLO -= 2.0 * beta2 / (Ml * Ml) * s1;
      delta_NLO += 12.0 * beta4 / (Ml * Ml) * s2;
      delta_NLO += 4.0 * beta4 / (Ml * Ml) * s1 * s1;
    }
  };

  for (int p : V_a_idx) {
    // Non-edges (p, m) with p smaller, m > p.
    for (int m = p + 1; m < q; ++m)
      if (G(p, m) == 0) process_E_pair(p, m);
    // Non-edges (l, p) with l < p smaller, l not in V_a.
    for (int l = 0; l < p; ++l)
      if (G(l, p) == 0 && !in_V[l]) process_E_pair(l, p);
  }

  // ---- N_v (alpha > 1 only): vertex term at v in V_a. The N_v formula
  // also depends on S1[v]/S2[v], which we already have for all v.
  if (alpha_nontrivial) {
    for (int p : V_a_idx) {
      double Mv = M_v[p];
      double bracket_4 = 2.0 * S2[p] + S1[p] * S1[p];
      double Mv2 = Mv * Mv, Mv3 = Mv2 * Mv;
      double Na = -am1 * ( 3.0 / (8.0 * Mv2)
                          - 3.0 * beta * S1[p] / Mv2
                          + 2.0 * beta2 * bracket_4 / Mv2 );
      double Nb = am1 * am1 * ( 5.0 / (12.0 * Mv3)
                                - 2.0 * beta * S1[p] / Mv3
                                + 4.0 * beta2 * bracket_4 / Mv3 );
      double Nc = am1 * (static_cast<double>(nu[p]) + 1.0)
                  * ( 5.0 / (12.0 * Mv3) - beta * S1[p] / Mv3 );
      delta_NLO += Na + Nb + Nc;
    }
  }

  return log_C0 + log_LO + delta_NLO;
}

// Public general-alpha incremental log Z_NLO ratio under single-edge toggle.
double log_Z_NLO_gamma_delta_incr_alphaN(
    const arma::imat& G_before, int i, int j,
    double alpha, double beta, double sigma, double delta,
    bool include_F
) {
  int q = G_before.n_rows;

  // Affected vertex set V_a = {i, j} U N_G_before(i). N_after(i) differs
  // from N_before(i) by at most {j}, and j is in V_a regardless.
  std::vector<bool> in_V(q, false);
  in_V[i] = true;
  in_V[j] = true;
  for (int v = 0; v < q; ++v)
    if (v != i && G_before(i, v) == 1)
      in_V[v] = true;
  std::vector<int> V_a_idx;
  V_a_idx.reserve(q);
  for (int v = 0; v < q; ++v) if (in_V[v]) V_a_idx.push_back(v);

  // Evaluate partial sum on G_before first so any accidental aliasing
  // between G_before and G_after (Armadillo copy-on-write) cannot corrupt
  // the before evaluation.
  double before = partial_log_Z_NLO_gamma_on_V(
      G_before, in_V, V_a_idx, alpha, beta, sigma, include_F, delta);
  arma::imat G_after(G_before);
  G_after(i, j) = 1 - G_before(i, j);
  G_after(j, i) = G_after(i, j);
  double after = partial_log_Z_NLO_gamma_on_V(
      G_after, in_V, V_a_idx, alpha, beta, sigma, include_F, delta);
  return after - before;
}
