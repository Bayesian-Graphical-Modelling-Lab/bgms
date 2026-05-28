#pragma once

#include <RcppArmadillo.h>

// Savage-Dickey 1D conditional density on l_ji (DEGORD-permuted trailing
// Cholesky entry), evaluated at a user-supplied point x_eval (typically the
// Roverato slave m_ij). The log-kernel is
//
//   f(x) = -A x² + B x + (alpha - 1) log(s_jj + x²)
//
// where, with K = L Lᵀ, (i, j) permuted to (p-2, p-1):
//
//   A    = (l_ii² / sigma² + 2 beta + S_jj_data) / 2
//   B    = l_ii² m_ij / sigma²  -  S_ij_data l_ii
//   s_jj = K_jj_old - l_ji_old²  (the l_ji-independent rest of K_jj)
//
// Encompassing prior: Normal slab N(0, sigma²) on every off-diagonal K_kl
// plus Gamma(alpha, beta) on the diagonal halves K_kk/2 and the determinant
// tilt |K|^delta. l_ji integrates over R (no PD truncation in L-space).
//
// At alpha = 1 the (alpha - 1) log term vanishes and f is a pure Gaussian
// in x; mode = B/(2A), curvature = 2A, Laplace is exact, NLO terms are
// identically zero. The primitive recovers the closed-form Gaussian.
//
// Status codes:
//   0  ok
//   2  curvature non-positive (Laplace invalid at mode; chain should revert)
//   3  Newton did not converge within newton_max_iter

namespace savage_dickey {

struct LSDResult {
    double log_density;   // log pi(x_eval | rest, Y) via Laplace + optional NLO
    double x_mode;
    double curvature;     // -f''(x_mode)
    int    status;
};

LSDResult density_at_l_ji_one(
    double x_eval,
    double A,
    double B,
    double s_jj,
    double alpha,
    bool   nlo            = true,
    int    newton_max_iter = 50,
    double newton_tol      = 1e-10);

}  // namespace savage_dickey
