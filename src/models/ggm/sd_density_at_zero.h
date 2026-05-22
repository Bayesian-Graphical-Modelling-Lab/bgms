#pragma once

#include <RcppArmadillo.h>

// Savage-Dickey 1D conditional density of K_ij at zero, given K_{-ij} and
// data sufficient stats (S, n). The encompassing posterior is the slab-only
// joint with determinant tilt delta and Gaussian slab N(0, sigma^2) on every
// off-diagonal. As a function of x = K_ij with K_{-ij} fixed:
//
//   L(x) = (delta + n/2) log[1 + 2 c1 x - c2 x^2] - S_ij x - x^2/(2 sigma^2),
//
// where c1 = Sigma_0[i,j], c2 = Sigma_0[i,i] Sigma_0[j,j] - Sigma_0[i,j]^2,
// Sigma_0 = (K with K_ij set to zero)^{-1}.
//
// log pi_enc(K_ij = 0 | K_{-ij}, Y) is then ell(0) - log Z_1D, evaluated by
// Laplace at the conditional mode with optional Tierney-Kadane 1/n NLO
// correction.
//
// PD-feasible interval for K_ij given K_{-ij}: D(x) = 1 + 2 c1 x - c2 x^2 > 0
// is a downward-opening quadratic (c2 > 0 since Sigma_0 is PD) with roots
//   x_pm = (c1 +/- sqrt(c1^2 + c2)) / c2,    x_- < 0 < x_+.
// Discriminant c1^2 + c2 = Sigma_0[i,i] Sigma_0[j,j] > 0, so roots are
// always real. K is PD iff K_ij in (x_-, x_+); cost O(1) given c1, c2.
//
// PD-truncation correction (apply_pd_truncation = true):
//   The LO Laplace Z_1D = exp(ell(x*)) sqrt(2 pi / kappa) assumes integration
//   over R, but the encompassing density is supported only on (x_-, x_+).
//   Multiplying by Z_trunc = P(N(x*, 1/kappa) in (x_-, x_+)) tightens the
//   approximation to Z_1D_truncated = Z_1D * Z_trunc, so
//     log pi_enc(K_ij = 0 | K_{-ij}, Y) = -[log Z_1D_LO + delta_NLO + log Z_trunc].
//   In the operational regime 1/sqrt(kappa) << x_+ - x_- this is below the
//   NLO residual; in stress cells (small q, large sigma) it is structural.
//   The flag is opt-in to preserve bit-for-bit parity with Z's pre-2026-05-22
//   primitive for the unit-test suite.
//
// Status codes:
//   0  ok
//   1  K_0 not positive definite
//   2  curvature non-positive (Laplace invalid at mode)
//   3  Newton did not converge within newton_max_iter
//
// Port of ~/SV/Z/R/src/sd_density_at_zero.cpp (companion-AI handoff
// 2026-05-22) + PD-truncation extension (handoff thread 2026-05-22).

namespace ggm_sd {

struct SDResult {
    double log_density;
    double x_mode;
    double curvature;
    double x_minus;       // PD-interval lower endpoint
    double x_plus;        // PD-interval upper endpoint
    double log_Z_trunc;   // log P(N(x_mode, 1/sqrt(curvature)) in (x_-, x_+))
    int    status;
};

SDResult density_at_zero_one(
    const arma::mat& K, int i, int j,
    const arma::mat& S, int n_obs,
    double delta, double sigma,
    bool nlo = true,
    bool apply_pd_truncation = false,
    int newton_max_iter = 50,
    double newton_tol = 1e-10);

}  // namespace ggm_sd
