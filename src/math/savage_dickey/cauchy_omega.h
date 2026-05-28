#pragma once

#include <RcppArmadillo.h>
#include "rng/rng_utils.h"

/**
 * @file cauchy_omega.h
 * @brief Per-edge auxiliary update for the Cauchy slab via scale mixture.
 *
 * The Cauchy slab on an off-diagonal precision entry,
 *     K_ij | gamma = 1  ~  Cauchy(0, sigma),
 * is represented via the scale-mixture-of-normals
 *     omega   ~  InvGamma(1/2, 1/2),
 *     K_ij | gamma=1, omega  ~  N(0, sigma^2 * omega).
 * Conditioning on the L-space framework that the Savage-Dickey
 * between-step uses (rest_L = (l_ii, c_3, K_other, ...) -- K_jj varies
 * with K_ij through Roverato slaving), the effective prior on K_ij
 * picks up a diag-prior contribution. In K-space the (slab x diag)
 * kernel is
 *     -A_K(omega) * K^2 + B_K * K,
 *     A_K(omega) = 1 / (2 sigma^2 omega) + A_diag_K,
 *     B_K        = (model-specific, omega-independent).
 * The conditional p(omega | K, gamma=1, rest_L) therefore depends on
 * omega both directly (through the slab quadratic) and through the
 * omega-dependent normaliser Z(omega) = integral of exp(-A_K(omega)
 * K^2 + B_K K) dK.  The pure-slab conjugate update
 * InvGamma(1, 1/2 + K^2/(2 sigma^2)) is bias-free only when A_diag_K
 * = B_K = 0; under the L-space framework with a nonzero diag prior
 * it is off by ~1% (see experiments/cauchy-slab/06-*).
 *
 * Slice sampling on u = log(omega) handles the omega-dependent
 * normaliser uniformly. The conditional is unimodal and smooth across
 * the realistic chain regime so stepping-out + shrinkage with a
 * default step width of one log-unit is robust.
 *
 * See the bgms manuscript / Savage-Dickey article for the derivation.
 *
 * Model-agnostic: each model that uses this primitive computes A_diag_K
 * and B_K from its own state. For the GGM with Roverato slaving:
 *     A_diag_K =  beta / (2 * l_ii^2),
 *     B_K      = -beta * m_ij / l_ii.
 * The Mixed MRF and other Cholesky-parameterised models would derive
 * analogous expressions from their own L-space framework.
 */

namespace savage_dickey {

/**
 * Slice-sample omega from p(omega | K, gamma=1, rest_L) under the Cauchy
 * scale-mixture representation in the L-space framework.
 *
 * @param K          K-space value of the off-diagonal precision entry.
 * @param sigma2     Nominal slab variance (= sigma^2, NOT sigma^2 * omega).
 * @param A_diag_K   Model-specific diag-prior contribution to the K-space
 *                   quadratic coefficient (e.g., beta/(2 l_ii^2) for GGM).
 * @param B_K        Model-specific diag-prior contribution to the K-space
 *                   linear coefficient (e.g., -beta m_ij / l_ii for GGM).
 * @param omega_curr Current value of omega; slice sampler initial state.
 * @param rng        Random number generator.
 * @return New omega draw; returns omega_curr on slice-sampler failure
 *         (max-expand or max-shrink reached without acceptance).
 */
double slice_sample_cauchy_omega_active(
    double K,
    double sigma2,
    double A_diag_K,
    double B_K,
    double omega_curr,
    SafeRNG& rng);

/**
 * Sample omega from its prior, InvGamma(1/2, 1/2). Used at gamma = 0
 * where K_ij = 0 and the slab carries no information about omega.
 *
 * @param rng Random number generator.
 * @return Draw from InvGamma(1/2, 1/2).
 */
double sample_cauchy_omega_prior(SafeRNG& rng);

}  // namespace savage_dickey
