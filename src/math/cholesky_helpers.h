#pragma once

/**
 * @file cholesky_helpers.h
 * @brief Shared algebraic helpers for Cholesky-based precision updates.
 *
 * Pure functions with no model-specific state.  Used by both GGMModel and
 * MixedMRFModel for proposal constant extraction and log-determinant
 * computation.
 */

#include <RcppArmadillo.h>
#include <cmath>

namespace cholesky_helpers {

/**
 * Log-determinant of a positive-definite matrix from its upper-triangular
 * Cholesky factor R (where Ω = R'R).
 *
 * @param R  Upper-triangular Cholesky factor.
 * @return   log|Ω| = 2 Σ log(R_ii).
 */
inline double get_log_det(const arma::mat& R) {
    return 2.0 * arma::accu(arma::log(R.diag()));
}

/**
 * Schur complement element: A(ii,jj) − A(ii,i) A(jj,i) / A(i,i).
 *
 * Used to compute entries of the inverse of a submatrix from the full
 * covariance matrix.
 *
 * @param A   Symmetric positive-definite matrix.
 * @param i   Conditioning index.
 * @param ii  Row index of the desired element.
 * @param jj  Column index of the desired element.
 * @return    Schur complement entry.
 */
inline double compute_inv_submatrix_i(const arma::mat& A, size_t i,
                                      size_t ii, size_t jj) {
    return A(ii, jj) - A(ii, i) * A(jj, i) / A(i, i);
}

} // namespace cholesky_helpers
