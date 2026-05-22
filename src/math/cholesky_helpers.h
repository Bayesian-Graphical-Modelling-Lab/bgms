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

/**
 * Trailing 2×2 block of the Cholesky factor under symmetric permutation that
 * places (i, j) at the last two positions.
 *
 * Given upper-triangular U with K = UᵀU and indices i < j < p, computes the
 * permuted lower-triangular factor L_perm satisfying L_perm L_permᵀ = P K Pᵀ
 * where P sends i → p−2 and j → p−1, and returns the four entries of the
 * trailing 2×2 (l_ii, l_ji, l_jj).
 *
 * Implementation: Bunch-style adjacent symmetric swaps with Givens rotations.
 * Each swap is O(p) and unconditionally stable; total cost O(p²) per call.
 * Does NOT go through Σ_BB = K⁻¹_BB, so avoids the 1/det amplification that
 * dominates error near the PD boundary.
 *
 * @param U   Upper-triangular Cholesky factor of K (K = UᵀU).
 * @param i   First edge index (0-based), with i < j.
 * @param j   Second edge index (0-based), with j < p.
 * @return    Struct with l_ii, l_ji, l_jj and ok (false if a non-positive
 *            diagonal is produced, i.e. K is numerically singular at the edge).
 */
struct Perm2x2Result {
    double l_ii;
    double l_ji;
    double l_jj;
    bool   ok;
};

Perm2x2Result perm_to_trailing_2x2(const arma::mat& U, size_t i, size_t j,
                                   arma::mat& L_scratch);

// Convenience overload that allocates its own scratch. Use in tests / one-off
// callers; hot loops should pass a pre-sized scratch to avoid allocation.
inline Perm2x2Result perm_to_trailing_2x2(const arma::mat& U, size_t i,
                                          size_t j) {
    arma::mat L_scratch;
    return perm_to_trailing_2x2(U, i, j, L_scratch);
}

} // namespace cholesky_helpers
