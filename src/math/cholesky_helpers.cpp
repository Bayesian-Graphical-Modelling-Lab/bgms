/**
 * @file cholesky_helpers.cpp
 * @brief Non-inline cholesky_helpers definitions.
 */

#include "cholesky_helpers.h"

#include <cmath>

namespace cholesky_helpers {

namespace {

// Symmetric adjacent swap of rows/cols (k, k+1) of K = L Lᵀ, maintaining L
// lower-triangular. Implements P L Pᵀ followed by a single Givens rotation on
// columns (k, k+1) that zeros the off-diagonal entry created by the
// permutation. Cost O(p − k).
//
// Returns false if the diagonal entry L(k, k) at the (k, k+1) position before
// rotation is non-positive (singular K). In that case L is left in an
// inconsistent state.
inline bool swap_rows_cols_adjacent(arma::mat& L, std::size_t k) {
    const std::size_t p = L.n_rows;
    L.swap_rows(k, k + 1);
    L.swap_cols(k, k + 1);

    // Row k now has a non-zero at column k+1 (the "violating" entry).
    // Apply a Givens rotation on columns (k, k+1) to zero L(k, k+1).
    const double a = L(k, k);
    const double b = L(k, k + 1);
    const double r = std::hypot(a, b);
    if (!(r > 0.0)) {
        return false;
    }
    const double c = a / r;
    const double s = b / r;

    // Apply the Givens rotation to columns (k, k+1) at rows k .. p-1.
    // Rows < k have zero in both columns (lower-triangular structure).
    L(k, k)     = r;
    L(k, k + 1) = 0.0;
    for (std::size_t row = k + 1; row < p; ++row) {
        const double a_row = L(row, k);
        const double b_row = L(row, k + 1);
        L(row, k)     =  c * a_row + s * b_row;
        L(row, k + 1) = -s * a_row + c * b_row;
    }
    return true;
}

}  // namespace

Perm2x2Result perm_to_trailing_2x2(const arma::mat& U, std::size_t i,
                                   std::size_t j, arma::mat& L_scratch) {
    Perm2x2Result out{0.0, 0.0, 0.0, false};
    const std::size_t p = U.n_rows;
    if (!(i < j && j < p && p >= 2)) {
        return out;
    }

    // Convert upper-triangular U → lower-triangular L (K = L Lᵀ = Uᵀ U). Use
    // the caller's scratch buffer to avoid per-call allocation; arma's
    // assignment-to-same-shape reuses the existing memory.
    if (L_scratch.n_rows != p || L_scratch.n_cols != p) {
        L_scratch.set_size(p, p);
    }
    L_scratch = U.t();
    arma::mat& L = L_scratch;

    // Move j → p-1 (via adjacent symmetric swaps). Since j > i, these swaps do
    // not touch position i.
    for (std::size_t k = j; k + 1 < p; ++k) {
        if (!swap_rows_cols_adjacent(L, k)) {
            return out;
        }
    }
    // Move i → p-2 (adjacent swaps). Position p-1 (holding original j) is not
    // touched by these swaps.
    for (std::size_t k = i; k + 2 < p; ++k) {
        if (!swap_rows_cols_adjacent(L, k)) {
            return out;
        }
    }

    out.l_ii = L(p - 2, p - 2);
    out.l_ji = L(p - 1, p - 2);
    out.l_jj = L(p - 1, p - 1);
    out.ok   = (out.l_ii > 0.0) && (out.l_jj > 0.0);
    return out;
}

}  // namespace cholesky_helpers
