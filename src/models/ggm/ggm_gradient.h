#pragma once

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <utility>
#include "models/ggm/graph_constraint_structure.h"

/**
 * One stored Givens rotation: rows (r1, r2), angle (c, s), column j.
 */
struct GivensRotation {
    double c, s;
    size_t r1, r2;
    size_t col;   ///< Column being zeroed (needed for backward pass step E).
};

/**
 * Forward map result: Phi, K, and cached Givens data for the backward pass.
 *
 * The free-element Cholesky parameterization maps an unconstrained
 * theta vector to a positive-definite precision matrix K = Phi^T Phi,
 * where Phi is upper-triangular and excluded edges are enforced via
 * null-space constraints.
 */
struct ForwardMapResult {
    /// Upper-triangular Cholesky factor (p x p).
    arma::mat Phi;
    /// Precision matrix K = Phi^T Phi (p x p).
    arma::mat K;
    /// Log-determinant of the Jacobian |det J|.
    double log_det_jacobian;

    /// Per-column null-space basis N_q (q-1 x d_q). Empty for q=0 or d_q=0.
    std::vector<arma::mat> Nq;
    /// Per-column QR R-factor diagonals (for Jacobian). Empty when m_q=0.
    std::vector<arma::vec> R_diag;
    /// Per-column stored Givens rotations from the QR of A_q^T.
    /// Used by the reverse-Givens backward pass. Empty when m_q=0.
    std::vector<std::vector<GivensRotation>> givens_rotations;
    /// Per-column Q matrix from Givens QR of A_q^T.
    /// Needed by the reverse-Givens backward pass for c_bar/s_bar.
    std::vector<arma::mat> Q_full;
    /// Per-column R matrix (working matrix after Givens QR). n x m_q.
    std::vector<arma::mat> R_full;
    /// Per-column psi values (log-diagonal of Phi).
    arma::vec psi;
};

/**
 * GGM gradient engine for the free-element Cholesky parameterization.
 *
 * Computes the forward map theta -> (Phi, K) with Jacobian and the
 * reverse-mode gradient of the log-posterior. The cross-column
 * adjoint (backward pass through the null-space basis N_q) uses
 * reverse-mode differentiation through stored Givens rotations,
 * giving an exact analytic gradient for all constraint dimensions.
 *
 * To avoid per-leapfrog-step allocation, all workspace matrices are
 * pre-allocated and reused. Call rebuild() when the graph changes.
 */
class GGMGradientEngine {
public:
    GGMGradientEngine() = default;

    /**
     * Initialize the engine for a given graph structure.
     *
     * Pre-allocates all workspace matrices sized to the graph.
     * Must be called before forward_map() or logp_and_gradient().
     *
     * @param structure   Precomputed graph constraint structure
     * @param n           Sample size
     * @param suf_stat    Sufficient statistic S = X^T X (p x p)
     * @param pairwise_scale  Cauchy slab scale for off-diagonal K entries
     */
    void rebuild(const GraphConstraintStructure& structure,
                 size_t n,
                 const arma::mat& suf_stat,
                 double pairwise_scale);

    /**
     * Forward map: theta -> (Phi, K, log|det J|).
     *
     * Processes columns left-to-right, building Phi column by column.
     * For each column q >= 2: builds A_q from earlier Phi columns,
     * computes Givens QR of A_q^T for the null-space basis N_q, sets
     * x_q = N_q f_q, and accumulates the Jacobian.
     *
     * @param theta  Parameter vector of length p + |E|
     * @return ForwardMapResult with Phi, K, log|det J|, and cached Givens data
     */
    ForwardMapResult forward_map(const arma::vec& theta) const;

    /**
     * Combined log-posterior and gradient evaluation.
     *
     * Runs the forward map, then the backward pass to compute the
     * gradient of:
     *   L = (n/2) log|K| - (1/2) tr(KS) + log_prior(K) + log|det J|
     *
     * @param theta  Parameter vector of length p + |E|
     * @return (log-posterior value, gradient vector)
     */
    std::pair<double, arma::vec> logp_and_gradient(const arma::vec& theta) const;

    /**
     * Full-space log-posterior and gradient for RATTLE integration.
     *
     * Operates on the full position vector x in R^{p(p+1)/2} —
     * the raw Cholesky entries (off-diagonal) and log-diagonal (psi).
     * No null-space transformation or QR decomposition is needed:
     * the gradient w.r.t. each x_{iq} is simply Phi_bar_{iq}.
     *
     * The Jacobian includes only the Cholesky-to-K and log-diagonal
     * terms. The QR Jacobian terms from the null-space basis are not
     * needed because RATTLE preserves the correct measure on the
     * constraint manifold.
     *
     * @param x  Full position vector of dimension p(p+1)/2
     * @return (log-posterior value, gradient vector of same dimension)
     */
    std::pair<double, arma::vec> logp_and_gradient_full(const arma::vec& x) const;

    /**
     * Givens QR of an n x m matrix M (n >= m).
     *
     * Computes M = Q R via bottom-to-top Givens rotations. Stores the
     * rotation sequence for reverse-mode differentiation.
     *
     * @param M       Input matrix (n x m)
     * @param Q       Output: orthogonal matrix (n x n)
     * @param R       Output: upper-trapezoidal matrix (n x m)
     * @param R_diag  Output: absolute diagonal of R (length min(n,m))
     * @param rots    Output: stored Givens rotations
     */
    static void givens_qr(
        const arma::mat& M,
        arma::mat& Q,
        arma::mat& R,
        arma::vec& R_diag,
        std::vector<GivensRotation>& rots);

    static void build_Aq(const arma::mat& Phi,
                         const ColumnConstraints& col,
                         size_t q,
                         arma::mat& Aq);

private:
    const GraphConstraintStructure* structure_ = nullptr;
    size_t n_ = 0;
    size_t p_ = 0;
    const arma::mat* suf_stat_ = nullptr;
    double pairwise_scale_ = 1.0;
};
