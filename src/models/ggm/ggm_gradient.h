#pragma once

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <utility>
#include "models/ggm/graph_constraint_structure.h"

/**
 * Forward map result: Phi, K, and cached QR data for the backward pass.
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
    /// Per-column psi values (log-diagonal of Phi).
    arma::vec psi;
};

/**
 * GGM gradient engine for the free-element Cholesky parameterization.
 *
 * Computes the forward map theta -> (Phi, K) with Jacobian and the
 * reverse-mode gradient of the log-posterior. The cross-column
 * adjoint (backward pass through the null-space basis N_q) uses
 * finite differences, matching the validated R prototype.
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
     * computes QR of A_q^T for the null-space basis N_q, sets
     * x_q = N_q f_q, and accumulates the Jacobian.
     *
     * @param theta  Parameter vector of length p + |E|
     * @return ForwardMapResult with Phi, K, log|det J|, and cached QR data
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

    static void compute_null_basis_and_rdiag(
        const arma::mat& Aq,
        arma::mat& Nq,
        arma::vec& R_diag);

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

    // --- Cached references (set by rebuild) ---
};
