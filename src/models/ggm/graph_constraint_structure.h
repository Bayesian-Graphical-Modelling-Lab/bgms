#pragma once

#include <RcppArmadillo.h>
#include <vector>

/**
 * Per-column constraint information for the free-element Cholesky
 * parameterization.
 *
 * For column q of the Cholesky factor Phi, the off-diagonal entries
 * x_q = Phi[1:q-1, q] must satisfy A_q x_q = 0 for excluded edges.
 * The null-space basis N_q has d_q = (q-1) - m_q columns, where m_q
 * is the number of excluded edges in column q.
 */
struct ColumnConstraints {
    /// 0-based row indices of excluded edges (i where edge (i,q) is absent).
    std::vector<size_t> excluded_indices;
    /// 0-based row indices of included edges (i where edge (i,q) is present).
    std::vector<size_t> included_indices;
    /// Number of excluded edges: m_q = |excluded_indices|.
    size_t m_q;
    /// Number of free parameters (included edges): d_q = (q-1) - m_q.
    size_t d_q;
};

/**
 * Precomputed constraint structure for a fixed graph.
 *
 * Built once when the graph changes (edge indicator toggle). Stores
 * the shape of constraints for each column — which rows are excluded
 * and included — but not the constraint matrix values themselves
 * (those depend on Phi and change during leapfrog integration).
 *
 * The theta vector layout is:
 *   (psi_1, [f_2, psi_2], [f_3, psi_3], ..., [f_p, psi_p])
 * where f_q has d_q elements and psi_q is the log-diagonal.
 */
struct GraphConstraintStructure {
    /// Per-column constraint information, indexed by q (0-based).
    std::vector<ColumnConstraints> columns;
    /// Active parameter dimension: p + sum(d_q) = p + |E|.
    size_t active_dim;
    /// Full parameter dimension: p + p(p-1)/2 (all possible off-diag slots).
    size_t full_dim;
    /// Offset of each column's block in the theta vector.
    /// theta_offsets[q] is where (f_q, psi_q) starts in theta.
    std::vector<size_t> theta_offsets;
    /// Offset of each column's block in the full (zero-padded) theta vector.
    /// full_theta_offsets[q] is where the (q-1) off-diag + 1 diag slots start.
    std::vector<size_t> full_theta_offsets;
    /// Number of variables.
    size_t p;

    /**
     * Build the constraint structure from an edge indicator matrix.
     *
     * @param edge_indicators  p x p symmetric 0/1 matrix
     */
    void build(const arma::imat& edge_indicators) {
        p = edge_indicators.n_rows;
        columns.resize(p);

        active_dim = 0;
        full_dim = 0;
        theta_offsets.resize(p);
        full_theta_offsets.resize(p);

        for (size_t q = 0; q < p; ++q) {
            auto& col = columns[q];
            col.excluded_indices.clear();
            col.included_indices.clear();

            for (size_t i = 0; i < q; ++i) {
                if (edge_indicators(i, q) == 1) {
                    col.included_indices.push_back(i);
                } else {
                    col.excluded_indices.push_back(i);
                }
            }

            col.m_q = col.excluded_indices.size();
            col.d_q = col.included_indices.size();

            theta_offsets[q] = active_dim;
            full_theta_offsets[q] = full_dim;

            // f_q (d_q entries) + psi_q (1 entry)
            active_dim += col.d_q + 1;
            // (q-1) slots for off-diag + 1 for diag in full vector
            full_dim += q + 1;
        }
    }

    /**
     * Offset of psi_q within the active theta vector.
     * psi_q is the last element of column q's block.
     */
    size_t psi_offset(size_t q) const {
        return theta_offsets[q] + columns[q].d_q;
    }

    /**
     * Offset of psi_q within the full (zero-padded) theta vector.
     * In the full vector, column q has (q-1) off-diag slots + 1 diag slot.
     * psi_q is at position full_theta_offsets[q] + q (the last slot).
     */
    size_t full_psi_offset(size_t q) const {
        return full_theta_offsets[q] + q;
    }
};
