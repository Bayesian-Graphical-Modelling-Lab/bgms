# Projected NUTS for Sparse Precision Matrices

## Concept
Run Hamiltonian dynamics in an unconstrained Euclidean space, then project intermediate states back onto the sparse positive-definite cone that encodes the graph. Projection enforces symmetry, sparsity, and PD by solving a convex program after each integrator step.

## Construction
1. Parameterize a dense symmetric matrix \(S\) with unconstrained entries.
2. During each leapfrog step, update \(S\) using standard NUTS integration.
3. After the position update, solve a projection problem:
   \[
   \min_{Q \succeq 0} \|Q - S\|_F^2 \quad \text{s.t.}\quad Q_{ij} = 0 \;\text{for}\; (i,j) \notin E.
   \]
   This is a log-det barrier or semidefinite program over the sparsity-constrained cone.
4. Use the projected \(Q\) to evaluate the log posterior and gradients (via implicit differentiation through the projection).

## Pros
- Parameter dimension equals the number of edges (after eliminating fixed zeros) if you operate in the sparse coordinate basis.
- Projection guarantees feasibility even if leapfrog steps wander outside the cone.
- Compatible with existing convex-optimization libraries (e.g., sparse log-det solvers).

## Cons
- Each projection may require solving a sizable convex problem (costly for large graphs).
- Differentiating through the projection demands custom adjoint code or implicit differentiation.
- Numerical errors from repeated projection can accumulate, making reversibility checks essential.

## References
- Dahl, J., Vandenberghe, L., & Roychowdhury, V. (2008). *Covariance selection for non-chordal graphs via chordal embedding.* Optimization and Engineering.
- Nishihara, R., et al. (2014). *Proximal algorithms for constrained HMC.* (Workshop paper; describes projection-based integrators.)
