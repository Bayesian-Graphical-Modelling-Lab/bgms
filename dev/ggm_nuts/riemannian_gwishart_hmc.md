# Riemannian Manifold HMC on the G-Wishart Cone

## Concept
Treat the space of graph-constrained positive-definite matrices as a Riemannian manifold endowed with the Fisher information metric of the G-Wishart distribution. Run Riemannian HMC (RMHMC) so that geodesic proposals respect the sparsity constraints without auxiliary completions.

## Construction
1. Define the target density as the posterior over the sparse precision matrix, typically proportional to a G-Wishart prior times the data likelihood.
2. Use the Fisher metric of the G-Wishart (or an approximation) as the position-dependent mass matrix.
3. At each leapfrog step solve the implicit equations for RMHMC (usually via fixed-point iterations) to update momentum and position.
4. Ensure gradients are evaluated only on free entries; metric tensors inherit the sparse structure, but linear solves involve fill-ins.

## Pros
- Moves directly on the constrained space; no auxiliary parameters or projections.
- Metric adapts to local curvature, improving exploration in stiff posteriors.
- Naturally accommodates different priors via the metric definition.

## Cons
- Each leapfrog step requires solving nonlinear equations plus sparse linear systems; computationally heavy.
- Implementing efficient sparse RMHMC is complex (automatic differentiation of metrics, symbolic factorization reuse).
- Requires tuning step sizes and number of fixed-point iterations to maintain stability.

## References
- Byrne, S., & Girolami, M. (2013). *Geodesic Monte Carlo on embedded manifolds.* Scandinavian Journal of Statistics 40(4).
- Lenkoski, A. (2013). *A direct sampler for G-Wishart variates.* Stat 2(1). (Appendix discusses the manifold structure.)
