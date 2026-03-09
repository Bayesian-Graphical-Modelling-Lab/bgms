# Elliptical Slice / Surrogate Transitions

## Concept
Retain the graph-respecting adaptive Metropolis kernel but wrap it with gradient-informed or Gaussian surrogate moves to boost local exploration without heavy Hamiltonian mechanics. Elliptical slice sampling (ESS) is attractive when the prior is Gaussian with sparse precision; surrogate transitions such as preconditioned Crank–Nicolson (pCN) or blocked ESS reuse sparse linear solves.

## Construction
1. Use the existing RJ step to update the edge set.
2. Conditional on a fixed graph, draw proposals from a Gaussian surrogate with precision equal (or close) to the current posterior curvature, e.g., via solving \(Q x = z\) with sparse methods.
3. Apply ESS or pCN to accept/reject along the elliptical trajectory without tuning step sizes.
4. Optionally interleave occasional gradient-informed Metropolis proposals (Langevin, quasi-Newton) that exploit the sparse precision as a preconditioner.

## Pros
- No need for chordal completions or manifold integrators.
- Exact ESS has acceptance probability 1 when the likelihood is Gaussian; in general it is robust to tuning.
- Sparse linear algebra can be reused across iterations, matching the structure already present for RJ moves.

## Cons
- Mixing can lag behind full NUTS, especially with highly correlated posteriors.
- Requires drawing from the surrogate Gaussian, which still needs sparse Cholesky solves each iteration.
- Does not exploit gradients beyond what is embedded in the surrogate; tuning the surrogate to match posterior curvature can be difficult.

## References
- Murray, I., Adams, R. P., & MacKay, D. J. C. (2010). *Elliptical slice sampling.* AISTATS.
- Cotter, S., et al. (2013). *MCMC methods for functions: modifying old algorithms to make them faster.* Statistical Science 28(3).
