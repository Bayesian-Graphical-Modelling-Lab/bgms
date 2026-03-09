# Chordal Completion Augmentation for GGM NUTS

## Goal
Use the No-U-Turn Sampler for Gaussian graphical model parameters while enforcing the sparsity pattern of a non-chordal precision matrix. The idea is to work on an augmented parameter space where positive definiteness is guaranteed, then map back to the original graph.

## Strategy Overview
1. **Chordal completion**: Given the current graph \(G\), build a chordal supergraph \(\tilde G\) by adding the minimal set of fill-in edges. Algorithms such as minimum degree, approximate minimum degree, or nested dissection supply reasonable triangulations.
2. **Perfect elimination order**: Compute a PEO on \(\tilde G\). This ordering ensures that sparse Cholesky factorization of any PD matrix respecting \(\tilde G\) introduces no further fill-ins.
3. **Sparse Cholesky parameterization**:
   - Parameterize a lower-triangular matrix \(L\) whose sparsity matches \(\tilde G\). Entries corresponding to edges (including fill-ins) are free parameters; non-edges are fixed at zero.
   - Diagonal entries live on the log scale to enforce positivity: \(L_{ii} = \exp(d_i)\) with unconstrained \(d_i \in \mathbb{R}\).
   - All other free entries are unconstrained reals. Together they form the vector that NUTS updates.
4. **Precision reconstruction**: Form \(Q = L L^\top\). By construction \(Q\) is positive definite for any unconstrained parameter values.
5. **Projection back to \(G\)**: Extract the entries of \(Q\) that correspond to edges and diagonals of the original graph, discarding the auxiliary fill-in elements. These extracted entries plug into the likelihood and priors defined on \(G\).

## Why This Works
- Cholesky factors of non-chordal sparse matrices inevitably introduce fill-ins. By expanding to \(\tilde G\) we accept those fill-ins explicitly, rather than letting them appear implicitly during factorization.
- NUTS operates on the free parameters of \(L\), which live on \(\mathbb{R}^{k}\) with no boundary constraints, so the Hamiltonian dynamics remain well-defined.
- The auxiliary parameters tied to fill-in edges never enter the likelihood directly; they ensure the reconstructed precision stays PD while respecting the observed sparsity once projected.

## Computational Cost
- Building \(\tilde G\) and its PEO costs roughly linear time in the number of edges plus fill-ins but can be dominated by triangulation for dense graphs.
- The dimension of the parameter space increases with the number of fill-in edges. If \(G\) is far from chordal, this inflation can slow NUTS steps and gradient evaluations.
- Sparse matrix routines (symbolic factorization, reuse of elimination trees) should be cached across iterations to keep costs manageable.

## Practical Tips
- Reuse the same chordal completion for all iterations while the graph structure stays fixed. If the edge set changes (RJ step), recompute the triangulation only when necessary.
- Choose heuristics (e.g., AMD) that minimize fill-ins, balancing numerical stability and efficiency.
- Gradients of \(Q = L L^\top\) w.r.t. the free entries of \(L\) follow standard reverse-mode rules; implement them with sparse structures to avoid dense \(O(p^2)\) work.

## References
- Rue, H., & Held, L. (2005). *Gaussian Markov Random Fields: Theory and Applications*. Chapters 2–3 discuss chordal completions and sparse Cholesky factors.
- Atay-Kayis, A., & Massam, H. (2005). "A Monte Carlo method for computing the marginal likelihood in nondecomposable Gaussian graphical models." *Journal of Computational and Graphical Statistics*.
- Dahl, J., Vandenberghe, L., & Roychowdhury, V. (2008). "Covariance selection for non-chordal graphs via chordal embedding." *Optimization and Engineering*.
