# Blockwise / Subspace HMC for Sparse GGMs

## Concept
Split the precision matrix parameters into clique-based subsets that admit local Cholesky factorizations without fill-ins, then run Hamiltonian proposals (HMC or NUTS) within each subset while conditioning on the rest. The approach mimics block Gibbs but uses gradient-informed trajectories to accelerate mixing inside large cliques.

## Construction
1. Decompose the graph into cliques and separators (junction tree or block-cut tree).
2. For each clique, identify the free precision entries (including its separator interface).
3. Parameterize each block via a local Cholesky factor (log-diagonals + unconstrained subdiagonal entries). Because the block is chordal by definition, PD is automatic.
4. Run NUTS inside the block, treating separator values from neighboring blocks as fixed during the leapfrog trajectory.
5. Accept/reject per block, then cycle through blocks or randomly scan them.

## Pros
- Keeps parameter dimensionality close to the true number of edges.
- Exploits sparsity: small cliques reuse cached symbolic factorizations.
- Amenable to parallel updates across disconnected components.

## Cons
- Requires dynamic bookkeeping of separators when RJ moves add/remove edges.
- Acceptance rates can drop if separators induce strong cross-block dependence.
- Implementation complexity: need reversible Metropolis-within-Gibbs logic and gradient code per block.

## References
- Green, P. J., & Thomas, A. (2013). *Sampling decomposable graphs using a Markov chain on junction trees.* Biometrika 100(1).
- Mohammadi, A., & Wit, E. C. (2015). *Bayesian structure learning in sparse Gaussian graphical models.* Bayesian Analysis 10(1).
