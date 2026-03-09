# Partial-correlation reparameterization for GGM NUTS

## Core idea
Map unconstrained parameters to admissible precision matrices by working in partial correlations (PCs). Each PC lives on \((-1, 1)\), so we represent it via \(\rho = \tanh(\eta)\) with \(\eta \in \mathbb{R}\). Diagonal precision entries use log-scales to ensure positivity. The mapping guarantees positive definiteness and enforces graph zeros, allowing NUTS to operate on \(\eta\) and the log-diagonals without boundary constraints.

## Construction steps
1. Choose a node ordering and form sequential regressions (modified Cholesky / Schur complements).
2. For every permitted edge \((i, j)\) with \(j < i\), introduce an unconstrained \(\eta_{ij}\); set the partial correlation \(\rho_{ij} = \tanh(\eta_{ij})\).
3. Assemble the unit-diagonal correlation matrix using the recursive PC-to-correlation formulas (Joe, 2006; Liechty et al., 2004). Missing edges receive \(\rho = 0\).
4. Introduce unconstrained \(\zeta_i\) for each node and set conditional precision scales as \(\exp(\zeta_i)\).
5. Form \(Q = D^{-1/2} R^{-1} D^{-1/2}\) where \(D = \operatorname{diag}(\exp(\zeta_i))\) and \(R\) is the correlation matrix from step 3.

## Graph enforcement
Only existing edges get \(\eta\)-parameters, so their PCs can vary; missing edges stay fixed at zero PCs, and the recursive construction keeps those entries zero in the final precision. Thus \(Q\) respects the graph after the similarity transform.

## Using with NUTS
- Parameter vector consists of \(\eta\) (one per edge) plus \(\zeta\) (one per node).
- Log posterior is evaluated through \(Q(\eta, \zeta)\).
- Gradients flow through the smooth tanh and matrix operations, so autodiff applies. Cache intermediate matrices to keep complexity reasonable.

## Practical notes
- Clamp \(|\rho|\) slightly below 1 (e.g., \(\rho = \tanh(\eta)\times (1-\epsilon)\)) for numerical stability.
- Reorder nodes (AMD, nested dissection) to limit numerical error and keep regressions well-conditioned.
- Jacobian adjustments: add \(\log \cosh(\eta_{ij})\) for each PC transform and \(\zeta_i\) for each log-diagonal.

## References
- Liechty, J., Liechty, M., & Müller, P. (2004). Bayesian correlation estimation. *Journal of the American Statistical Association*, 99(447), 1021–1035.
- Joe, H. (2006). Generating random correlation matrices based on partial correlations. *Journal of Multivariate Analysis*, 97(10), 2177–2189. (See also UBC Technical Report 2005-11.)
