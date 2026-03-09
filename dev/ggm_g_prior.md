# g-type Slab Prior for bgms Edge Weights

## Motivation
Manual tuning of the Cauchy slab scale is awkward and can obscure whether discrepancies with comparators stem from modeling or shrinkage choices. In regression-style spike-and-slab models, Zellner's $g$-prior offers automatic scaling by tying coefficients to the design matrix. We can mimic that idea for individual edge parameters $\theta_{ij}$ in the mixed GGM sampler by anchoring the slab variance to local Fisher information, giving consistent units across nodes and datasets.

## Proposed Formulation
For each candidate edge $(i,j)$ with corresponding conditional regression coefficient $\theta_{ij}$:

- Use the usual spike mixture: $\theta_{ij} \sim (1-\gamma_{ij}) \delta_0 + \gamma_{ij} \cdot p_{\text{slab}}(\theta_{ij})$ with inclusion indicator $\gamma_{ij} \in \{0,1\}$.
- Replace the free-scale Cauchy slab with a conditionally normal slab whose variance adapts via a $g$ parameter:
  $$p_{\text{slab}}(\theta_{ij} \mid g, s_{ij}) = \mathcal N\bigl(0, g \cdot s_{ij}^{-1}\bigr).$$
- Here $s_{ij}$ approximates the Fisher information for $\theta_{ij}$ when the rest of the graph is fixed. For Gaussian nodes this could be the conditional residual variance; for mixed nodes we can use the negative Hessian of the log-likelihood restricted to $\theta_{ij}$.

### Hyperprior on $g$
Adopt a standard Beta-prime / inverse-gamma style:
- $g \sim \text{BetaPrime}(a_g, b_g)$ with density proportional to $g^{a_g-1}(1+g)^{-(a_g+b_g)}$.
- Default choice $a_g = 1$, $b_g = 1$ mirrors the hyper-$g$ prior, favoring moderate shrinkage while leaving heavy tails.
- Alternative: $g/(1+g) \sim \text{Beta}(a_g, b_g)$ to keep $g$ finite with more intuitive tuning.

## Implementation Sketch
1. **Compute $s_{ij}$ efficiently**: during local edge updates we already build score/Hessian terms. Cache $s_{ij}$ (or its approximation) when evaluating the conditional likelihood.
2. **Augment sampler state**: treat $g$ as either global (shared across edges) or edge-specific $g_{ij}$. Start with a global $g$ for simplicity.
3. **Gibbs/Metropolis steps**:
   - Conditional on $g$, sampling $\theta_{ij}$ from the Gaussian slab is straightforward (or use MH if the mixed-node likelihood prevents conjugacy).
   - Update $g$ via MH or slice sampling using the Beta-prime prior and the product of Gaussian slabs.
4. **Compatibility with spike indicator**: no change—the inclusion probability update uses the new slab density.

## Evaluation Plan
- **Full-graph scenario**: compare bgms with Cauchy slabs vs g-type slabs on $p=20$ dense graphs; track RMSE and coverage relative to truth and bggm.
- **Sparse scenario**: measure inclusion ROC/AUC and calibration when using g-type slabs, ensuring sparsity is not overly penalized.
- **Sensitivity**: vary $(a_g, b_g)$ to see if results are robust; consider placing a hyperprior on $s_{ij}$ approximation error if needed.

## Open Questions
- How stable is $s_{ij}$ for mixed discrete nodes? Need to confirm the Hessian remains positive definite; otherwise use a sandwich estimator.
- Global vs local $g$: does sharing one $g$ across edges over-shrink high-degree nodes? Potential extension with $g_i$ per node.
- Interaction with existing spike penalties: should inclusion priors be adjusted when switching slab families to keep edge probabilities comparable?

## Next Steps
1. Prototype $s_{ij}$ computation in the continuous-only case to validate scaling.
2. Implement a global-$g$ sampler branch and benchmark on small graphs.
3. Decide whether to expose the g-type slab as a user option or keep it experimental during the performance study.
