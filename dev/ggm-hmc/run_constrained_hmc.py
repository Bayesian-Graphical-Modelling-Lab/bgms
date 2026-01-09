"""
Constrained HMC sampler for GGM precision matrices.

This script is called from R via subprocess to sample precision matrices
with specified zero constraints using the mici library.

Usage: python run_constrained_hmc.py <input_json_file>
"""

import json
import sys
import time
import os
import numpy as np
import mici
import jax
from numpyro.distributions.transforms import LowerCholeskyTransform
import arviz  # For diagnostics

# Force CPU backend to avoid CUDA autotuner issues with small matrices
# Comment out this line to use GPU (may require CUDA driver updates)
jax.config.update("jax_default_device", jax.devices("cpu")[0])

jax.config.update("jax_enable_x64", True)


def main():
    # Load input data from R
    with open(sys.argv[1], "r") as f:
        data = json.load(f)

    n_variable = data["n_variable"]
    n_obs = data["n_obs"]
    S = np.array(data["scatter_matrix"], dtype=np.float64)
    zero_indices = np.array(data["zero_indices"], dtype=np.int64)
    n_warm_up_iter = data["n_warm_up_iter"]
    n_main_iter = data["n_main_iter"]
    n_chain = data["n_chain"]
    seed = data["seed"]
    output_file = data["output_file"]
    samples_file = data["samples_file"]

    print(f"Loaded data: n_obs={n_obs}, n_variable={n_variable}, n_zero_pairs={len(zero_indices)}")
    print(f"Chains: {n_chain}, Warmup: {n_warm_up_iter}, Samples: {n_main_iter}")
    sys.stdout.flush()

    # Precompute Cholesky factor of scatter matrix for efficient trace computation
    # tr(Omega S) = tr(L L^T S) = ||S_chol^T @ L||_F^2
    S_chol = np.linalg.cholesky(S)

    # Set up transformations
    vector_to_cholesky = LowerCholeskyTransform()
    cholesky_to_vector = vector_to_cholesky.inv

    def constr(u, zero_indices):
        L = vector_to_cholesky(u)
        return jax.vmap(lambda i, j: L[i] @ L[j])(*zero_indices.T)

    def neg_log_dens(u, n_obs, S_chol):
        """
        Negative log posterior for GGM precision matrix.

        log p(Omega | X) ∝ (n/2) * log|Omega| - (1/2) * tr(Omega * S) + log p(Omega)

        where S = X'X is the scatter matrix.
        We use a standard normal prior on the unconstrained parameters u.

        Optimizations:
        - Avoid forming Omega = L @ L.T explicitly
        - Use tr(L L^T S) = ||S_chol^T @ L||_F^2 where S = S_chol @ S_chol^T
        - This reduces from 2 O(p³) ops to 1 O(p³) op
        """
        L = vector_to_cholesky(u)

        # Log determinant: log|Omega| = 2 * sum(log(diag(L)))
        log_det_Omega = 2 * jax.numpy.log(jax.numpy.diag(L)).sum()

        # Trace term: tr(Omega * S) = tr(L L^T S) = ||S_chol^T @ L||_F^2
        # This avoids explicitly forming Omega = L @ L.T
        S_chol_T_L = S_chol.T @ L
        trace_term = jax.numpy.sum(S_chol_T_L ** 2)

        # Gaussian likelihood: (n/2) * log|Omega| - (1/2) * tr(Omega * S)
        log_likelihood = (n_obs / 2) * log_det_Omega - 0.5 * trace_term

        # Prior on unconstrained parameters (standard normal)
        log_prior = -0.5 * (u**2).sum()

        # Return negative log posterior
        return -(log_likelihood + log_prior)

    def trace_func(state):
        L = vector_to_cholesky(state.pos)
        return {"u": state.pos, "P": L @ L.T}

    # Create initial states
    rng = np.random.default_rng(seed)
    scale = 0.01

    init_states = []
    for c in range(n_chain):
        random_matrix = rng.standard_normal((n_variable, n_variable))
        P_init = np.identity(n_variable) + scale * random_matrix @ random_matrix.T
        P_init[zero_indices[:, 0], zero_indices[:, 1]] = 0.0
        P_init[zero_indices[:, 1], zero_indices[:, 0]] = 0.0
        L_init = np.linalg.cholesky(P_init)
        u_init = np.asarray(cholesky_to_vector(L_init))
        assert not np.any(np.isnan(u_init)), "NaN in initial state"
        assert abs(constr(u_init, zero_indices).max()) < 1e-8, "Constraint violation"
        init_states.append(u_init)

    print(f"Created {len(init_states)} initial states")
    print("Running constrained HMC sampling...")
    sys.stdout.flush()

    # Time the sampling
    start_time = time.perf_counter()

    # Run sampling
    results = mici.sample_constrained_hmc_chains(
        n_warm_up_iter=n_warm_up_iter,
        n_main_iter=n_main_iter,
        init_states=init_states,
        neg_log_dens=lambda u: neg_log_dens(u, n_obs, S_chol),
        constr=lambda u: constr(u, zero_indices),
        backend="jax",
        seed=rng,
        monitor_stats=("accept_stat", "n_step", "step_size"),
        trace_funcs=[trace_func],
        n_worker=1,
        use_thread_pool=False,
    )

    end_time = time.perf_counter()
    sampling_duration = end_time - start_time

    print(f"\nSampling complete! Duration: {sampling_duration:.2f} seconds")
    sys.stdout.flush()


    ess = arviz.ess(results.traces, var_names=["u"])
    r_hat = arviz.rhat(results.traces, var_names=["u"])

    min_ess = float(ess.min().u.data)
    max_rhat = float(r_hat.max().u.data)

    # Extract posterior samples of precision matrix
    # Shape: (n_chain, n_main_iter, n_variable, n_variable)
    P_samples = np.array(results.traces["P"])
    P_mean = P_samples.mean(axis=(0, 1))

    print(f"Min ESS: {min_ess:.1f}, Max R-hat: {max_rhat:.3f}")
    print(f"P_samples shape: {P_samples.shape}")
    sys.stdout.flush()

    # Save full posterior samples as numpy array
    # Reshape to 2D for RcppCNPy compatibility (doesn't support 4D arrays)
    # Original shape: (n_chain, n_main_iter, n_variable, n_variable)
    # Saved shape: (n_chain * n_main_iter, n_variable * n_variable)
    P_samples_flat = P_samples.reshape(-1, n_variable * n_variable)
    np.save(samples_file, P_samples_flat)
    print(f"Posterior samples saved to: {samples_file} (flattened shape: {P_samples_flat.shape})")

    # Helper to convert NaN/Inf to None for JSON compatibility
    def sanitize_for_json(val):
        if isinstance(val, float) and (np.isnan(val) or np.isinf(val)):
            return None
        return val

    # Save summary results as JSON
    output_data = {
        "min_ess": sanitize_for_json(min_ess),
        "max_rhat": sanitize_for_json(max_rhat),
        "P_mean": P_mean.tolist(),
        "n_chain": n_chain,
        "n_iter": n_main_iter,
        "samples_file": samples_file,
        "sampling_duration_seconds": sampling_duration,
    }

    with open(output_file, "w") as f:
        json.dump(output_data, f)

    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
