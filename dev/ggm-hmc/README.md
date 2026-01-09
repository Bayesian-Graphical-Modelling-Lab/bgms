# GGM Constrained HMC Sampling

This folder contains scripts for sampling precision matrices for Gaussian Graphical Models (GGMs) using constrained Hamiltonian Monte Carlo (HMC). The constrained HMC approach enforces exact zeros in the precision matrix (corresponding to missing edges in the graph) as hard constraints.

## Overview

- `run_constrained_hmc.py` - Python script that performs the constrained HMC sampling using JAX and [mici](https://github.com/matt-graham/mici)
- `run_constrained_hmc_subprocess.R` - R script that generates data, calls the Python sampler via subprocess, and processes results

## Prerequisites

### 1. Clone the mici example repository

The Python environment is defined in a separate repository. Clone it as a sibling folder:

```bash
cd /path/to/bgms/dev
git clone https://github.com/matt-graham/ggm-precision-constrained-hmc
```

Your folder structure should look like:
```
bgms/dev/
├── ggm-hmc/                           # This folder
│   ├── README.md
│   ├── run_constrained_hmc.py
│   └── run_constrained_hmc_subprocess.R
└── ggm-precision-constrained-hmc/     # Cloned repo with uv environment
    ├── .venv/
    ├── pyproject.toml
    └── ...
```

**Alternative:** If you prefer a different location, adjust the `python_dir` path in `run_constrained_hmc_subprocess.R`.

### 2. Set up the Python environment with uv

Install [uv](https://docs.astral.sh/uv/) if you don't have it:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then create and sync the Python environment:

```bash
cd /path/to/bgms/dev/ggm-precision-constrained-hmc
uv sync
```

This will create a `.venv` folder with Python 3.14 (free-threaded) and all required dependencies.

## R Package Dependencies

The R script requires:
- `bgms` (this package)
- `mvtnorm`
- `jsonlite`
- `RcppCNPy` (for reading numpy arrays)
- `BDgraph` (for generating G-Wishart samples)

```r
install.packages(c("mvtnorm", "jsonlite", "RcppCNPy", "BDgraph"))
```

## Usage

```r
# Set working directory to bgms root
setwd("/path/to/bgms")

# Run the script
source("dev/ggm-hmc/run_constrained_hmc_subprocess.R")
```

The script will:
1. Generate simulated data from a sparse GGM
2. Run the bgms MH sampler for comparison
3. Call the Python constrained HMC sampler
4. Load and summarize the results

