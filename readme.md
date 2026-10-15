# vicsek-cluster-interactions

[![Julia](https://img.shields.io/badge/Julia-1.12-9558B2?style=flat&logo=julia)](https://julialang.org)
[![Paper](https://img.shields.io/badge/Paper-PRR_Submitted-blue)](https://github.com/dynamicalsystemslaboratory/vicsek-cluster-interactions)

This repository contains the numerical implementation, data generation scripts, and analysis pipeline for the paper:

> **Cluster behavior determines global order in the Vicsek model**\
> Michael Napoli, Ofek Lauber Bonomo, Ashiqur Rahman, Nicole Abaid, and Maurizio Porfiri\
> *New York University & Virginia Tech*

---

## Overview

The Vicsek model describes self-propelled particles updating their heading based on local neighbor interactions and noise. This repository provides tools to analyze Vicsek dynamics through a Smoluchowski-like coagulation-fragmentation kinetic framework. By investigating cluster-cluster collision, internal splitting (fragmentation), and geometric scaling (convex hull perimeters and areas), this codebase bridges microscopic interaction rules and macroscopic phase transitions.

---

## Repository Structure

```text
vicsek-cluster-interactions/
├── args.jl                         # Command-line argument parsing for simulations
├── prob.jl                         # Core problem definition for Vicsek dynamics & ODEs
├── 0-params/                       # Parameter sweep configuration lists
│   ├── nos-list.jl                 # Noise amplitude (η) sweeps
│   ├── rho-list.jl                 # Density (ρ) sweeps
│   └── spd-list.jl                 # Particle speed (s) sweeps
├── 1-order/                        # Global order dynamics & phase transition analysis
│   ├── order-0-gen.ipynb           # Initial Vicsek state generation
│   ├── order-1-gen-*/              # Order parameter data generation scripts
│   ├── order-2-crit-*.ipynb        # Critical noise analysis (η_c)
│   └── order-3-cmd-*.ipynb         # Command-line batch execution wrappers
├── 2-cluster/                      # Coagulation-fragmentation rate & kernel extraction
│   ├── cluster-1-gen*/             # Cluster trajectory & interaction data generation
│   ├── cluster-2-calib*.ipynb      # Parameter triplet (κ, λ_f, λ_c) calibration to Ω
│   └── cluster-3-freq*.ipynb       # Empirical transition matrix Ω plotting and raw measurements
├── 3-geom/                         # Cluster geometry & scale-invariance computations
│   ├── geom-1-gen/                 # Cluster perimeter & area data generation at critical noise
│   ├── geom-2-calib*.ipynb         # Poly-logarithmic function calibration to geometry trends
│   └── geom-3-freq*.ipynb          # Geometry calibration parameter result plotting
├── 4-master/                       # Coagulation-fragmentation equation testing on Vicsek data
├── 5-size/                         # Master equation calibration, size scaling, and extrapolation
│   ├── size-0-params.jl            # Speed and density subset for extrapolation analysis
│   ├── size-1-gen-order/           # Order parameter data generation across system sizes (N)
│   ├── size-2-crit-*.ipynb         # Critical noise evaluation across system sizes
│   ├── size-3-gen-cluster/         # Cluster interaction data generation across system sizes
│   ├── size-4-calib-cluster.ipynb  # Transition probability map (Ω) calibration across system sizes
│   ├── size-5-gen-geom/            # Cluster geometry data generation across system sizes
│   ├── size-6-calib-geom.ipynb     # Geometry scaling parameter calibration across system sizes
│   ├── size-7-mean-cluster.ipynb   # Extraction of size-change probability maps
│   ├── size-8-calib-kernel.ipynb   # Coagulation-fragmentation kernel calibration
│   ├── size-9-kernel.ipynb         # System-size scaling and decay strength analysis of kernel
│   ├── size-10-master.ipynb        # Kinetic master equation simulation & cluster distribution extraction
│   ├── size-11-extrap*/            # Extrapolation simulations for larger system size (N=1024)
│   └── size-12-compile.ipynb       # Extrapolation results compilation and kernel scaling plots
└── README.md                       # Repository documentation
```

---

## Computational Framework

* Vicsek Particle Simulation: Simulates $N$ particles moving at constant speed $s$ with periodic boundary conditions across density $\rho$ and noise amplitude $\eta$.
* Kinetic ODE Integration: Integrates the closed coagulation-fragmentation rate equations: $$\dot{n}_i(t) = \sum_{j=1}^{i-1} \frac{j}{i} C_{j,i-j} n_j(t) n_{i-j}(t) + \sum_{j=i+1}^N F_{j,i} n_j(t) - \left( \sum_{j=1}^{N-i} C_{i,j} n_j(t) + \frac{1}{2}\sum_{j=1}^{i-1} F_{i,j} \right) n_i(t).$$
Solved numerically to steady state using the 9th-order adaptive Runge-Kutta algorithm (Vern9) from Julia's DifferentialEquations.jl.

---

## Prerequisites and Installation

* Requirements: Julia `v1.12` or higher.

---

## Reproducing Paper Results

* All folders and their contents are labeled with indeces denoting the order of execution when reproducing results reported in the paper. More detailed discussion are below.

1.  Order-disorder transition (`1-order/`)

    1. `order-1-gen*/`: Simulate system dynamics to evaluate the mean order parameter $\langle \varphi \rangle$ across noise $\eta$, density $\rho$, and speed $s$. Variations of this sub-process also exist for sweeps in the density and speed in isolation.
    2. `order-2-crit-*.ipynb`: Using the order parameter data simulated in step $1$, compute the critical noise (estimated via the point where the trend cross the $\langle\varphi\rangle=0.5$ threshold).
    3. `order-3-cmd-*.ipynb`: Plot the order parameter results and compare for changes in the cluster size distribution.
        * Creates Figure 2 in the main text.

2. Cluster interaction analysis (`2-cluster/`)

    1. `cluster-1-gen*/`: Generate cluster interaction simulations necessary to compute the probability $\Omega$ from the data.
    2. `cluster-2-calib*.ipynb`: Calibrate the parameter triplet $(\kappa,\lambda_f,\lambda_c)$ to the probability mass functions describing $\Omega$ for each transition probability element $\omega_{i,j}$.
    3. `cluster-3-freq-rho.ipynb`: Plot the transition probability results, including its raw measurement and the calibrated fits.
        * Creates Figure 3 in the main text.

3. Cluster geometry analysis (`3-geom/`)

    1. `geom-1-gen/`: Generate cluster geometry data at critical noise across many pairs of the density and speed, $(\rho,s)$.
    2. `geom-2-calib*.ipynb`: Calibrate the poly-logarithmic function to the cluster perimeter and area trends.
    3. `geom-3-freq*.ipynb`: Plot the calibration parameter results for variations in both speed and density, $(\rho,s)$.
        * Creates Figure 4 in the main text.

4. Testing the coagulation-fragmentation equations (`4-master/`)
    * This folder contains files and scripts necessary for testing the coagulation-fragmentation equations on data extracted from the Vicsek model. It does not contribute any analysis directly necessary to compiling the results of the paper.

5. Calibration of master equations and projection in size (`5-size/`)

    0. `size-0-params.jl`: Speed and density subset to consider to extrapolation analysis.
    1. `size-1-gen-order/`: Generate data across the parameter set for variations in smaller system sizes (denoted by the variable Nlist).
    2. `size-2-crit-*.ipynb`: Compute the critical noise across the simulated changes in the system size, $N$.
    3. `size-3-gen-cluster/`: Generate cluster interaction data for smaller instances of the system size.
    4. `size-4-calib-cluster.ipynb`: Calibrate the probability map, $\Omega$, parameters across range of smaller system sizes.
    5. `size-5-gen-geom/`: Generate cluster perimeter and area data for smaller instances of the system size.
    6. `size-6-calib-geom.ipynb`: Calibrate the probability map, $\Omega$, parameters across range of smaller system sizes.
    7. `size-7-mean-cluster.ipynb`: Extract the approximation of the change-in-size probability maps across the cluster range in systems of varying sizes.
    8. `size-8-calib-kernel.ipynb`: Calibrate the coagulation-fragmentation kernel for a system of a given size. This script is repeated for the range of small sizes considered in the paper.
    9. `size-9-kernel.ipynb`: Compare the kernel parameter across the computed system sizes, and extract the decay strength of the parameter.
    10. `size-10-master.ipynb`: Run simulations of the master equations for a focal system size and save the extract cluster size distribution to the data folder. The `runsim` parameter at the top of the file allows one to choose whether to re-run the master equation simulation, or to use existing data for plotting.
    11. `size-11-extrap*`: Run simulations of the master equations for the extrapolated system size ($N=1024$). Data is then saved aas an approximate cluster size distribution to the data folder of the larger system size.
    12. `size-12-compile.ipynb`: Plot the results of the extrpolation, as well as the kernel parameter scaling.
        * Creates Figure 5 in the main text.
