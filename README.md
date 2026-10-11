# 20532 Macroeconometrics — Problem Sets

This repository contains MATLAB code and generated outputs for the problem sets of 20532 Macroeconometrics.

## Repository structure

```text
20532-macroeconometrics-ps/
├── .gitattributes
├── .gitignore
├── README.md
└── ps1/
    ├── code/
    │   ├── PS1_complete_solution.m
    │   └── arma_generator_ps1.m
    └── outputs/
        ├── PS1_all_results.mat
        ├── question1_loop_vs_filter.png
        ├── question1_simulated_series.csv
        ├── question1_summary.csv
        ├── question2_initial_condition_and_burnin.png
        ├── question2_stationary_path.csv
        ├── question2_summary.csv
        ├── question3_monte_carlo_draws.csv
        ├── question3_ols_distribution.png
        ├── question3_summary.csv
        ├── question4_distributions_by_sample_size.png
        ├── question4_monte_carlo_draws.csv
        ├── question4_summary.csv
        ├── question5_distributions_by_sample_size.png
        ├── question5_mean_convergence.png
        ├── question5_monte_carlo_draws.csv
        ├── question5_summary.csv
        ├── question6_arma_example.csv
        ├── question6_arma_example.png
        └── question6_summary.csv
```

## Problem Set 1

The main script is `ps1/code/PS1_complete_solution.m`.

### Question 1 — Simulation of an AR(1)

An AR(1) process is simulated using both a recursive `for` loop and MATLAB's `filter`. The two methods are compared using the same innovations and initial condition.

### Question 2 — Initial conditions and burn-in

The exercise studies an AR(1) process starting far from its unconditional mean and uses a burn-in period to obtain a proper stationary realization.

### Question 3 — Monte Carlo distribution of the OLS estimator

Monte Carlo simulations are run for an AR(1) process with `phi = 0.4` and `T = 250`. Each simulation estimates an intercept and the AR coefficient by OLS and tests `H0: phi = 0`.

### Question 4 — Effect of sample size

The OLS distribution is studied for `phi = 0.9` using `T = 50, 100, 200, 1000` and 1,000 simulations for each sample size.

### Question 5 — Fitting an AR(1) to MA(1) data

MA(1) data with `theta = 0.6` are generated, and the no-intercept regression `x_t = a x_(t-1) + v_t` is estimated. The Monte Carlo mean approaches `theta / (1 + theta^2)`.

### Question 6 — General ARMA(p,q) generator

`arma_generator_ps1.m` implements a general loop-based ARMA(p,q) simulator. It accepts polynomial coefficients or roots, returns the simulated series and innovations, and is validated against MATLAB's `filter`.

## Supporting function

The supporting function is `ps1/code/arma_generator_ps1.m`.

## Outputs

Generated figures, CSV files, and MATLAB result files are saved in `ps1/outputs/`. All paths are relative, so the project does not depend on user-specific absolute paths.

## How to run PS1

1. Clone or download the repository.
2. Open MATLAB.
3. Open `ps1/code/PS1_complete_solution.m`.
4. Run the script.
5. Outputs are saved automatically in `ps1/outputs/`.

## Collaboration workflow

The `main` branch contains the accepted version. Changes should be proposed through a fork and Pull Request, then reviewed before they are merged into `main`.

## Software

- MATLAB
- Git
- GitHub
