# 20532 Macroeconometrics — Problem Sets

This repository contains MATLAB solutions, explanatory notes, and generated results for **20532 Macroeconometrics**. It currently includes Problem Sets 1 and 2; the PS2 solution follows the **Fall 2026** assignment, including Question 1(d).

## Repository structure

```text
20532-macroeconometrics-ps/
├── .gitattributes
├── .gitignore
├── README.md
├── ps1/
│   ├── explanatory_note.pdf
│   ├── code/
│   │   ├── PS1_complete_solution.m
│   │   └── arma_generator_ps1.m
│   └── output/                 # Included figures, CSV tables and MAT results
└── ps2/
    ├── explanatory_note.pdf
    ├── code/
    │   └── PS2_complete_solution.m
    └── output/                 # Input data, figures, CSV tables and MAT results
```

The tree summarizes the repository structure. See [Running the scripts](#running-the-scripts) for the locations used when regenerating results.

## Problem Set 1

Main script: [PS1_complete_solution.m](ps1/code/PS1_complete_solution.m). The [explanatory note](ps1/explanatory_note.pdf) discusses the methods and results.

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

The supporting function, [arma_generator_ps1.m](ps1/code/arma_generator_ps1.m), is stored alongside the main script. Included results are in [ps1/output](ps1/output/).

## Problem Set 2 — Fall 2026

Main script: [PS2_complete_solution.m](ps2/code/PS2_complete_solution.m). It contains four question sections and one local function for OLS estimation of a VAR. The [explanatory note](ps2/explanatory_note.pdf) discusses the economic interpretation.

### Question 1 — Unit-root testing

A driftless random walk with `T = 250` and innovation variance 1 is simulated 10,000 times. The exercise covers:

- **(a)** Empirical OLS slope distributions, with and without an intercept.
- **(b)** The Dickey–Fuller t-statistic and rejection frequencies using a one-sided standard Normal critical value.
- **(c)** Empirical t-statistic percentiles compared with the Dickey–Fuller values tabulated in Hamilton.
- **(d)** The joint null `alpha = 0` and `rho = 0`, with `rho = phi - 1`. The two-restriction Wald statistic `W = 2F` is compared with chi-square(2), conventional F, and the non-standard Dickey–Fuller `Phi_1` reference distribution.

Part (d) implements the joint test in the **Fall 2026** assignment. Its 5% `Phi_1` critical value is 4.63 for the F statistic. The saved [rejection table](ps2/output/question1_rejections.csv) reports rejection frequencies and Monte Carlo standard errors.

### Question 2 — Spurious regression

The four cases discussed in Enders are simulated with 5,000 replications per case at `T = 100, 500, 1000`:

1. Two independent stationary series.
2. An I(1) dependent variable and an independent I(0) regressor.
3. Two independent I(1) series.
4. Two cointegrated I(1) series sharing a common stochastic trend.

The script stores regression coefficients, conventional t-statistics and R-squared values. In the cointegrated case the true slope is **one**: rejection of a zero slope measures power, so the test around the true slope is reported separately. The [summary table](ps2/output/question2_spurious_summary.csv) compares all cases and sample sizes.

### Question 3 — Invertibility and impulse responses

The assignment's two-variable DGP is simulated with `beta = 0.6`, shock variances 1 and 0.8, and 500 observations per dataset. Each of 1,000 replications estimates a VAR(4) with an intercept and Cholesky ordering `[x, y]`.

Estimated mean IRFs and their 2.5th and 97.5th sampling percentiles are compared with analytical true IRFs over horizons 0–20. Both use one-standard-deviation shocks. The script also checks the equivalence of the recursive two-equation representation and reduced-form OLS plus Cholesky identification.

The determinant of the DGP's MA polynomial has roots 1 and -0.6. The interior root explains why the structural shocks are non-fundamental and why recursive VAR responses differ from the true responses. See the [IRF comparison](ps2/output/question3_irf_comparison.png) and [numerical summary](ps2/output/question3_irf_summary.csv).

### Question 4 — Granger causality

The included [Romer–Romer workbook](ps2/output/Romer_Romer.xlsx) contains US inflation, unemployment, the federal funds rate and Romer–Romer monetary policy shocks over **1969Q1–1996Q4**.

A four-variable VAR(4) with an intercept is estimated over 1970Q1–1996Q4, after accounting for the lags. Tests assess whether RR shocks predict each macro variable and whether each macro variable predicts RR shocks, conditional on the full VAR. A further joint test considers all three macro variables in the RR equation.

The [Granger-causality table](ps2/output/question4_granger_tests.csv) reports Wald and F statistics, p-values and rejection decisions. These tests measure predictive content; failure to reject does not establish shock exogeneity.

### Included PS2 outputs

The [ps2/output](ps2/output/) directory contains:

| Format | Contents |
| --- | --- |
| 3 PNG figures | Unit-root distributions, the four spurious-regression cases, and estimated versus true IRFs. |
| 7 CSV tables | Four Q1 summaries, one Q2 summary, one Q3 IRF summary and one Q4 test table. |
| 4 MAT files | Q1–Q3 Monte Carlo draws and summaries, plus Q4 data and VAR estimates. |
| 1 XLSX workbook | Input data for Q4, `Romer_Romer.xlsx`. |

## Requirements and reproducibility

- **MATLAB R2021a or newer** for PS2.
- **Statistics and Machine Learning Toolbox** for the distribution functions and percentiles used in PS2.
- PS2 estimates VARs with ordinary matrix operations; it does not require Econometrics Toolbox.

PS2 uses separate `twister` seeds for each Monte Carlo section: `2053201`, `2053202` and `2053203` for Q1, Q2 and Q3, respectively. Sample sizes and replication counts are defined in the corresponding sections. Scripts derive paths from their own location rather than a user-specific absolute path.

## Running the scripts

Clone or download the repository, open MATLAB, and set the Current Folder to the repository root.

### Run PS1

```matlab
run(fullfile('ps1', 'code', 'PS1_complete_solution.m'));
```

The current PS1 script creates `ps1/outputs/` for newly generated results. The results already included in this repository are stored in `ps1/output/` (singular).

### Run PS2

The current PS2 script reads the workbook from the directory containing the script and also saves regenerated outputs there. Copy the included input workbook alongside the script before running:

```matlab
copyfile(fullfile('ps2', 'output', 'Romer_Romer.xlsx'), ...
         fullfile('ps2', 'code', 'Romer_Romer.xlsx'));
run(fullfile('ps2', 'code', 'PS2_complete_solution.m'));
```

The run creates the figures, CSV tables, MAT files and `PS2_command_window_log.txt` in `ps2/code/`. The included result copies remain in `ps2/output/`. The adjacent workbook makes the run independent of the original course folder; the script also supports a fallback search for that folder when available.

Rerunning either script replaces its generated files with the same names. Local command-window logs are excluded by `.gitignore`.

## Collaboration workflow

The `main` branch contains the accepted version. Changes should be proposed through a fork and Pull Request, then reviewed before they are merged into `main`.

When preparing an update for GitHub, include the relevant MATLAB source, explanatory note and intended result files. Keep relative links and directory capitalization consistent with the committed paths, since GitHub paths are case-sensitive.
