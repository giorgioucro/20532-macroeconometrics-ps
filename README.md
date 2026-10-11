\# 20532 Macroeconometrics — Problem Sets



This repository contains the MATLAB code and generated outputs for the Problem Sets of \*\*20532 Macroeconometrics\*\*.



\## Repository structure



```text

20532-macroeconometrics-ps/

├── PS1/

│   ├── code/

│   │   ├── PS1\_complete\_solution.m

│   │   └── arma\_generator\_ps1.m

│   └── outputs/

├── PS2/

├── PS3/

├── README.md

└── .gitignore

```



The repository will be expanded as the subsequent Problem Sets are completed.



\## Problem Set 1



The main script is:



`PS1/code/PS1\_complete\_solution.m`



It contains the complete solution to Questions 1–6.



\### Question 1 — Simulation of an AR(1)



Simulates an AR(1) process with two alternative methods:



\- an explicit recursive `for` loop;

\- MATLAB's `filter` function.



Using the same innovation sequence and the same initial condition, the two approaches generate the same realization up to numerical precision.



\### Question 2 — Initial conditions and burn-in



Studies the effect of initializing a stationary AR(1) process far from its unconditional mean.



Since the effect of the initial condition disappears over time when `|phi| < 1`, additional observations are simulated and discarded as a burn-in period. The retained sample is therefore effectively independent of the arbitrary starting value.



\### Question 3 — Monte Carlo distribution of the OLS estimator



Studies the empirical distribution of the OLS estimator for an AR(1) process with `phi = 0.4` and `T = 250`.



The experiment uses 5,000 Monte Carlo replications. In each replication, an intercept and the AR coefficient are estimated by OLS. A two-sided test of



`H0: phi = 0`



is also performed and its empirical rejection frequency is computed.



\### Question 4 — Effect of sample size



Studies the finite-sample distribution of the OLS estimator for a persistent AR(1) process with `phi = 0.9`.



The experiment uses 1,000 simulations for each sample size:



`T = 50, 100, 200, 1000`.



The resulting distributions illustrate how the bias and dispersion of the estimator change as the sample size increases.



\### Question 5 — Fitting an AR(1) to MA(1) data



Generates data from an MA(1) process with `theta = 0.6` and estimates the regression



`x\_t = a x\_(t-1) + v\_t`



without an intercept.



The Monte Carlo experiment studies the distribution of the OLS estimator for different sample sizes and shows that its mean converges to the population linear-projection coefficient



`theta / (1 + theta^2)`.



\### Question 6 — General ARMA(p,q) generator



Implements a general loop-based function for simulating ARMA(p,q) processes.



The function accepts:



\- the number of observations;

\- the variance of the white-noise innovations;

\- AR and MA polynomials supplied either as coefficients or as roots.



It returns the simulated ARMA realization and the corresponding innovations. The implementation is tested on an ARMA(1,2) example and validated both across the coefficient/root parameterizations and against MATLAB's `filter` function.



\## Supporting function



`PS1/code/arma\_generator\_ps1.m`



This file contains the general ARMA(p,q) simulation function used in Question 6.



\## Outputs



Running the main script automatically generates the relevant tables, figures, CSV files, and MATLAB result files in:



`PS1/outputs/`



The code uses paths relative to the location of the main script, so the repository can be cloned and run without requiring user-specific absolute paths.



\## How to run PS1



1\. Clone or download the repository.

2\. Open MATLAB.

3\. Open `PS1/code/PS1\_complete\_solution.m`.

4\. Run the script.

5\. The generated outputs will be saved automatically in `PS1/outputs/`.



\## Collaboration workflow



The `main` branch contains the current accepted version of the project.



Contributions should be proposed through Pull Requests:



1\. Fork the repository.

2\. Make the desired changes in the fork.

3\. Commit and push the changes.

4\. Open a Pull Request into the original repository.

5\. Review the proposed changes before merging them into `main`.



This workflow keeps the accepted version of the project under review while allowing collaborators to propose modifications.



\## Software



\- MATLAB

\- Git

\- GitHub

