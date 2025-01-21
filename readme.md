#### Readme

This file describes all the files needed for the replication of the Monte Carlo simulations and Boston Housing dataset application presented in  manuscript-Tables 1 and 3, and supplement-Tables 2-15 in the paper by Barde, Cherodian and Tchuente "Moran's I Lasso for models with spatially correlated data".

The replication package contains 1 Matlab and 3 R code files:
- 'MiLasso_bounds_sens.m' replicates the sensitivity analysis table (supplement-Table 1). Expected computational time is $<30$ seconds. Tables printed as LaTex code in the command window.
- 'MiLasso_selcomp.R' replicates the main Monte Carlo simulation results which compare Mi-Lasso, CV-Lasso and FstepZ for both set-up A and B (manuscript-Table 1 and supplement-Table 2-13). Expected computational time is 32 hours. Tables printed as LaTex code in console.
- 'MiLasso_comptime.R' replicates the computational time table (supplement-Table 14). Expected computational time is 18 hours. Table printed as in console.
- 'MiLasso_bhousing_app.R' replicates the regression results of the Boston housing application (manuscript-Table 3 and supplement-Table 15). Please note the data is from the package 'spdep'. Expected computational time is $<30$ seconds. Regression output printed in console.

The replication package contains a 'renv} lockfile enabling replicators to install the exact packages (glmnet, penalized, xtable, spdep, lmtest, sandwich) and related dependencies used in the generation of the main results, with details provided below. The exact package versions used in the 'renv} lockfile are listed in \url{/logs/sessionInfo_win.txt}.

Finally, the replication package contains a subfolder \url{/logs} which contains session information and replication logs for 'MiLasso_selcomp.R} and 'MiLasso_bhousing_app.R} on 3 different platforms:
- Windows 10 build 19045
- macOS Ventura 13.7.2
- Linux elementary OS 6.1 (Jolnir), which is based on Ubuntu 20.04

**Note:** This multiplatform replication exercise reveals small and occasional variations in the Monte Carlo tables produced by the R scripts even when the exact R environment is locked using 'renv'. Examination of the three 'sessionInfo' log files reveals that this is due to differing versions of the BLAS and LAPACK utilities provided by the different platforms. The results displayed in the main paper and online supplement are the ones obtained using the Windows 10 platform.

**Data Availability Statement**: No data is provided with this replication package, the data used in  the Boston Housing application (results in manuscript-Table 3 and supplement-Table 15) is from the R package 'spdep' \citep{bw18}. The data is loaded by the command 'data(boston)'}' and is from \cite{hr78_bh} corrected for a few minor errors and augmented with the latitude and longitude by \cite{gp96_sbh}.

### Instructions for Replicators

Using the 'renv' lockfile requires using R version 4.4.1 or 4.4.2 and installing the 'renv' package using 'install.packages("renv"))', then running the command 'renv::restore()'' to install the correct dependencies. Once the environment is setup, the 3 R files can be run. The structure of the four scripts is as follows:

- 'MiLasso_bounds_sens.m' replicates the sensitivity analysis table (supplement-Table 1). The file structure is as follows:
  1. Set parameter values for sensitivity analysis.
  2. Preallocate output matrix and iterate over sensitivity dimensions
  3. Latex code for supplement Tables 1 is then presented in console.
  
- 'MiLasso_selcomp.R' reproduces the main simulation results presented in the appendix for both set-ups A and B (supplement Tables 2-13). Table 1 in the main paper appends the results presented in Supplement Tables 3, 6 and 9, with the results for `Naïve Mi-Lasso' and `Naïve CV-Lasso' omitted. The file structure is as follows:
  1. Results for set-up A are estimated and saved.
  2. Latex code for supplement Tables 2-10 is then presented in console.
  3. Latex code for main paper Table 1 is then printed in console.
  4. Results for set-up B are estimated and saved.
  5. Latex code for supplement Tables 11-13 is then printed in console.

- 'MiLasso_comptime.R' replicates the computational time table (supplement-Table 14). The file structure is as follows:
  1. Functions and simulation set-up is specified.
  2. Computational times are saved in a matrix
  3. Computational times printed in console.
  4. Computational times relative to Mi-Lasso printed in console.

- 'MiLasso_bhousing_app.R' replicates the regression results of the Boston housing application (manuscript-Table 3 and supplement-Table 15). Manuscript-Table 3 is the same as supplement-Table 15 but with rows corresponding to parameters not discussed in the main paper omitted. The file structure is as follows:
1. Data is loaded and cleaned.
2. Results for each of the estimator considered estimated.
3. Regression results for each of the estimators printed in console.
