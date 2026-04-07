
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Systematic measurement changes: A simulation study

<!-- badges: start -->

<!-- badges: end -->

This repository contains the code to reproduce the simulation study in
the paper  
*“Assessing systematic measurement changes in single-wave study data: A
simulation study”*.

The study compares seven statistical methods for quantifying systematic
changes in a sequence of measurements across different settings:  
- Autoregressive integrated moving average (ARIMA)  
- Fused lasso signal approximator (FLSA)  
- Generalized additive model (GAM)  
- Locally weighted scatterplot smoothing (LOWESS)  
- Moving average (MA)  
- Pruned exact linear time (PELT)  
- Piecewise regression (PR)

The simulation results can be interactively explored via an accompanying
[Shiny application](https://sys-measurement-changes.bips.eu/). The
source code for the app is available in a separate [GitHub
repository](https://github.com/bips-hb/shiny-sys-measurement-changes).

## 📁 Repository structure

- `functions/`: Contains function definitions for
  [01_simulation_study.R](01_simulation_study.R):
  - [algorithms.R](functions/algorithms.R): Statistical methods
  - [problems.R](functions/problems.R): Simulation settings  
- `renv/`: Contains project-specific environment configuration files.
- `results/`: Contains all outputs of the simulation study, including
  figures, tables, and additional files. Raw results are stored in
  [res.rds](results/res.rds).  
- `01_simulation_study.R`: Runs the simulation study using the
  `batchtools` R package. Results are saved in
  [res.rds](results/res.rds)
- `02_visualizations.R`: Generates [figures](results/figures) and
  [tables](results/tables) for the main manuscript.
- `03_additional_file_1.R`: Creates the additional Excel file
  [additional_file_1.xlsx](results/additional_file_1.xlsx) containing
  bias and mean squared error for different settings.
- `04_additional_file_2.Rmd`: R Markdown to create
  [additional_file_2.pdf](results/addtional_file_2.pdf) with the output
  from the [Shiny
  application](https://sys-measurement-changes.bips.eu/).
- `renv.lock`: Records package versions.

## 🚀 Reproducing the results

To reproduce the simulation study, run
[01_simulation_study.R](01_simulation_study.R). This executes the
simulations and saves the results to `results/res.rds`.

To recreate the figures and tables for the main manuscript, run
[02_visualizations.R](02_visualizations.R). This saves the figures under
`results/figures` and the tables under `results/tables`.

To recreate the additional files
[additional_file_1.xlsx](results/additional_file_1.xlsx) and
[additional_file_2.pdf](results/additional_file_2.pdf), run
[03_additional_file_1.R](03_additional_file_1.R) and
[04_additional_file_2.Rmd](04_additional_file_2.Rmd), respectively.

> Note: It is **not necessary** to run the simulation study to create
> the figures, tables, or additional files, as the results are already
> stored in `results/res.rds`.

## 📚 R environment

This project uses the `renv` R package to ensure reproducible R
environments. The exact package versions used for the simulation study
are recorded in [renv.lock](renv.lock). To restore the environment and
install all required packages:

### 1. Clone the repository

Clone the repository (or download the ZIP) to save
[renv.lock](renv.lock) locally.

### 2. Open the project

Open the project in R (e.g. in RStudio). Make sure
[renv.lock](renv.lock) is in the project root directory.

### 3. Install renv

Install `renv` if not already installed.

``` r
install.packages("renv")
```

### 4. Restore project library

Restore the project library using [renv.lock](renv.lock)

``` r
renv::restore()
```
