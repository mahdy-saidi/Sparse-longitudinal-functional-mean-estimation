# Estimation of the mean function in sparse longitudinal functional data

This repository contains R code for Monte Carlo experiments related to the estimation of the mean function in sparse longitudinal functional data.

The project studies observations of the form

$$
Y_{ij}(u) = \mu(u, T_{ij}) + X_i(u, T_{ij}) + \varepsilon_{ij}(u),
$$

where $Y_{ij}(u)$ is a functional observation for subject $i$ at visit time $T_{ij}$, $\mu(u,t)$ is the mean function of interest, $X_i(u,t)$ is a subject-specific stochastic process, and $\varepsilon_{ij}(u)$ is an additive noise process.

The main goal is to compare several estimators of the mean function $\mu(u,t)$ in a sparse longitudinal setting, using simulation experiments.


## Repository structure

```text
.
├── README.md
├── script.R
└── outputs/
```

The main R script is `script.R`.

## Requirements

The code is written in R.

The package `tikzDevice` is optional. If it is installed and a working LaTeX distribution with the `pgf` package is available, figures are exported as TikZ files. Otherwise, the code falls back to PDF output.

To install the optional package:

```r
install.packages("tikzDevice")
```

For TikZ output, a LaTeX installation is also required. With TinyTeX, the required LaTeX packages can be installed using:

```r
tinytex::tlmgr_install(c("pgf", "preview"))
```

## How to run the code

### Quick run

For a faster test run, use:

```bash
SIM_SETTINGS=QUICK SIM_OUTDIR=./outputs Rscript script.R
```

This mode uses smaller Monte Carlo sample sizes and is recommended for checking that the code runs correctly.

### Full run

For the full simulation study, use:

```bash
SIM_SETTINGS=FULL SIM_OUTDIR=./outputs Rscript script.R
```

This mode is computationally heavier and may take significantly longer.

## Output files

The output directory is controlled by the environment variable `SIM_OUTDIR`.

For example:

```bash
SIM_OUTDIR=./outputs
```

The script can generate:

* `.tex` files containing TikZ figures or LaTeX tables;
* `.pdf` figures when TikZ is unavailable;
* `.csv` files with raw numerical results if CSV export is enabled.

By default, the script uses:

```r
WRITE_PDF <- FALSE
WRITE_CSV <- FALSE
```

To save additional PDF figures or CSV files, change these options directly in the R script:

```r
WRITE_PDF <- TRUE
WRITE_CSV <- TRUE
```

## Simulation settings

The simulation size is controlled by:

```r
SETTINGS <- Sys.getenv("SIM_SETTINGS", unset = "FULL")
```

Available modes are:

* `QUICK`: small simulation, useful for testing;
* `FULL`: full Monte Carlo simulation, used for final results.

The master seed is set by:

```r
BASESEED <- 2026L
```

This ensures reproducibility of the simulation experiments.

## Statistical methods

The code compares several estimators of the mean function:

1. uniform weights;
2. subject-balanced weights;
3. Monte Carlo control-neighbors weights;
4. spacing-based estimator;
5. plug-in estimator using an estimated design density;
6. bivariate local-linear benchmark.

The estimation strategy is based on a series expansion in the time direction:

$$
\mu(u,t) \approx \sum_{k=1}^K \beta_k(u)\phi_k(t).
$$

The simulation study investigates the effect of the truncation level $K$, the sparse longitudinal design, the visit-time distribution, and the weighting scheme on the integrated squared error.

## Notes

The full simulation may require substantial computation time. It is recommended to start with `SIM_SETTINGS=QUICK` before running the full version.

The generated results are intended to support a methodological study on mean function estimation for sparse longitudinal functional data.
