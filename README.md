# tseLCA

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

## Overview

**tseLCA** (*Three-Step Estimation for Latent Class Analysis*) introduces bias-adjusted three-step estimation for structural latent class models in R. The package provides a comprehensive framework for estimating latent class models with covariates and distal outcomes while preserving the measurement structure established during class formation.

Building upon the efficient measurement-model estimation procedures implemented in **multilevLCA**, **tseLCA** extends existing functionality through modern three-step estimators, classification-error corrections, and variance estimation procedures that appropriately account for uncertainty from the latent class measurement stage.

## Key Features

### Bias-Adjusted Three-step Estimation

**tseLCA** is the first R package to provide a unified implementation of modern bias-adjusted three-step estimators for latent class analysis. In contrast to traditional one-step approaches (implemented by the popular package, **poLCA**), where the inclusion of covariates may alter the underlying latent class definitions, three-step methods preserve the measurement model estimated in the first stage and subsequently adjust for classification error when estimating structural relationships.

The package implements both BCH- and ML-based three-step estimators with sandwich variance estimators that propagate uncertainty from the measurement model through the classification-error correction process.

### Flexible Measurement and Structural Samples

Unlike conventional latent class software that uses a one-step estimation approach, **tseLCA** allows measurement and structural models to be estimated using different datasets. This flexibility enables researchers to calibrate a measurement model on a primary or reference sample and subsequently apply the resulting class definitions to an external dataset.

### Support for Multiple Distal Outcome Types

**tseLCA** provides native support for a broad range of distal outcome distributions, including:

-   Continuous outcomes (Gaussian)
-   Count outcomes (Poisson)
-   Binary outcomes (Bernoulli).

### Automated Model Optimization

Latent class estimation is often susceptible to local maxima and convergence to suboptimal solutions. To improve estimation reliability, **tseLCA** incorporates automated diagnostic procedures that monitor model quality during measurement-model estimation.

### Missing Data Handling

Following a similar approach as **multilevLCA**, **tseLCA** employs Full-Information Maximum Likelihood (FIML) estimation to accommodate partially observed response patterns without discarding incomplete observations.

## Installation

You can install the development version of tseLCA from GitHub like so:

``` r
# Install developmental tseLCA from the GitHub repository
if (!require("pak")) {
  install.packages("pak")
}

pak::pak("SamLeeBYU/tseLCA")
```

Once tseLCA is on CRAN, then you can simply install it from a CRAN server.

``` r
# Install tseLCA from CRAN
install.packages("tseLCA")
```

Then read the introductory vignette on this package's webpage here: <https://SamLeeBYU.github.io/tseLCA/articles/tseLCA-workflow.html>

## Example

This is a basic example which shows you how to simulate data and run a three-step LCA with a covariate in a single function call:

``` r
library(tseLCA)

# 1. Generate synthetic data 
# (3 classes, 6 dichotomous items, and a multinomial logit covariate 'Zp')
d <- generate_data(
  n = 500, 
  separation = "high", 
  scenario = "covariate", 
  seed = 1
)

# 2. Estimate the three-step model
# This automatically fits the measurement model and estimates covariate effects
fit <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  # Proportional assignment is recommended for better uncertainty propagation
  use.modal.assignment = FALSE
)

# 3. View the measurement and structural model estimates
summary(fit)
```