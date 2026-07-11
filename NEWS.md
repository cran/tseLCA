# tseLCA 1.0.0

-   Initial submission to CRAN.

## Core Estimation Framework

-   Implemented BCH and ML bias-adjusted three-step estimators for latent class analysis (LCA).
-   Added support for structural models containing covariates ($Z_p$), distal outcomes ($Z_o$), and combined models (estimating the relationship between $Z_p$ and the latent class first, followed by the distal outcome adjusting for covariate-adjusted posteriors).
-   Implemented analytic sandwich variance estimation to correctly propagate measurement uncertainty from the first-step LCA through classification-error correction in the final step.
-   Added a robust standard error option (`use.simple.cov = TRUE`) that bypasses the measurement-uncertainty correction for faster computation in large, well-separated samples.

## Measurement Model (Step 1) Integration

-   Integrated with the 'multilevLCA' package for efficient Step-1 measurement model estimation.
-   Added support for polytomous indicator items (0-based integer coding).
-   Implemented Full Information Maximum Likelihood (FIML) to handle missing data in the measurement model via the `incomplete = TRUE` argument (using a two-pass row-filtering strategy).
-   Added the ability to pass a pre-fitted measurement model (via the `step1` argument) to reuse across multiple structural models or apply to different sample subsets.
-   Implemented automated random restarts for the measurement model triggered when entropy $R^2$ falls below a user-specified threshold.

## Algorithmic Flexibility & Structural Models

-   Added support for both modal and proportional (soft) posterior class assignment (`use.modal.assignment`).
-   Integrated Gaussian, Poisson, and binomial families for distal outcome estimation.
-   Added the `rebase` argument to allow users to easily change the reference latent class for the multinomial logit parameterization while maintaining invariant log-likelihoods.
-   Implemented two-step EM estimation (`fitZ_from_fit0()`) to generate stable starting values for the three-step structural model.

## Utilities and Methods

-   Included standard S3 methods for `tseLCA` objects: `summary()`, `coef()`, `vcov()`, and `plot()` (which delegates to 'multilevLCA' for item-profile visualization).
-   Built a data-generating process (`generate_data()`) that replicates the Bakk & Kuha (2018) simulation study design for both covariates and distal outcomes under varying separation conditions.

# tseLCA 1.0.1

## CRAN resubmission
- Removed single quotes around acronyms in DESCRIPTION; added explanations
  of BCH, ML, and LCA.
- Replaced `T`/`F` with `TRUE`/`FALSE` throughout internal codebase
- Added `\value` tags to all exported functions missing them, including `bk2018_params`.
- `inst/examples`: examples now write to `tempdir()` instead of the home
  filespace.
- `inst/examples`: commented out `rm(list = ls())` calls.
- `inst/examples`: commented out `install.packages()` calls.

# tseLCA 1.0.2

## CRAN resubmission
- Title: corrected hyphenated compound capitalization to "Three-Step" per
  reviewer request.