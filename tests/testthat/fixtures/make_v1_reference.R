# Regenerate tests/testthat/fixtures/v1_reference.rds.
#
# Run from the package root, against the version whose numbers should be
# pinned:
#   Rscript --vanilla tests/testthat/fixtures/make_v1_reference.R
#
# Originally generated from tseLCA 1.1.1 (commit b65d12d). Regenerate only for
# deliberate numerical changes, and document each one in NEWS.md.
#
# Regenerations:
#   - 2.0 Step 1b bug fixes. AIC/BIC of covariate and combined models change
#     (corrected parameter counts). Measurement-only posteriors are now
#     pinned (previously excluded: not in data-row order). New configurations
#     cov_poly_corrected and cov_sparse pin polytomous corrected SEs (NA in
#     1.1.1). All estimates and variances otherwise unchanged vs 1.1.1.
#   - 2.0 Step 1c: gaussian distal outcomes estimate the within-class
#     variance sigma2 (fixed at 1 in 1.1.1). Gaussian ML distal estimates and
#     variances, and gaussian log-likelihoods/AIC/BIC (ML and BCH), change;
#     sigma2 is now pinned. Other families unchanged.
#   - 2.0 Step 1d: ML distal likelihood over the expanded data (proportional
#     assignment weights outside the log, Bakk, Tekle & Vermunt 2013).
#     Proportional-assignment ML distal fits (dis_gauss_ml, both_ml_prop)
#     change; modal-assignment, BCH, and covariate fits are unchanged.
#   - 2.0 Step 3: covariate designs built with model.matrix(); the intercept
#     is named "(Intercept)" (was "Intercept"). Names only; all numbers
#     unchanged.
#   - 2.0 Step 11b: Step-2 Jacobian J.2 = d theta2 / d theta1 in the
#     corrected (Step-1 uncertainty) variance. Its item-parameter columns
#     were ordered item by item while the Step-1 variance is ordered class by
#     class, and under modal assignment it differentiated the assignment
#     weights as if they were posteriors. Both checked against numerical
#     Jacobians. Corrected ML variances change (covariate SEs +0-12% here;
#     distal barely); estimates, robust/BCH variances, and fit statistics are
#     unchanged.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-v1-reference.R")

data_list <- v1_data()
ref <- lapply(names(v1_configs), function(nm) {
  message("Fitting ", nm)
  v1_extract(v1_fit(v1_configs[[nm]], data_list))
})
names(ref) <- names(v1_configs)
attr(ref, "tseLCA_version") <- as.character(utils::packageVersion("tseLCA"))
attr(ref, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

saveRDS(ref, "tests/testthat/fixtures/v1_reference.rds", version = 2)
message("Wrote ", length(ref), " reference fits")
