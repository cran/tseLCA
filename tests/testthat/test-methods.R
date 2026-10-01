# tests/testthat/test-methods.R
#
# Every S3 method on every kind of fitted object (the v1 fixture
# configurations in helper-v1-reference.R): consistent names and dimensions,
# and the generic stats tooling (logLik/AIC/BIC/nobs/confint/printCoefmat)
# working out of the box.

# Every specification here is also covered, more cheaply, by other test files;
# this file cross-checks them exhaustively and is skipped on CRAN for time.
skip_on_cran()

fits <- local({
  dl <- v1_data()
  lapply(v1_configs, v1_fit, data_list = dl)
})

test_that("the parent class carries shared methods", {
  m <- attr(methods(class = "tseLCA"), "info")$generic
  for (g in c("logLik", "nobs", "plot", "posterior", "classes",
              "class_sizes", "item_probs")) {
    expect_true(g %in% m, info = g)
  }
  s <- attr(methods(class = "tseLCA_structural"), "info")$generic
  for (g in c("coef", "vcov", "summary", "print")) {
    expect_true(g %in% s, info = g)
  }
})

test_that("class hierarchy", {
  expect_equal(class(fits$meas_high), c("tseLCA_measurement", "tseLCA"))
  expect_equal(class(fits$cov_ml_modal),
               c("tseLCA_covariate", "tseLCA_structural", "tseLCA"))
  expect_equal(class(fits$dis_gauss_ml),
               c("tseLCA_distal", "tseLCA_structural", "tseLCA"))
  expect_equal(class(fits$both_ml_prop),
               c("tseLCA_both", "tseLCA_structural", "tseLCA"))
})

for (nm in names(fits)) {
  test_that(paste("generic methods work:", nm), {
    fit <- fits[[nm]]
    fp <- if (inherits(fit, "tseLCA_both")) fit$distal else fit

    # logLik / AIC / BIC / nobs agree with the stored fit statistics
    ll <- logLik(fit)
    expect_s3_class(ll, "logLik")
    expect_equal(as.numeric(ll), fp$llik)
    expect_equal(AIC(fit), fp$AIC)
    expect_equal(BIC(fit), fp$BIC)
    expect_equal(nobs(fit), fp$nobs)

    # coef names match vcov; confint is non-empty and equals est +/- z * SE
    est <- coef(fit)
    V <- vcov(fit)
    expect_equal(names(est), rownames(V))
    expect_equal(rownames(V), colnames(V))
    ci <- confint(fit)
    expect_equal(dim(ci), c(length(est), 2L))
    expect_equal(rownames(ci), names(est))
    se <- sqrt(diag(V))
    expect_equal(unname(ci[, 1]), unname(est - qnorm(0.975) * se))
    expect_equal(unname(ci[, 2]), unname(est + qnorm(0.975) * se))

    # posteriors / classes
    P <- posterior(fit)
    expect_equal(ncol(P), 3L)
    expect_equal(length(classes(fit)), nrow(P))
    expect_true(all(classes(fit) %in% 1:3))

    # measurement-model accessors are available on every fit
    expect_equal(sum(class_sizes(fit)), 1, tolerance = 1e-6)
    expect_equal(ncol(item_probs(fit)), 3L)

    # print / summary / plot run
    expect_output(print(fit))
    expect_output(print(summary(fit)))
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(suppressWarnings(plot(fit)))
  })
}

test_that("coef(summary()) is printCoefmat-ready on structural fits", {
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (!inherits(fit, "tseLCA_structural")) next
    cm <- coef(summary(fit))
    expect_equal(colnames(cm), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    expect_equal(unname(cm[, "Estimate"]), unname(coef(fit)))
    expect_output(printCoefmat(cm, has.Pvalue = TRUE), "Estimate")
  }
})
