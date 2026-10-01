# tests/testthat/test-wrapper.R
#
# One-call interface tseLCA() and component accessors.

# Every specification here is also covered, more cheaply, by other test files;
# this file cross-checks them exhaustively and is skipped on CRAN for time.
skip_on_cran()

dl <- v1_data()
f_m <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
chain <- function(d, cov = NULL, outcome = NULL, family = "gaussian", method = "ML",
                  se = "corrected", assignment = "modal", ref = 1, missing = "listwise") {
  set.seed(1L)
  m <- tse_lca(f_m, data = d, nclass = 3, missing = missing)
  if (is.null(cov) && is.null(outcome)) return(m)
  fit <- tse_classify(m, assignment = assignment)
  if (!is.null(cov)) fit <- tse_covariate(fit, cov, method = method, se = se, ref = ref)
  if (!is.null(outcome)) {
    fit <- tse_distal(fit, stats::reformulate("1", response = outcome), family = family,
                      method = method, se = se)
  }
  fit
}
wrap <- function(formula, d, ...) {
  set.seed(1L)
  tseLCA(formula, data = d, nclass = 3, ...)
}
same <- function(a, b, info) {
  expect_equal(class(a), class(b), info = info)
  if (inherits(a, "tseLCA_structural")) {
    expect_equal(coef(a), coef(b), info = info)
    expect_equal(vcov(a), vcov(b), info = info)
  }
  expect_equal(logLik(a), logLik(b), info = info)
}

test_that("tseLCA() equals the step-wise chain", {
  same(wrap(f_m, dl$cov_high), chain(dl$cov_high), "measurement")
  same(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp, dl$cov_mid, method = "BCH",
            assignment = "proportional"),
       chain(dl$cov_mid, ~ Zp, method = "BCH", assignment = "proportional"), "covariate BCH")
  same(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp, dl$cov_high, ref = "C2", se = "robust"),
       chain(dl$cov_high, ~ Zp, ref = "C2", se = "robust"), "covariate ref")
  same(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1 | Zpois, dl$dis, family = "poisson"),
       chain(dl$dis, outcome = "Zpois", family = "poisson"), "distal poisson")
  same(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1 | Zcat, dl$dis, family = "multinomial"),
       chain(dl$dis, outcome = "Zcat", family = "multinomial"), "distal multinomial")
  same(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, dl$both, assignment = "proportional"),
       chain(dl$both, ~ Zp, "Zo", assignment = "proportional"), "combined")
  same(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp, dl$miss, missing = "fiml"),
       chain(dl$miss, ~ Zp, missing = "fiml"), "fiml")
})

test_that("the formula determines the kind of model", {
  d <- dl$both
  expect_s3_class(wrap(f_m, d), "tseLCA_measurement")
  expect_s3_class(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp, d), "tseLCA_covariate")
  expect_s3_class(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1 | Zo, d), "tseLCA_distal")
  expect_s3_class(wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, d), "tseLCA_both")
  fit <- wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ log(Zp) + I(Zp^2) | Zo, d)
  expect_equal(rownames(coef(fit, component = "covariate", matrix = TRUE)),
               c("(Intercept)", "log(Zp)", "I(Zp^2)"))
  expect_equal(fit$call[[1L]], as.name("tseLCA"))
})

test_that("accessors return the step-wise components", {
  d <- dl$both
  fit <- wrap(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, d)
  m <- measurement(fit)
  expect_s3_class(m, "tseLCA_measurement")
  expect_equal(logLik(m), logLik(chain(d)))
  clf <- classification(fit)
  expect_s3_class(clf, "tseLCA_classify")
  expect_equal(posterior(clf), posterior(fit))
  expect_equal(coef(covariate(fit)), coef(fit, component = "covariate"))
  ds <- distal(fit)
  expect_s3_class(ds, "tseLCA_distal")
  expect_equal(coef(ds), coef(fit, component = "distal"))
  expect_equal(vcov(ds), vcov(fit, component = "distal"))
  expect_equal(omnibus_test(ds)$statistic, omnibus_test(fit)$statistic)
  expect_output(print(ds), "distal outcome")

  fc <- covariate(fit)
  expect_identical(covariate(fc), fc)
  expect_identical(measurement(m), m)
  expect_identical(classification(clf), clf)
  expect_s3_class(measurement(clf), "tseLCA_measurement")
  expect_error(distal(fc), "no distal outcome")
  expect_error(covariate(ds), "no covariate")
  set.seed(1L)
  t3 <- three_step(d, paste0("Y", 1:6), 3L, Zp.names = "Zp")
  expect_error(classification(t3), "does not keep its classification")
})

test_that("input checks", {
  d <- dl$both
  expect_error(tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp, d, 2:3), "single number")
  expect_error(tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo | Zp, d, 3), "at most two")
  expect_error(tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | log(Zo), d, 3), "single column")
  expect_error(tseLCA("formula", d, 3), "must be a formula")
  expect_error(tseLCA(Y1 ~ Zp, d, 3), "cbind")
})
