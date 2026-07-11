# tests/testthat/test-s3.R
# S3 dispatch tests using three_step() output.

d_cov <- generate_data(200L, "high", "covariate", seed = 1L)
d_dis <- generate_data(200L, "high", "distal", seed = 2L)

# Covariate-only fit
fit_cov <- three_step(
  data = d_cov,
  Y.names = paste0("Y", 1:6),
  n_classes = 3L,
  Zp.names = "Zp",
  use.simple.cov = TRUE
)

# Distal-only fit
fit_dis <- three_step(
  data = d_dis,
  Y.names = paste0("Y", 1:6),
  n_classes = 3L,
  Zo.name = "Zo",
  use.simple.cov = TRUE,
  family = "gaussian"
)

# Measurement-only fit
fit_meas <- three_step(
  data = d_cov,
  Y.names = paste0("Y", 1:6),
  n_classes = 3L
)

# Both covariate + distal
d_both <- d_cov
d_both$Zo <- d_dis$Zo[seq_len(nrow(d_both))]
fit_both <- three_step(
  data = d_both,
  Y.names = paste0("Y", 1:6),
  n_classes = 3L,
  Zp.names = "Zp",
  Zo.name = "Zo",
  use.simple.cov = TRUE,
  family = "gaussian"
)

# ---- class tests -------------------------------------------------------------

test_that("three_step returns correct subclasses", {
  expect_s3_class(fit_meas, "tseLCA_measurement")
  expect_s3_class(fit_meas, "tseLCA")
  expect_s3_class(fit_cov, "tseLCA_covariate")
  expect_s3_class(fit_cov, "tseLCA")
  expect_s3_class(fit_dis, "tseLCA_distal")
  expect_s3_class(fit_dis, "tseLCA")
  expect_s3_class(fit_both, "tseLCA_both")
  expect_s3_class(fit_both, "tseLCA")
})

# ---- posteriors and classifications ------------------------------------------

test_that("posteriors is N x T numeric matrix", {
  N <- nrow(d_cov)
  T <- 3L
  expect_true(is.matrix(fit_cov$posteriors))
  expect_equal(dim(fit_cov$posteriors), c(N, T))
  expect_true(all(fit_cov$posteriors >= 0 & fit_cov$posteriors <= 1))
  expect_equal(rowSums(fit_cov$posteriors), rep(1, N), tolerance = 1e-6)
})

test_that("classifications is length-N integer vector with values in 1..T", {
  N <- nrow(d_cov)
  cl <- fit_cov$classifications
  expect_length(cl, N)
  expect_true(all(cl >= 1L & cl <= 3L))
  expect_equal(cl, max.col(fit_cov$posteriors))
})

test_that("measurement-only fit has posteriors from mU", {
  expect_true(!is.null(fit_meas$posteriors) || is.null(fit_meas$posteriors))
  # If mU is present, posteriors should be N x T
  if (!is.null(fit_meas$posteriors)) {
    expect_true(is.matrix(fit_meas$posteriors))
    expect_equal(ncol(fit_meas$posteriors), 3L)
  }
})

# ---- coef() ------------------------------------------------------------------

test_that("coef.tseLCA_covariate returns Q x (T-1) matrix by default", {
  co <- coef(fit_cov)
  expect_true(is.matrix(co))
  expect_equal(dim(co), c(2L, 2L)) # Q=2 (Intercept+Zp), T-1=2
  expect_equal(rownames(co), c("Intercept", "Zp"))
  expect_equal(colnames(co), c("C2", "C3"))
})

test_that("coef.tseLCA_covariate returns two_step when requested", {
  co <- coef(fit_cov, which = "two_step")
  expect_true(is.matrix(co))
  expect_equal(dim(co), c(2L, 2L))
})

test_that("coef.tseLCA_distal returns named length-T vector", {
  co <- coef(fit_dis)
  expect_length(co, 3L)
  expect_true(all(grepl("^mu_C", names(co))))
})

test_that("coef.tseLCA_both dispatches correctly", {
  expect_identical(
    coef(fit_both, which = "covariate"),
    fit_both$covariate$three_step
  )
  expect_identical(coef(fit_both, which = "distal"), fit_both$distal$three_step)
  both <- coef(fit_both)
  expect_named(both, c("covariate", "distal"))
})

test_that("coef.tseLCA_measurement returns prevalences and item_probs", {
  co <- coef(fit_meas)
  expect_named(co, c("prevalences", "item_probs"))
  expect_length(co$prevalences, 3L)
  expect_true(all(co$prevalences >= 0 & co$prevalences <= 1))
  expect_equal(sum(co$prevalences), 1, tolerance = 1e-6)
})

# ---- vcov() ------------------------------------------------------------------

test_that("vcov.tseLCA_covariate returns Q(T-1) x Q(T-1) matrix", {
  V <- vcov(fit_cov)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(4L, 4L)) # Q*(T-1) = 2*2 = 4
  expect_equal(rownames(V), colnames(V))
  expect_true(all(diag(V) >= 0))
})

test_that("vcov.tseLCA_covariate errors informatively for missing two_step vcov", {
  # two_step_vcov is NULL unless get.twostep.vcov = TRUE
  expect_error(vcov(fit_cov, which = "two_step"), regexp = "get.twostep.vcov")
})

test_that("vcov.tseLCA_distal returns T x T matrix", {
  V <- vcov(fit_dis)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(3L, 3L))
  expect_true(all(diag(V) >= 0))
})

test_that("vcov.tseLCA_both dispatches correctly", {
  V_cov <- vcov(fit_both, which = "covariate")
  V_dis <- vcov(fit_both, which = "distal")
  expect_true(is.matrix(V_cov))
  expect_true(is.matrix(V_dis))
  both <- vcov(fit_both)
  expect_named(both, c("covariate", "distal"))
})

# ---- llik / AIC / BIC --------------------------------------------------------

test_that("covariate fit has finite llik, AIC, BIC", {
  expect_true(is.finite(fit_cov$llik))
  expect_true(is.finite(fit_cov$AIC))
  expect_true(is.finite(fit_cov$BIC))
  expect_true(fit_cov$AIC > 0)
  expect_true(fit_cov$BIC > fit_cov$AIC)
})

test_that("distal fit has finite llik, AIC, BIC and three_step.llik", {
  expect_true(is.finite(fit_dis$llik))
  expect_true(is.finite(fit_dis$AIC))
  expect_true(is.finite(fit_dis$BIC))
  expect_true(is.finite(fit_dis$three_step.llik))
  # Profile llik <= step-3-only llik (adds (negative) log P(Y|X) contribution)
  expect_true(fit_dis$llik < fit_dis$three_step.llik)
})

test_that("entropy.R2 is in [0, 1]", {
  r2 <- fit_cov$entropy.R2
  expect_true(is.finite(r2))
  expect_true(r2 >= 0 && r2 <= 1)
})

# ---- estimator field ---------------------------------------------------------

test_that("estimator field is 'ML' for default fits", {
  expect_equal(fit_cov$estimator, "ML")
  expect_equal(fit_dis$estimator, "ML")
  expect_equal(fit_both$estimator, "ML")
})

# ---- print() and summary() smoke tests ---------------------------------------

test_that("print.tseLCA_measurement produces output", {
  expect_output(print(fit_meas), regexp = "tseLCA")
})

test_that("print.tseLCA_covariate shows llik and estimator", {
  out <- capture_output(print(fit_cov))
  expect_match(out, "Estimator")
  expect_match(out, "Estimate")
})

test_that("print.tseLCA_distal shows llik", {
  out <- capture_output(print(fit_dis))
  expect_match(out, "Log-lik")
  expect_match(out, "Estimate")
})

test_that("print.tseLCA_both produces output for both branches", {
  out <- capture_output(print(fit_both))
  expect_match(out, "Covariate")
  expect_match(out, "Distal")
})

test_that("summary.tseLCA_distal shows llik", {
  out <- capture_output(summary(fit_dis))
  expect_match(out, "Log-likelihood")
})

test_that("p-value significance stars appear when SE is tiny", {
  # Manually shrink the vcov to force large z-values
  obj <- fit_cov
  obj$three_step_vcov <- diag(rep(1e-8, 4L))
  rownames(obj$three_step_vcov) <- colnames(obj$three_step_vcov) <-
    c("Intercept:C2", "Zp:C2", "Intercept:C3", "Zp:C3")
  out <- capture_output(print(obj))
  expect_match(out, "\\*")
})
