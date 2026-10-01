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
  expect_true(is.matrix(posterior(fit_cov)))
  expect_equal(dim(posterior(fit_cov)), c(N, T))
  expect_true(all(posterior(fit_cov) >= 0 & posterior(fit_cov) <= 1))
  expect_equal(rowSums(posterior(fit_cov)), rep(1, N), tolerance = 1e-6)
})

test_that("classifications is length-N integer vector with values in 1..T", {
  N <- nrow(d_cov)
  cl <- classes(fit_cov)
  expect_length(cl, N)
  expect_true(all(cl >= 1L & cl <= 3L))
  expect_equal(cl, max.col(posterior(fit_cov)))
})

test_that("measurement-only fit has an N x T posterior matrix", {
  expect_true(is.matrix(posterior(fit_meas)))
  expect_equal(ncol(posterior(fit_meas)), 3L)
})

# ---- coef() ------------------------------------------------------------------

test_that("coef.tseLCA_covariate returns a vector named like vcov()", {
  co <- coef(fit_cov)
  expect_false(is.matrix(co))
  expect_length(co, 4L) # Q=2 (Intercept+Zp) x T-1=2
  expect_equal(names(co), rownames(vcov(fit_cov)))
  expect_equal(names(co), c("(Intercept):C2", "Zp:C2", "(Intercept):C3", "Zp:C3"))
})

test_that("coef(matrix = TRUE) returns the Q x (T-1) matrix", {
  co <- coef(fit_cov, matrix = TRUE)
  expect_true(is.matrix(co))
  expect_equal(dimnames(co), list(c("(Intercept)", "Zp"), c("C2", "C3")))
  expect_equal(as.vector(co), unname(coef(fit_cov)))
})

test_that("coef.tseLCA_covariate returns two-step estimates when requested", {
  co <- coef(fit_cov, step = "two_step")
  expect_length(co, 4L)
  expect_equal(names(co), names(coef(fit_cov)))
  expect_equal(dim(coef(fit_cov, step = "two_step", matrix = TRUE)), c(2L, 2L))
})

test_that("coef.tseLCA_distal returns named length-T vector", {
  co <- coef(fit_dis)
  expect_length(co, 3L)
  expect_true(all(grepl("^mu_C", names(co))))
})

test_that("coef.tseLCA_both selects components", {
  expect_equal(names(coef(fit_both, component = "covariate")), names(coef(fit_cov)))
  expect_equal(unname(coef(fit_both, component = "distal")), unname(coef(distal(fit_both))))
  both <- coef(fit_both)
  expect_length(both, 4L + 3L)
  expect_equal(names(both), rownames(vcov(fit_both)))
  expect_named(coef(fit_both, matrix = TRUE), c("covariate", "distal"))
})

test_that("coef() rejects the tseLCA 1.x `which` argument", {
  expect_error(coef(fit_cov, which = "two_step"), regexp = "replaced in tseLCA 2.0")
  expect_error(vcov(fit_both, which = "distal"), regexp = "replaced in tseLCA 2.0")
})

test_that("coef.tseLCA_measurement returns log-ratio parameters named like vcov()", {
  co <- coef(fit_meas)
  expect_length(co, 2L + 3L * 6L)
  expect_equal(names(co), rownames(vcov(fit_meas)))
  pi <- class_sizes(fit_meas)
  expect_equal(unname(co[1:2]), log(unname(pi[2:3]) / unname(pi[1])))
})

test_that("class_sizes() and item_probs() return the measurement model", {
  pi <- class_sizes(fit_meas)
  expect_named(pi, c("C1", "C2", "C3"))
  expect_equal(sum(pi), 1, tolerance = 1e-6)
  expect_equal(dim(item_probs(fit_meas)), c(6L, 3L))
  expect_length(class_sizes(fit_cov), 3L) # available on structural fits
  expect_equal(dim(item_probs(fit_dis)), c(6L, 3L))
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
  expect_error(vcov(fit_cov, step = "two_step"), regexp = "get.twostep.vcov")
})

test_that("vcov.tseLCA_distal returns T x T matrix", {
  V <- vcov(fit_dis)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(3L, 3L))
  expect_true(all(diag(V) >= 0))
})

test_that("vcov.tseLCA_both returns components and a block matrix", {
  V_cov <- vcov(fit_both, component = "covariate")
  V_dis <- vcov(fit_both, component = "distal")
  expect_equal(dim(V_cov), c(4L, 4L))
  expect_equal(dim(V_dis), c(3L, 3L))
  V <- vcov(fit_both)
  expect_equal(dim(V), c(7L, 7L))
  expect_equal(V[1:4, 1:4], V_cov, ignore_attr = TRUE)
  expect_equal(V[5:7, 5:7], V_dis, ignore_attr = TRUE)
  expect_true(all(is.na(V[1:4, 5:7]))) # cross-covariances not computed
})

# ---- llik / AIC / BIC --------------------------------------------------------

test_that("covariate fit has finite llik, AIC, BIC", {
  expect_true(is.finite(as.numeric(logLik(fit_cov))))
  expect_true(is.finite(AIC(fit_cov)))
  expect_true(is.finite(BIC(fit_cov)))
  expect_true(AIC(fit_cov) > 0)
  expect_true(BIC(fit_cov) > AIC(fit_cov))
})

test_that("distal fit has finite llik, AIC, BIC and three_step.llik", {
  expect_true(is.finite(as.numeric(logLik(fit_dis))))
  expect_true(is.finite(AIC(fit_dis)))
  expect_true(is.finite(BIC(fit_dis)))
  expect_true(is.finite(fit_dis$three_step.llik))
  # Profile llik <= step-3-only llik (adds (negative) log P(Y|X) contribution)
  expect_true(as.numeric(logLik(fit_dis)) < fit_dis$three_step.llik)
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

# ---- print() and summary() -----------------------------------------------------

test_that("print.tseLCA_measurement produces output", {
  expect_output(print(fit_meas), regexp = "measurement model")
})

test_that("print.tseLCA_covariate shows fit and coefficient table", {
  out <- capture_output(print(fit_cov))
  expect_match(out, "Estimator")
  expect_match(out, "Estimate")
  expect_match(out, "Pr(>|z|)", fixed = TRUE)
})

test_that("print.tseLCA_distal shows llik", {
  out <- capture_output(print(fit_dis))
  expect_match(out, "Log-lik")
  expect_match(out, "Estimate")
})

test_that("print.tseLCA_both produces output for both components", {
  out <- capture_output(print(fit_both))
  expect_match(out, "Covariate")
  expect_match(out, "Distal")
})

test_that("summary() returns a summary object with a printCoefmat-ready table", {
  s <- summary(fit_dis)
  expect_s3_class(s, "summary.tseLCA_structural")
  expect_match(capture_output(print(s)), "Log-lik")
  cm <- coef(summary(fit_cov))
  expect_equal(colnames(cm), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(unname(cm[, "Estimate"]), unname(coef(fit_cov)))
  expect_output(printCoefmat(cm), "Estimate")
  expect_equal(nrow(coef(summary(fit_both))), 7L)
  expect_s3_class(summary(fit_meas), "summary.tseLCA_measurement")
  expect_output(print(summary(fit_meas)), "Item-response probabilities")
})

test_that("p-value significance stars appear when SE is tiny", {
  # Manually shrink the vcov to force large z-values
  obj <- fit_cov
  obj$three_step_vcov <- diag(rep(1e-8, 4L))
  rownames(obj$three_step_vcov) <- colnames(obj$three_step_vcov) <-
    c("(Intercept):C2", "Zp:C2", "(Intercept):C3", "Zp:C3")
  out <- capture_output(print(obj))
  expect_match(out, "*", fixed = TRUE)
})
