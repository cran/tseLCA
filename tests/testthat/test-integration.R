# tests/testthat/test-integration.R
# Full pipeline tests using multilevLCA.
# All tests run on CRAN

# ---- Step 1 ----------------------------------------------------------------------------------------------------------------------------------------

test_that("lca_step1 returns correctly-shaped fit0", {
  d <- generate_data(200L, "high", "distal", seed = 42L)
  s1 <- lca_step1(d, paste0("Y", 1:6), n_classes = 3L, verbose = FALSE)

  fit0 <- s1$fit0
  expect_type(fit0, "list")
  #mPhi: n_free x T = 6 x 3 (one row per binary item)
  expect_equal(dim(fit0$mPhi), c(6L, 3L))
  #vPi: length T, sums to 1
  expect_length(fit0$vPi, 3L)
  expect_equal(sum(fit0$vPi), 1, tolerance = 1e-10)
  #All phi values in (0, 1)
  expect_true(all(fit0$mPhi > 0 & fit0$mPhi < 1))
  #LLKSeries present
  expect_true(nrow(fit0$LLKSeries) >= 1L)
})

test_that("lca_step1 with Zp returns fitZ with named mGamma", {
  d <- generate_data(200L, "high", "covariate", seed = 43L)
  s1 <- lca_step1(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.two.step = TRUE,
    verbose = FALSE
  )

  expect_false(is.null(s1$fitZ))
  #mGamma: Q x (T-1) = 2 x 2
  expect_equal(dim(s1$fitZ$mGamma), c(2L, 2L))
  expect_equal(rownames(s1$fitZ$mGamma)[1L], "Intercept")
  expect_equal(colnames(s1$fitZ$mGamma), c("C2", "C3"))
})

# ---- Measurement only --------------------------------------------------------------------------------------------------------------------

test_that("three_step measurement-only returns tseLCA_measurement", {
  d <- generate_data(200L, "high", "distal", seed = 50L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 2L,
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "tseLCA_measurement")
  expect_s3_class(fit, "tseLCA")
  expect_true(is.finite(fit$AIC))
  expect_true(is.finite(fit$BIC))
  expect_equal(fit$n_classes, 2L)
})

# ---- Covariate model ----------------------------------------------------------------------------------------------------------------------

test_that("three_step covariate returns tseLCA_covariate with correct structure", {
  d <- generate_data(250L, "high", "covariate", seed = 100L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "tseLCA_covariate")
  expect_equal(fit$n_classes, 3L)

  #Coefficient matrix: Q x (T-1) = 2 x 2
  co <- fit$three_step
  expect_equal(dim(co), c(2L, 2L))
  expect_equal(rownames(co)[1L], "Intercept")
  expect_equal(colnames(co), c("C2", "C3"))

  #Vcov: Q(T-1) x Q(T-1) = 4 x 4
  vc <- fit$three_step_vcov
  expect_equal(dim(vc), c(4L, 4L))
  expect_equal(vc, t(vc)) # symmetric
  expect_true(all(sqrt(diag(vc)) > 0))

  #Two-step starting values present
  expect_equal(dim(fit$two_step), c(2L, 2L))

  #Model fit
  expect_true(is.finite(fit$llik))
  expect_true(is.finite(fit$AIC))
  expect_true(is.finite(fit$BIC))
})

test_that("three_step BCH covariate runs and returns finite SEs", {
  d <- generate_data(250L, "high", "covariate", seed = 101L)
  fit <- suppressWarnings(
    three_step(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      Zp.names = "Zp",
      use.bch = TRUE,
      use.simple.cov = TRUE,
      verbose = FALSE
    )
  )
  expect_s3_class(fit, "tseLCA_covariate")
  expect_true(all(is.finite(fit$three_step)))
})

# ---- Distal model ----------------------------------------------------------------------------------------------------------------------------

test_that("three_step gaussian distal returns tseLCA_distal with named estimates", {
  d <- generate_data(250L, "high", "distal", seed = 200L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zo",
    family = "gaussian",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "tseLCA_distal")
  expect_length(fit$three_step, 3L)
  expect_named(fit$three_step, paste0("mu_C", 1:3))
  expect_equal(dim(fit$three_step_vcov), c(3L, 3L))
  expect_equal(rownames(fit$three_step_vcov), paste0("mu_C", 1:3))
  expect_true(all(sqrt(diag(fit$three_step_vcov)) > 0))
  #True mu = (-1, 0, 1) up to class labelling; range should span negatives and positives
  expect_true(min(fit$three_step) < 0)
  expect_true(max(fit$three_step) > 0)
})

test_that("three_step with both Zp and Zo returns tseLCA_both", {
  d <- generate_data(250L, "high", "covariate", seed = 300L)
  #Add a synthetic distal outcome
  d$Zo <- rnorm(nrow(d), mean = d$X - 2, sd = 0.5)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    Zo.name = "Zo",
    family = "gaussian",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "tseLCA_both")
  expect_false(is.null(fit$covariate))
  expect_false(is.null(fit$distal))
  expect_equal(fit$n_classes, 3L)
  expect_equal(fit$family, "gaussian")
})

# ---- Missing data ----------------------------------------------------------------------------------------------------------------------------

test_that("three_step uses all Y rows when Z has missing values", {
  set.seed(1L)
  d <- generate_data(250L, "high", "covariate", seed = 400L)
  #Introduce 20 missing Zp values
  d$Zp[sample(250L, 20L)] <- NA

  fit_full <- three_step(
    generate_data(250L, "high", "covariate", seed = 400L),
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.simple.cov = TRUE,
    verbose = FALSE
  )
  fit_miss <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  #Estimates should be close (same measurement model, ~20 fewer Z obs)
  expect_equal(fit_miss$three_step, fit_full$three_step, tolerance = 0.5)
  #Both converge
  expect_true(all(is.finite(fit_miss$three_step)))
})

# ---- Coverage tests: estimates within 2 SEs of truth ----------------------------------------------------
#
# True DGP parameters (Bakk & Kuha 2018):
#   Covariate: slopes b = (0, -1, 1) for classes (ref, C2, C3)
#              intercepts b0 = (0, 2.3446, -3.6554)
#   Distal:    class means mu = (-1, 1, 0)
#
# Class labels from the estimator may differ from the DGP labelling.
# We align by Zp slope sign: most negative slope -> DGP class 2 (b=-1),
# most positive -> DGP class 3 (b=+1).

test_that("covariate estimates are within 2 SEs of true slopes and intercepts", {
  d <- generate_data(2000L, "high", "covariate", seed = 42L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  #True non-reference class params, sorted by slope ascending: (-1, 1)
  #Align estimated classes to true classes by Zp slope sign
  slopes <- fit$three_step["Zp", ]
  ord <- order(slopes)

  true_intercepts <- c(2.3446, -3.6554) # b0 for (C2, C3) in DGP ordering
  true_slopes <- c(-1, 1)

  ses <- sqrt(diag(fit$three_step_vcov))

  est_int <- fit$three_step["Intercept", ord]
  est_slope <- fit$three_step["Zp", ord]
  se_int <- ses[c(1L, 3L)][ord] # Intercept SEs
  se_slope <- ses[c(2L, 4L)][ord] # Zp SEs

  for (j in 1:2) {
    expect_true(
      abs(est_slope[j] - true_slopes[j]) <= 2 * se_slope[j],
      label = sprintf(
        "slope[%d]: est=%.3f, true=%.3f, 2SE=%.3f",
        j,
        est_slope[j],
        true_slopes[j],
        2 * se_slope[j]
      )
    )
    expect_true(
      abs(est_int[j] - true_intercepts[j]) <= 2 * se_int[j],
      label = sprintf(
        "intercept[%d]: est=%.3f, true=%.3f, 2SE=%.3f",
        j,
        est_int[j],
        true_intercepts[j],
        2 * se_int[j]
      )
    )
  }
})

test_that("covariate slope signs match DGP (negative and positive)", {
  d <- generate_data(300L, "high", "covariate", seed = 77L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  slopes <- fit$three_step["Zp", ]
  #The two estimated slopes should have opposite signs
  expect_true(any(slopes < 0), label = "at least one negative Zp slope")
  expect_true(any(slopes > 0), label = "at least one positive Zp slope")
})

test_that("distal estimates recover true class mean ordering", {
  d <- generate_data(300L, "high", "distal", seed = 88L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zo",
    family = "gaussian",
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  true_mu_sorted <- sort(c(-1, 0, 1))
  est_mu_sorted <- sort(fit$three_step)

  #Each sorted estimate should be within 2 SEs of the sorted truth
  ses_sorted <- sort(sqrt(diag(fit$three_step_vcov)))
  for (j in seq_along(true_mu_sorted)) {
    expect_true(
      abs(est_mu_sorted[j] - true_mu_sorted[j]) <= 2 * ses_sorted[j],
      label = sprintf(
        "mu[%d]: est=%.3f, true=%.3f, se=%.3f",
        j,
        est_mu_sorted[j],
        true_mu_sorted[j],
        ses_sorted[j]
      )
    )
  }
})

test_that("measurement model recovers high-separation phi structure", {
  d <- generate_data(300L, "high", "distal", seed = 55L)
  s1 <- lca_step1(d, paste0("Y", 1:6), n_classes = 3L, verbose = FALSE)
  phi <- s1$fit0$mPhi # 6 x 3

  #With high separation (phi_true = 0.9/0.1), each item should have
  # at least one class with phi > 0.7 and at least one with phi < 0.3
  expect_true(
    all(apply(phi, 1L, max) > 0.65),
    label = "each item has a high-probability class"
  )
  expect_true(
    all(apply(phi, 1L, min) < 0.35),
    label = "each item has a low-probability class"
  )
})

test_that("corrected SEs (use.simple.cov=FALSE) are >= simple SEs", {
  d <- generate_data(300L, "high", "covariate", seed = 500L)
  fit1 <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.modal.assignment = FALSE,
    use.simple.cov = TRUE,
    verbose = FALSE
  )
  fit2 <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    use.modal.assignment = FALSE,
    use.simple.cov = FALSE,
    verbose = FALSE
  )
  se1 <- sqrt(diag(fit1$three_step_vcov))
  se2 <- sqrt(diag(fit2$three_step_vcov))
  expect_true(all(se2 >= se1))
})
