# tests/testthat/test-structural.R
#
# Step 3 building blocks: tse_covariate(), tse_distal(), tse_twostep().
# Each specification is checked against three_step() with the same
# measurement model (passed as step1), which the v1 regression fixtures pin.

# Every specification here is also covered, more cheaply, by other test files;
# this file cross-checks them exhaustively and is skipped on CRAN for time.
skip_on_cran()

items <- paste0("Y", 1:6)
f_items <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
dl <- v1_data()
fit_m <- function(d, ...) {
  set.seed(1L)
  tse_lca(f_items, data = d, nclass = 3, ...)
}
same_fit <- function(a, b, info) {
  expect_equal(coef(a), coef(b), tolerance = 1e-7, info = info)
  expect_equal(vcov(a), vcov(b), tolerance = 1e-7, ignore_attr = TRUE, info = info)
  expect_equal(logLik(a), logLik(b), tolerance = 1e-7, info = info)
}

test_that("tse_covariate() matches three_step() for every estimator setting", {
  specs <- list(
    ml_modal = list(d = "cov_high", method = "ML", assign = "modal", se = "corrected"),
    ml_prop = list(d = "cov_mid", method = "ML", assign = "proportional", se = "corrected"),
    bch_modal = list(d = "cov_high", method = "BCH", assign = "modal", se = "corrected"),
    bch_prop = list(d = "cov_mid", method = "BCH", assign = "proportional", se = "corrected"),
    ml_robust = list(d = "cov_high", method = "ML", assign = "modal", se = "robust"),
    poly = list(d = "poly", method = "ML", assign = "modal", se = "corrected"),
    sparse = list(d = "sparse", method = "ML", assign = "modal", se = "corrected")
  )
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    d <- dl[[sp$d]]
    m <- fit_m(d)
    fit <- tse_covariate(tse_classify(m, assignment = sp$assign), ~ Zp,
                         method = sp$method, se = sp$se)
    ref <- three_step(d, items, 3L, Zp.names = "Zp", step1 = m,
                      use.bch = sp$method == "BCH",
                      use.modal.assignment = sp$assign == "modal",
                      use.simple.cov = sp$se == "robust")
    same_fit(fit, ref, nm)
    expect_equal(fit$estimator, sp$method, info = nm)
    expect_equal(coef(fit, step = "two_step"), coef(ref, step = "two_step"), info = nm)
  }
})

test_that("FIML measurement models carry over to the covariate model", {
  d <- dl$miss
  m <- fit_m(d, missing = "fiml")
  fit <- tse_covariate(tse_classify(m), ~ Zp)
  ref <- three_step(d, items, 3L, Zp.names = "Zp", step1 = m, incomplete = TRUE)
  same_fit(fit, ref, "fiml")
})

test_that("the reference class matches three_step(rebase =) and relevel()", {
  d <- dl$cov_high
  m <- fit_m(d)
  cl <- tse_classify(m)
  fit1 <- tse_covariate(cl, ~ Zp)
  for (r in c("C2", "C3")) {
    fit <- tse_covariate(cl, ~ Zp, ref = r)
    same_fit(fit, three_step(d, items, 3L, Zp.names = "Zp", step1 = m, rebase = r), r)
    same_fit(stats::relevel(fit1, ref = r), fit, paste("relevel", r))
    expect_equal(fit$ref, r)
  }
  # a change of reference class is a linear reparameterization
  g1 <- coef(fit1, matrix = TRUE)
  g2 <- coef(tse_covariate(cl, ~ Zp, ref = "C2"), matrix = TRUE)
  expect_equal(unname(g2[, "C1"]), unname(-g1[, "C2"]), tolerance = 1e-4)
  expect_equal(unname(g2[, "C3"]), unname(g1[, "C3"] - g1[, "C2"]), tolerance = 1e-4)
})

test_that("a measurement model from another sample (multi-sample)", {
  big <- dl$cov_high
  small <- dl$cov_mid[1:250, ]
  m <- fit_m(big)
  fit <- tse_covariate(tse_classify(m, newdata = small), ~ Zp)
  same_fit(fit, three_step(small, items, 3L, Zp.names = "Zp", step1 = m), "multi-sample")
  expect_equal(nobs(fit), 250L)
})

test_that("the uncorrected estimator is a weighted logit of the assigned classes", {
  skip_if_not_installed("nnet")
  d <- dl$cov_mid
  m <- fit_m(d)
  cl <- tse_classify(m)
  fit <- tse_covariate(cl, ~ Zp, method = "none")
  mn <- nnet::multinom(factor(classes(cl)) ~ Zp, data = d, trace = FALSE,
                       reltol = 1e-14, abstol = 1e-14, maxit = 5000)
  expect_equal(unname(coef(fit, matrix = TRUE)), unname(t(coef(mn))), tolerance = 1e-6)
  expect_equal(fit$estimator, "uncorrected")
  expect_equal(fit$se, "robust")
  # the correction moves the estimates away from zero
  ml <- tse_covariate(cl, ~ Zp)
  expect_true(all(abs(coef(ml)[c(2, 4)]) > abs(coef(fit)[c(2, 4)])))
})

test_that("tse_distal() matches three_step() for each family and estimator", {
  d <- dl$dis
  m <- fit_m(d)
  specs <- list(
    gauss_ml_prop = list(z = "Zo", family = "gaussian", method = "ML", assign = "proportional"),
    gauss_bch = list(z = "Zo", family = "gaussian", method = "BCH", assign = "modal"),
    poisson = list(z = "Zpois", family = "poisson", method = "ML", assign = "modal"),
    binomial = list(z = "Zbin", family = "binomial", method = "ML", assign = "modal"),
    multinomial = list(z = "Zcat", family = "multinomial", method = "ML", assign = "modal")
  )
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    fit <- tse_distal(tse_classify(m, assignment = sp$assign),
                      stats::as.formula(paste(sp$z, "~ 1")),
                      family = sp$family, method = sp$method)
    ref <- three_step(d, items, 3L, Zo.name = sp$z, family = sp$family, step1 = m,
                      use.bch = sp$method == "BCH",
                      use.modal.assignment = sp$assign == "modal")
    expect_s3_class(fit, "tseLCA_distal")
    same_fit(fit, ref, nm)
    expect_equal(fit$outcome, sp$z)
    if (sp$family == "gaussian") expect_equal(fit$sigma2, ref$sigma2, info = nm)
  }
  # family objects
  f1 <- tse_distal(tse_classify(m), Zpois ~ 1, family = stats::poisson())
  f2 <- tse_distal(tse_classify(m), Zpois ~ 1, family = "poisson")
  expect_equal(coef(f1), coef(f2))
})

test_that("uncorrected distal means are the class means of the assigned classes", {
  d <- dl$dis
  m <- fit_m(d)
  cl <- tse_classify(m)
  fit <- tse_distal(cl, Zo ~ 1, method = "none")
  expect_equal(unname(coef(fit)), as.vector(tapply(d$Zo, classes(cl), mean)), tolerance = 1e-6)
})

test_that("a covariate model passed to tse_distal() gives the combined model", {
  d <- dl$both
  m <- fit_m(d)
  for (a in c("modal", "proportional")) {
    fc <- tse_covariate(tse_classify(m, assignment = a), ~ Zp)
    fb <- tse_distal(fc, Zo ~ 1)
    ref <- three_step(d, items, 3L, Zp.names = "Zp", Zo.name = "Zo", step1 = m,
                      use.modal.assignment = a == "modal")
    expect_s3_class(fb, "tseLCA_both")
    same_fit(fb, ref, a)
    expect_equal(coef(fb, component = "covariate"), coef(fc))
  }
  fc <- tse_covariate(tse_classify(m), ~ Zp, method = "BCH")
  expect_equal(tse_distal(fc, Zo ~ 1)$estimator, "BCH")
  expect_error(tse_distal(fc, Zo ~ 1, method = "ML"), "one estimator")
})

test_that("tse_twostep() gives the two-step estimates", {
  d <- dl$cov_high
  m <- fit_m(d)
  ft <- tse_twostep(m, ~ Zp)
  ref <- three_step(d, items, 3L, Zp.names = "Zp", step1 = m)
  expect_s3_class(ft, c("tseLCA_twostep", "tseLCA_covariate"))
  expect_equal(coef(ft), coef(ref, step = "two_step"))
  expect_true(all(is.na(vcov(ft))))
  expect_output(print(ft), "Two-step latent class model")
  expect_equal(unname(rowSums(predict(ft, newdata = d[1:5, ]))), rep(1, 5))
})

test_that("tse_twostep(se = TRUE) uses multilevLCA's corrected variance", {
  skip_on_cran()
  d <- dl$cov_high
  m <- fit_m(d)
  ft <- suppressWarnings(tse_twostep(m, ~ Zp, se = TRUE))
  V <- vcov(ft)
  expect_false(anyNA(V))
  expect_true(all(diag(V) > 0))
  expect_equal(coef(ft), coef(tse_twostep(m, ~ Zp)), tolerance = 0.05)

  # another reference class re-parameterizes the estimates and their variance
  f3 <- expect_no_warning(tse_twostep(m, ~ Zp, ref = 3, se = TRUE))
  expect_equal(coef(f3), coef(tse_twostep(m, ~ Zp, ref = 3)), tolerance = 0.05)
  b <- coef(ft)
  expect_equal(unname(coef(f3)), unname(c(-b[3:4], b[1:2] - b[3:4])), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(f3)))[1:2]), unname(sqrt(diag(V))[3:4]), tolerance = 1e-8)
  expect_equal(unname(vcov(f3)[4, 4]), V[2, 2] + V[4, 4] - 2 * V[2, 4], tolerance = 1e-8)
})

test_that("predict(), anova(), and omnibus_test() on Step-3 models", {
  d <- dl$both
  m <- fit_m(d)
  fc <- tse_covariate(tse_classify(m), ~ Zp)
  p <- predict(fc, newdata = data.frame(Zp = c(1, 3, 5, NA)))
  expect_equal(dim(p), c(4L, 3L))
  expect_equal(colnames(p), c("C1", "C2", "C3"))
  expect_equal(unname(rowSums(p[1:3, ])), rep(1, 3))
  expect_true(all(is.na(p[4, ])))
  g <- coef(fc, matrix = TRUE)
  eta <- c(0, g[1, ] + 3 * g[2, ])
  expect_equal(unname(p[2, ]), unname(exp(eta) / sum(exp(eta))))
  expect_equal(predict(fc, newdata = d[1:5, ], type = "class"),
               max.col(predict(fc, newdata = d[1:5, ])))
  # the prediction honors a non-default reference class
  p3 <- predict(tse_covariate(tse_classify(m), ~ Zp, ref = "C3"), newdata = d[1:5, ])
  expect_equal(p3, predict(fc, newdata = d[1:5, ]), tolerance = 1e-4)

  a <- anova(fc)
  expect_s3_class(a, "anova")
  b <- coef(fc)[c("Zp:C2", "Zp:C3")]
  V <- vcov(fc)[names(b), names(b)]
  expect_equal(a["Zp", "Chisq"], as.numeric(t(b) %*% solve(V, b)))
  expect_equal(a["Zp", "Df"], 2)

  fb <- tse_distal(fc, Zo ~ 1)
  ot <- omnibus_test(fb)
  expect_s3_class(ot, "htest")
  expect_equal(unname(ot$parameter), 2)
  expect_output(print(ot), "Wald test")
})

test_that("factor covariates and interactions in tse_covariate()", {
  d <- dl$cov_high
  set.seed(9L)
  d$g <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  m <- fit_m(d)
  fit <- tse_covariate(tse_classify(m), ~ Zp * g)
  expect_equal(rownames(coef(fit, matrix = TRUE)), c("(Intercept)", "Zp", "gb", "Zp:gb"))
  a <- anova(fit)
  expect_equal(rownames(a), c("Zp", "g", "Zp:g"))
  expect_equal(unname(a$Df), c(2, 2, 2))
  expect_equal(dim(predict(fit, newdata = data.frame(Zp = 1:2, g = c("a", "b")))), c(2L, 3L))
})

test_that("input checks", {
  d <- dl$both
  m <- fit_m(d)
  cl <- tse_classify(m)
  expect_error(tse_covariate(m, ~ Zp), "tse_classify")
  expect_error(tse_covariate(cl, Zo ~ Zp), "one-sided")
  expect_error(tse_covariate(cl, ~ Zp, start = matrix(0, 3, 2)), "2 x 2 matrix")
  expect_error(tse_distal(m, Zo ~ 1), "tse_classify")
  expect_error(tse_distal(cl, Zo ~ Zp), "tse_covariate")
  expect_error(tse_distal(cl, ~ Zo), "outcome ~ 1")
  expect_error(tse_distal(cl, Zo ~ 1, family = "gamma"), "must be one of")
  expect_error(tse_distal(cl, Zo ~ 1, family = stats::poisson("identity")), "family objects")
  expect_error(tse_twostep(cl, ~ Zp), "measurement model")
  # user starting values reach the same optimum
  fit <- tse_covariate(cl, ~ Zp)
  fit_s <- tse_covariate(cl, ~ Zp, start = matrix(0, 2, 2))
  expect_equal(coef(fit_s), coef(fit), tolerance = 1e-4)
  expect_null(fit_s$two_step)
})
