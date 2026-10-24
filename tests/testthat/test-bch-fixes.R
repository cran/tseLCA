# tests/testthat/test-bch-fixes.R
#
# Regression tests for two bugs in the BCH third step of tseLCA 2.0.0:
#   1. For Poisson and Binomial distal outcomes, the inverse Hessian used for
#      the standard errors was the Gaussian one (a within-class variance
#      divided by the BCH class weights), not the one of the family.
#   2. The safeguard of the BCH covariate Newton-Raphson loop checked that the
#      Hessian, which is negative definite, was positive semi-definite, and
#      therefore stopped whenever the loop ran past a fifth of the maximum
#      number of iterations.

d <- generate_data(2000L, "high", "covariate", seed = 7L)
m <- local({
  set.seed(1L)
  tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
})
cl <- tse_classify(m)

test_that("BCH standard errors of Poisson and Binomial distal outcomes follow the family", {
  set.seed(11L)
  x <- classes(cl)
  d$cnt <- stats::rpois(nrow(d), exp(c(0.5, 1, 1.5)[x]))
  d$bin <- stats::rbinom(nrow(d), 1, c(0.2, 0.5, 0.8)[x])
  cl2 <- tse_classify(m, assignment = "modal")
  cl2$data <- d
  for (fam in c("poisson", "binomial")) {
    y <- if (fam == "poisson") "cnt" else "bin"
    bch <- suppressWarnings(tse_distal(cl2, stats::reformulate("1", y),
                                       family = fam, method = "BCH"))
    ml <- suppressWarnings(tse_distal(cl2, stats::reformulate("1", y),
                                      family = fam, method = "ML", se = "robust"))
    se_b <- sqrt(diag(vcov(bch)))[seq_len(3)]
    se_m <- sqrt(diag(vcov(ml)))[seq_len(3)]
    # same data and a well-separated model: the two sets of standard errors agree
    expect_equal(unname(se_b), unname(se_m), tolerance = 0.25, info = fam)
  }
})

test_that("the BCH covariate safeguard does not stop a well-defined fit", {
  ctl <- tse_control(step3.maxit = 5L)
  fit <- tryCatch(
    suppressWarnings(tse_covariate(cl, ~ Zp, method = "BCH", control = ctl)),
    error = function(e) e
  )
  expect_false(inherits(fit, "error"), info = if (inherits(fit, "error")) conditionMessage(fit))
})
