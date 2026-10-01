# tests/testthat/test-deprecated.R
#
# The tseLCA 1.x interface keeps working, warns once per session when called
# directly, and gives the same results as the 2.0 interface.

items <- paste0("Y", 1:6)
d <- generate_data(300L, "high", "covariate", seed = 71L)

with_deprecation_warnings <- function(code) {
  old <- options(tseLCA.warn.deprecated = TRUE)
  on.exit(options(old), add = TRUE)
  assign("warned", NULL, envir = .tse_state)
  on.exit(assign("warned", NULL, envir = .tse_state), add = TRUE)
  force(code)
}

test_that("three_step() warns once per session", {
  with_deprecation_warnings({
    set.seed(1L)
    expect_warning(three_step(d, items, 3L), "three_step\\(\\) is deprecated")
    set.seed(1L)
    expect_no_warning(three_step(d, items, 3L))
  })
})

test_that("1.x helpers warn when called directly, not when used internally", {
  with_deprecation_warnings({
    set.seed(1L)
    expect_warning(s1 <- lca_step1(d, items, 3L), "lca_step1\\(\\) is deprecated")
    expect_warning(
      fitZ_from_fit0(s1$fit0, d, items, "Zp"),
      "fitZ_from_fit0\\(\\) is deprecated"
    )
    # the new interface uses the same routines internally without warning
    set.seed(1L)
    expect_no_warning({
      m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
      tse_covariate(tse_classify(m), ~ Zp)
      tse_twostep(m, ~ Zp)
    })
  })
})

test_that("the warnings can be switched off", {
  old <- options(tseLCA.warn.deprecated = FALSE)
  on.exit(options(old), add = TRUE)
  assign("warned", NULL, envir = .tse_state)
  set.seed(1L)
  expect_no_warning(three_step(d, items, 3L))
})

test_that("three_step() and the 2.0 interface give the same fit", {
  set.seed(1L)
  old <- three_step(d, items, 3L, Zp.names = "Zp", use.modal.assignment = FALSE)
  set.seed(1L)
  new <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp, data = d, nclass = 3,
                assignment = "proportional")
  expect_equal(coef(old), coef(new))
  expect_equal(vcov(old), vcov(new))
  expect_equal(logLik(old), logLik(new))
})
