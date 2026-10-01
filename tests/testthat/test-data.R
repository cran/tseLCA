# tests/testthat/test-data.R
#
# Input handling (tseLCA 2.0): indicator recoding, model-matrix covariate
# designs, distal outcome validation, and tse_control().

items <- paste0("Y", 1:6)
d0 <- generate_data(400L, "high", "covariate", seed = 21L)

fit_meas <- function(d, ...) {
  set.seed(1L)
  three_step(d, items, 3L, ...)
}

# -- indicators -------------------------------------------------------------------

test_that("indicators coded 1..K, as factors, characters, or logicals fit identically", {
  ref <- fit_meas(d0)
  variants <- list(
    one_based = within(d0, for (v in items) assign(v, get(v) + 1L)),
    factor = within(d0, for (v in items) assign(v, factor(get(v), labels = c("no", "yes")))),
    character = within(d0, for (v in items) assign(v, c("no", "yes")[get(v) + 1L])),
    logical = within(d0, for (v in items) assign(v, get(v) == 1L))
  )
  for (nm in names(variants)) {
    f <- fit_meas(variants[[nm]][, c(items, "Zp")])
    expect_equal(logLik(f), logLik(ref), info = nm)
    expect_equal(class_sizes(f), class_sizes(ref), info = nm)
    expect_equal(posterior(f), posterior(ref), info = nm)
  }
})

test_that("indicator categories are stored with the measurement model", {
  d <- d0
  d$Y1 <- factor(d$Y1, labels = c("no", "yes"))
  f <- fit_meas(d)
  expect_equal(f$measurement_model$Y.levels$Y1, c("no", "yes"))
  expect_equal(f$measurement_model$Y.levels$Y2, c(0, 1))
})

test_that("a reused measurement model applies its categories to new data", {
  d <- d0
  for (v in items) d[[v]] <- factor(d[[v]], labels = c("no", "yes"))
  m <- fit_meas(d)
  # new sample in which Y1 only takes one of its two categories
  new <- d[d$Y1 == "yes", ]
  f <- three_step(new, items, 3L, Zp.names = "Zp", step1 = m, use.simple.cov = TRUE)
  expect_equal(unname(f$measurement_model$ivItemcat), rep(2L, 6))
  expect_equal(nrow(posterior(f)), nrow(new))
  # values outside the stored categories are an error
  bad <- d
  bad$Y2 <- factor(ifelse(bad$Y2 == "yes", "maybe", "no"))
  expect_error(
    three_step(bad, items, 3L, Zp.names = "Zp", step1 = m),
    "outside its categories"
  )
})

test_that("indicators need at least two categories", {
  d <- d0
  d$Y3 <- 1L
  expect_error(fit_meas(d), "fewer than two observed categories")
})

# -- covariates -------------------------------------------------------------------

test_that("factor covariates are dummy coded like explicit indicator columns", {
  d <- d0
  set.seed(5L)
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$gb <- as.numeric(d$g == "b")
  d$gc <- as.numeric(d$g == "c")
  m <- fit_meas(d)
  f_fac <- three_step(d, items, 3L, Zp.names = c("Zp", "g"), step1 = m)
  f_dum <- three_step(d, items, 3L, Zp.names = c("Zp", "gb", "gc"), step1 = m)
  expect_equal(unname(coef(f_fac)), unname(coef(f_dum)))
  expect_equal(unname(vcov(f_fac)), unname(vcov(f_dum)))
  expect_equal(rownames(coef(f_fac, matrix = TRUE)), c("(Intercept)", "Zp", "gb", "gc"))
})

test_that("covariate formulas support interactions and transformations", {
  d <- d0
  set.seed(6L)
  d$g <- factor(sample(c("a", "b"), nrow(d), replace = TRUE))
  cd <- clean_data(d, items, Zp.formula = ~ log(Zp) * g)
  expect_equal(colnames(cd$Z_mat), c("(Intercept)", "log(Zp)", "gb", "log(Zp):gb"))
  expect_equal(unname(cd$Z_mat[, "log(Zp):gb"]), log(d$Zp) * (d$g == "b"))
  expect_error(clean_data(d, items, Zp.formula = y ~ Zp), "one-sided")
  expect_error(clean_data(d, items, Zp.names = "nope"), "not found")
})

test_that("unused factor levels are dropped from the covariate design", {
  d <- d0
  d$g <- factor(rep(c("a", "b"), length.out = nrow(d)), levels = c("a", "b", "c"))
  cd <- clean_data(d, items, Zp.names = "g")
  expect_equal(colnames(cd$Z_mat), c("(Intercept)", "gb"))
})

test_that("rows with missing covariates are excluded from the covariate step only", {
  d <- d0
  d$Zp[1:20] <- NA
  cd <- clean_data(d, items, Zp.names = "Zp")
  expect_equal(nrow(cd$Y.obs), nrow(d))
  expect_equal(nrow(cd$Z_mat), nrow(d) - 20L)
})

test_that("combined models use distal rows with complete covariates", {
  d <- generate_data(500L, "high", "covariate", seed = 22L)
  d$Zo <- draw_Zo(d$X, bk2018_params$distal_params)
  d$Zp[1:25] <- NA
  set.seed(1L)
  m <- three_step(d, items, 3L)
  f <- three_step(d, items, 3L, Zp.names = "Zp", Zo.name = "Zo", step1 = m,
                  use.simple.cov = TRUE)
  # (1.1.1 kept those rows in the distal model with NA class priors)
  expect_false(anyNA(coef(f)))
  expect_false(anyNA(vcov(f)[1:4, 1:4]))
  expect_equal(nobs(f), nrow(d) - 25L)
})

# -- distal outcomes --------------------------------------------------------------

test_that("distal outcomes are validated for their family", {
  d <- generate_data(300L, "high", "distal", seed = 23L)
  set.seed(1L)
  m <- three_step(d, items, 3L)
  d$Zbin01 <- as.integer(d$Zo > 0)
  d$Zbin_f <- factor(ifelse(d$Zo > 0, "high", "low"), levels = c("low", "high"))
  d$Zbin_l <- d$Zo > 0
  fits <- lapply(c("Zbin01", "Zbin_f", "Zbin_l"), function(z) {
    three_step(d, items, 3L, Zo.name = z, family = "binomial", step1 = m,
               use.simple.cov = TRUE)
  })
  expect_equal(coef(fits[[2]]), coef(fits[[1]]))
  expect_equal(coef(fits[[3]]), coef(fits[[1]]))

  d$Zchr <- as.character(d$X)
  expect_error(three_step(d, items, 3L, Zo.name = "Zchr", step1 = m), "numeric")
  expect_error(
    three_step(d, items, 3L, Zo.name = "Zo", family = "poisson", step1 = m),
    "count"
  )
  expect_error(
    three_step(d, items, 3L, Zo.name = "Zo", family = "binomial", step1 = m),
    "binary"
  )
  expect_error(three_step(d, items, 3L, Zo.name = "nope", step1 = m), "not found")
})

# -- tse_control ------------------------------------------------------------------

test_that("tse_control() validates and maps settings", {
  ctl <- tse_control()
  expect_s3_class(ctl, "tse_control")
  expect_equal(ctl$step1.maxit, 5000L)
  expect_equal(ctl$hessian, "observed")
  expect_output(print(ctl), "step1.maxit")

  opts <- .opts_from_control(tse_control(step3.maxit = 50, hessian = "opg", n_init = 5),
                             use.bch = TRUE)
  expect_equal(opts$em.maxIter, 50L)
  expect_true(opts$correct.spec)
  expect_equal(opts$n_init, 5L)
  expect_true(opts$use.bch)

  expect_error(tse_control(step1.maxit = 0), "positive whole number")
  expect_error(tse_control(step1.restarts = -1), "non-negative")
  expect_error(tse_control(step1.restart.R2 = 2), "between 0 and 1")
  expect_error(tse_control(hessian = "exact"))
  expect_error(.opts_from_control(list()), "tse_control")
})
