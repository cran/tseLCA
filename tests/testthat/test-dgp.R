# tests/testthat/test-dgp.R
# DGP functions -- no external deps, all run on CRAN.

test_that("make_rho returns a 3 x 6 matrix with correct values", {
  rho <- make_rho(0.8)
  expect_equal(dim(rho), c(3L, 6L))
  expect_true(all(rho > 0 & rho < 1))

  # Class 1 (high responders): all items at pi = 0.8
  expect_equal(rho[1, ], rep(0.8, 6), tolerance = 1e-12)

  # Class 3 (low responders): all items at 1 - pi = 0.2
  expect_equal(rho[3, ], rep(0.2, 6), tolerance = 1e-12)

  # Class 2 (mixed): first 3 high, last 3 low
  expect_equal(rho[2, 1:3], rep(0.8, 3), tolerance = 1e-12)
  expect_equal(rho[2, 4:6], rep(0.2, 6)[4:6], tolerance = 1e-12)
})

test_that("make_rho rejects invalid inputs", {
  expect_error(make_rho(0.4))
  expect_error(make_rho(1.0))
  expect_error(make_rho(c(0.7, 0.8)))
})

test_that("mnl_probs rows sum to 1", {
  probs <- mnl_probs(1:5, bk2018_params$covariate_params)
  expect_equal(dim(probs), c(5L, 3L))
  expect_equal(rowSums(probs), rep(1, 5), tolerance = 1e-12)
})

test_that("mnl_probs marginal prevalences are ~1/3 each", {
  probs <- mnl_probs(1:5, bk2018_params$covariate_params)
  expect_equal(colMeans(probs), c(1 / 3, 1 / 3, 1 / 3), tolerance = 1e-4)
})

test_that("generate_data covariate scenario returns correct structure", {
  d <- generate_data(
    n = 200L,
    separation = "high",
    scenario = "covariate",
    seed = 1L
  )
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 200L)
  expect_true(all(paste0("Y", 1:6) %in% names(d)))
  expect_true("Zp" %in% names(d))
  expect_false("Zo" %in% names(d))
  expect_true(all(d$Zp %in% 1:5))
  expect_true(all(as.matrix(d[, paste0("Y", 1:6)]) %in% c(0L, 1L)))
  expect_true(all(d$X %in% 1:3))
})

test_that("generate_data distal scenario returns correct structure", {
  d <- generate_data(
    n = 200L,
    separation = "mid",
    scenario = "distal",
    seed = 2L
  )
  expect_true("Zo" %in% names(d))
  expect_false("Zp" %in% names(d))
  expect_true(is.numeric(d$Zo))
  expect_equal(nrow(d), 200L)
})

test_that("generate_data is reproducible with same seed", {
  d1 <- generate_data(200L, "high", "covariate", seed = 99L)
  d2 <- generate_data(200L, "high", "covariate", seed = 99L)
  expect_identical(d1, d2)
})

test_that("generate_data differs across seeds", {
  d1 <- generate_data(200L, "high", "covariate", seed = 1L)
  d2 <- generate_data(200L, "high", "covariate", seed = 2L)
  expect_false(identical(d1$Y1, d2$Y1))
})

test_that("generate_data rejects invalid n", {
  expect_error(generate_data(-1L, "high", "covariate"))
  expect_error(generate_data(0L, "high", "covariate"))
  expect_error(generate_data(1.5, "high", "covariate"))
})

test_that("generate_data rejects invalid separation and scenario via match.arg", {
  expect_error(generate_data(100L, "bad", "covariate"))
  expect_error(generate_data(100L, "high", "bad"))
})

test_that("bk2018_params has all required fields with correct dimensions", {
  p <- bk2018_params
  # Top-level fields
  expect_true(all(
    c("separation_levels", "covariate_params", "distal_params") %in% names(p)
  ))
  # Covariate params
  expect_true(all(c("b0", "b") %in% names(p$covariate_params)))
  expect_length(p$covariate_params$b0, 3L)
  expect_length(p$covariate_params$b, 3L)
  # Distal params
  expect_true(all(c("mu", "sigma") %in% names(p$distal_params)))
  expect_length(p$distal_params$mu, 3L)
  # separation_levels has low/mid/high
  expect_true(all(c("low", "mid", "high") %in% names(p$separation_levels)))
})
