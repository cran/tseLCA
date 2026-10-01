# tests/testthat/test-classify.R
#
# Step 2 building block: tse_classify().

items <- paste0("Y", 1:6)
f_items <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
d <- generate_data(500L, "high", "covariate", seed = 41L)
m <- local({
  set.seed(1L)
  tse_lca(f_items, data = d, nclass = 3)
})

test_that("tse_classify() reproduces three_step()'s Step 2", {
  cl <- tse_classify(m)
  expect_s3_class(cl, c("tseLCA_classify", "tseLCA"))
  s1 <- m$measurement_model
  s2 <- lca_step2(s1$Y.exp, s1$fit0, 3L, TRUE, 1e-2, FALSE,
                  ivItemcat = s1$ivItemcat, mDesign = s1$mDesign.exp)
  expect_equal(unname(cl$D), unname(t(s2$p.wx_mat)))
  expect_equal(unname(cl$weights), unname(s2$w.is))
  f <- three_step(d, items, 3L, Zp.names = "Zp", step1 = m, use.simple.cov = TRUE)
  expect_equal(posterior(cl), posterior(f), ignore_attr = TRUE)
  expect_equal(classes(cl), classes(f))
})

test_that("assignment weights and classification errors are well formed", {
  for (a in c("modal", "proportional")) {
    cl <- tse_classify(m, assignment = a)
    expect_equal(cl$assignment, a)
    expect_equal(unname(rowSums(cl$D)), rep(1, 3), info = a)
    expect_true(all(cl$D >= 0 & cl$D <= 1), info = a)
    expect_equal(rowSums(cl$weights), rep(1, nrow(d)), info = a)
    expect_equal(dimnames(cl$D), list(paste0("X=C", 1:3), paste0("W=C", 1:3)))
  }
  modal <- tse_classify(m)
  expect_true(all(modal$weights %in% c(0, 1)))
  expect_equal(max.col(modal$weights), classes(modal))
  prop <- tse_classify(m, assignment = "proportional")
  expect_equal(unname(prop$weights), unname(posterior(prop)))
  # proportional assignment spreads more weight off the diagonal
  expect_true(all(diag(prop$D) < diag(modal$D)))
  p <- posterior(modal)
  expect_equal(modal$entropy.R2, 1 - sum(-p * log(p)) / (nrow(p) * log(3)))
})

test_that("newdata is classified with the fixed measurement model", {
  new <- d[1:200, ]
  cl <- tse_classify(m, newdata = new)
  expect_equal(nobs(cl), 200L)
  expect_equal(posterior(cl), predict(m, newdata = new), ignore_attr = TRUE)
  expect_identical(cl$data, new)
  expect_equal(class_sizes(cl), class_sizes(m))
  expect_false(isTRUE(all.equal(cl$D, tse_classify(m)$D)))
})

test_that("missing-data handling is inherited from the measurement model", {
  dm <- d
  dm$Y1[1:40] <- NA
  set.seed(1L)
  mf <- tse_lca(f_items, data = dm, nclass = 3, missing = "fiml")
  expect_equal(nobs(tse_classify(mf)), nrow(dm))
  set.seed(1L)
  ml <- tse_lca(f_items, data = dm, nclass = 3)
  cl <- tse_classify(ml)
  expect_equal(nobs(cl), nrow(dm) - 40L)
  expect_equal(cl$rows, 41:nrow(dm))
})

test_that("measurement models from three_step() need newdata", {
  set.seed(1L)
  t3 <- three_step(d, items, 3L)
  expect_error(tse_classify(t3), "does not store its data")
  cl <- tse_classify(t3, newdata = d)
  expect_equal(posterior(cl), posterior(tse_classify(m)), ignore_attr = TRUE)
})

test_that("methods and input checks", {
  cl <- tse_classify(m)
  expect_output(print(cl), "Classification error probabilities")
  expect_error(logLik(cl), "no log-likelihood")
  expect_equal(nrow(posterior(cl)), nrow(d))
  expect_equal(item_probs(cl), item_probs(m))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(suppressWarnings(plot(cl)))

  expect_error(tse_classify(cl), "measurement model from tse_lca")
  set.seed(1L)
  expect_error(tse_classify(tse_lca(f_items, d, 1)), "at least two classes")
  expect_error(tse_classify(m, newdata = as.matrix(d)), "data frame")
  expect_error(tse_classify(m, assignment = "random"))
})

test_that("Step-3 models accept the classified data with added columns", {
  cl <- tse_classify(m)
  d2 <- d
  d2$group <- factor(ifelse(d2$Zp > 3, "high", "low"))
  expect_error(tse_covariate(cl, ~ Zp + group), "not found")
  fit <- tse_covariate(cl, ~ Zp + group, data = d2)
  expect_equal(rownames(coef(fit, matrix = TRUE)), c("(Intercept)", "Zp", "grouplow"))
  expect_equal(nrow(predict(fit)), nrow(d2))
  d2$Zo2 <- rnorm(nrow(d2))
  expect_s3_class(tse_distal(cl, Zo2 ~ 1, data = d2), "tseLCA_distal")
  expect_s3_class(tse_distal(fit, Zo2 ~ 1, data = d2), "tseLCA_both")
  # other rows are rejected
  expect_error(tse_covariate(cl, ~ Zp, data = d2[-1, ]), "same rows")
  bad <- d2
  bad$Y1 <- 1 - bad$Y1
  expect_error(tse_covariate(cl, ~ Zp, data = bad), "same rows")
})
