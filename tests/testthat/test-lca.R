# tests/testthat/test-lca.R
#
# Step 1 building block: tse_lca(), class enumeration, prediction.

items <- paste0("Y", 1:6)
f_items <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
d <- generate_data(500L, "high", "covariate", seed = 31L)

sel <- local({
  set.seed(1L)
  tse_lca(f_items, data = d, nclass = 1:4)
})

test_that("tse_lca() reproduces three_step()'s measurement model", {
  set.seed(1L)
  m <- tse_lca(f_items, data = d, nclass = 3)
  set.seed(1L)
  t3 <- three_step(d, items, 3L)
  expect_s3_class(m, "tseLCA_measurement")
  expect_equal(logLik(m), logLik(t3))
  expect_equal(posterior(m), posterior(t3))
  expect_equal(item_probs(m), item_probs(t3))

  dp <- v1_poly_data()
  dp$Y2[1:30] <- NA
  set.seed(1L)
  mp <- tse_lca(f_items, data = dp, nclass = 3, missing = "fiml")
  set.seed(1L)
  tp <- three_step(dp, items, 3L, incomplete = TRUE)
  expect_equal(logLik(mp), logLik(tp))
  expect_equal(nobs(mp), nrow(dp))
})

test_that("tse_lca() passes its control settings to Step 1", {
  set.seed(2L)
  m <- tse_lca(f_items, data = d, nclass = 3, control = tse_control(n_init = 3))
  set.seed(2L)
  t3 <- three_step(d, items, 3L, n_init = 3L)
  expect_equal(logLik(m), logLik(t3))
  expect_equal(m$control$n_init, 3L)
})

test_that("the enumeration table matches the individual fits", {
  expect_s3_class(sel, "tseLCA_select")
  tab <- as.data.frame(sel)
  expect_equal(tab$nclass, 1:4)
  expect_named(tab, c("nclass", "logLik", "npar", "AIC", "BIC", "SABIC",
                      "entropy.R2", "min.class", "nobs"))
  for (k in 1:4) {
    fk <- sel[[k]]
    expect_equal(fk$n_classes, k)
    expect_equal(tab$logLik[k], as.numeric(logLik(fk)))
    expect_equal(tab$AIC[k], AIC(fk))
    expect_equal(tab$BIC[k], BIC(fk))
  }
  n <- nrow(d)
  expect_equal(tab$SABIC, -2 * tab$logLik + tab$npar * log((n + 2) / 24))
  expect_equal(tab$npar, (tab$nclass - 1) + tab$nclass * 6)
  expect_true(is.na(tab$entropy.R2[1]))
  expect_true(all(tab$min.class > 0 & tab$min.class <= 1))
})

test_that("BIC recovers the true number of classes under high separation", {
  expect_equal(best_model(sel, "BIC")$n_classes, 3L)
  expect_equal(best_model(sel)$n_classes, 3L)
  expect_s3_class(best_model(sel, "SABIC"), "tseLCA_measurement")
  expect_error(sel[[7]], "No 7-class model")
  expect_output(print(sel), "smallest value")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(sel))
  expect_invisible(plot(sel, which = "BIC"))
})

test_that("the one-class model is the independence model", {
  m1 <- sel[[1]]
  p <- colMeans(d[, items])
  ll <- sum(colSums(d[, items]) * log(p) + colSums(1 - d[, items]) * log(1 - p))
  expect_equal(as.numeric(logLik(m1)), ll)
  expect_equal(unname(item_probs(m1)[, 1]), unname(p))
  expect_equal(unname(class_sizes(m1)), 1)
  expect_equal(unname(diag(vcov(m1))), unname(1 / (nrow(d) * p * (1 - p))))
  expect_equal(unname(coef(m1)), unname(qlogis(p)))
  expect_true(all(posterior(m1) == 1))
  expect_error(plot(m1), "at least two classes")
})

test_that("predict() gives posteriors for new data with stored categories", {
  m <- sel[[3]]
  expect_equal(predict(m), posterior(m))
  expect_equal(fitted(m), posterior(m))
  new <- d[1:10, items]
  expect_equal(predict(m, newdata = new), posterior(m)[1:10, ], ignore_attr = TRUE)
  expect_equal(predict(m, newdata = new, type = "class"), classes(m)[1:10])

  # missing values: skipped; an all-missing row is NA
  new2 <- new
  new2$Y1[2] <- NA
  new2[3, ] <- NA
  p2 <- predict(m, newdata = new2)
  expect_true(all(is.na(p2[3, ])))
  expect_equal(rowSums(p2[-3, ]), rep(1, 9), ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(p2[2, ], posterior(m)[2, ], check.attributes = FALSE)))

  # factor-coded model applied to factor-coded new data
  df <- d
  for (v in items) df[[v]] <- factor(df[[v]], labels = c("no", "yes"))
  set.seed(1L)
  mf <- tse_lca(f_items, data = df, nclass = 3)
  expect_equal(predict(mf, newdata = df[1:10, ]), posterior(mf)[1:10, ], ignore_attr = TRUE)
  bad <- df[1:3, ]
  bad$Y1 <- factor(c("no", "maybe", "yes"))
  expect_error(predict(mf, newdata = bad), "outside its categories")
})

test_that("a tse_lca() model can be passed to later steps", {
  m <- sel[[3]]
  set.seed(1L)
  t3 <- three_step(d, items, 3L)
  f1 <- three_step(d, items, 3L, Zp.names = "Zp", step1 = m, use.simple.cov = TRUE)
  f2 <- three_step(d, items, 3L, Zp.names = "Zp", step1 = t3, use.simple.cov = TRUE)
  expect_equal(coef(f1), coef(f2))
})

test_that("formula(), update(), and input checks", {
  m <- sel[[3]]
  expect_equal(formula(m), f_items)
  expect_equal(m$call$nclass, 3L)
  set.seed(1L)
  m2 <- update(m, nclass = 2)
  expect_equal(m2$n_classes, 2L)

  expect_error(tse_lca(cbind(Y1, Y2) ~ Zp, d, 2), "no covariates")
  expect_error(tse_lca(Y1 ~ 1, d, 2), "cbind")
  expect_error(tse_lca(cbind(Y1) ~ 1, d, 2), "at least two indicators")
  expect_error(tse_lca(cbind(Y1, Y1, Y2) ~ 1, d, 2), "Duplicated")
  expect_error(tse_lca(cbind(Y1, log(Y2)) ~ 1, d, 2), "column names")
  expect_error(tse_lca(~ Y1, d, 2), "two-sided")
  expect_error(tse_lca(f_items, d, 0), "positive whole")
  expect_error(tse_lca(f_items, as.matrix(d), 2), "data frame")
  expect_error(tse_lca(f_items, d, 2:3, start = d$X), "single `nclass`")
  expect_error(
    tse_lca(f_items, d, 3, start = d$X, control = tse_control(n_init = 2)),
    "not both"
  )
})

test_that("item_probs() of a fitted model can be used as `start`", {
  m <- sel[[3]]
  m2 <- tse_lca(f_items, data = d, nclass = 3, start = item_probs(m))
  expect_equal(item_probs(m2), item_probs(m), tolerance = 1e-5)
  dp <- v1_poly_data()
  set.seed(1)
  mp <- tse_lca(f_items, data = dp, nclass = 3)
  mp2 <- tse_lca(f_items, data = dp, nclass = 3, start = item_probs(mp))
  expect_equal(item_probs(mp2), item_probs(mp), tolerance = 1e-5)
})

test_that("as_tse_lca() rebuilds a measurement model from its parameters", {
  m <- sel[[3]]
  m2 <- as_tse_lca(f_items, data = d, class_sizes = class_sizes(m), item_probs = item_probs(m))
  expect_s3_class(m2, "tseLCA_measurement")
  expect_equal(logLik(m2), logLik(m), tolerance = 1e-6)
  expect_equal(posterior(m2), posterior(m), tolerance = 1e-6)
  expect_equal(vcov(m2), vcov(m), tolerance = 1e-4)
  fc <- tse_covariate(tse_classify(m), ~ Zp)
  fc2 <- tse_covariate(tse_classify(m2), ~ Zp)
  expect_equal(coef(fc2), coef(fc), tolerance = 1e-5)
  expect_equal(vcov(fc2), vcov(fc), tolerance = 1e-5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(m2))

  # polytomous indicators, and input checks
  dp <- v1_poly_data()
  set.seed(1L)
  mp <- tse_lca(f_items, data = dp, nclass = 3)
  mp2 <- as_tse_lca(f_items, dp, class_sizes(mp), item_probs(mp))
  expect_equal(logLik(mp2), logLik(mp), tolerance = 1e-6)
  expect_error(as_tse_lca(f_items, d, c(.5, .5), item_probs(m)), "6 x 2 matrix")
  expect_error(as_tse_lca(f_items, d, 1, item_probs(m)[, 1, drop = FALSE]), "two or more")
  bad <- item_probs(mp)
  bad[1, 1] <- bad[1, 1] + 0.1
  expect_error(as_tse_lca(f_items, dp, class_sizes(mp), bad), "sum to one")
})

test_that("item_probs() and class_sizes() give delta-method standard errors", {
  numeric_se <- function(m) {
    # the same quantities from finite differences of the softmax maps
    b <- coef(m)
    V <- vcov(m)
    f <- function(beta) {
      K <- length(class_sizes(m))
      pi_ <- exp(c(0, beta[seq_len(K - 1)]))
      c(pi_ / sum(pi_))
    }
    J <- sapply(seq_along(b), function(j) {
      e <- replace(numeric(length(b)), j, 1e-6)
      (f(b + e) - f(b - e)) / 2e-6
    })
    sqrt(diag(J %*% V %*% t(J)))
  }
  for (m in list(sel[[3]], local({
    set.seed(1)
    tse_lca(f_items, data = v1_poly_data(), nclass = 3)
  }))) {
    cs <- class_sizes(m, se = TRUE)
    expect_named(cs, c("estimate", "se"))
    expect_equal(cs$estimate, class_sizes(m))
    expect_equal(unname(cs$se), unname(numeric_se(m)), tolerance = 1e-6)
    ip <- item_probs(m, se = TRUE)
    expect_equal(ip$estimate, item_probs(m))
    expect_equal(dim(ip$se), dim(item_probs(m)))
    expect_true(all(ip$se >= 0))
  }

  # binary item: SE of P(Y = 1 | X = t) is p (1 - p) times the SE of its logit
  m <- sel[[3]]
  ip <- item_probs(m, se = TRUE)
  p <- ip$estimate["P(Y1|C)", 2]
  se_logit <- sqrt(vcov(m)["log(P(Y1=1|C2)/P(Y1=0|C2))", "log(P(Y1=1|C2)/P(Y1=0|C2))"])
  expect_equal(unname(ip$se["P(Y1|C)", 2]), unname(p * (1 - p) * se_logit))

  # structural models report those of their measurement model
  fc <- tse_covariate(tse_classify(m), ~ Zp, ref = 2)
  expect_equal(item_probs(fc, se = TRUE), ip)
  expect_equal(class_sizes(fc, se = TRUE), class_sizes(m, se = TRUE))
})

test_that("FIML with missing indicator values uses the default initialization", {
  dp <- d
  dp$Y2[1:40] <- NA
  expect_warning(
    m <- tse_lca(f_items, data = dp, nclass = 3, missing = "fiml",
                 control = tse_control(n_init = 3)),
    "not available"
  )
  expect_s3_class(m, "tseLCA_measurement")
  expect_equal(nobs(m), nrow(dp))
  expect_error(tse_lca(f_items, data = dp, nclass = 3, missing = "fiml",
                       start = item_probs(sel[[3]])), "cannot be used")
})
