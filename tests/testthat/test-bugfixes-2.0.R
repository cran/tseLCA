# tests/testthat/test-bugfixes-2.0.R
#
# Regression tests for bugs in tseLCA 1.1.1 fixed in 2.0:
#   1. Measurement-only posteriors/classifications were not in data-row order
#      (taken from multilevLCA's fit0$mU, which is sorted by response pattern).
#   2. Polytomous items were decoded from fit0$mU as category 0 on the
#      listwise path, making the Step-1 information matrix singular: the
#      measurement vcov and the corrected Step-3 vcov were all NA.
#   3. AIC/BIC parameter counts counted one-hot indicator columns (covariate
#      models) and ignored covariate coefficients (combined models).

dl <- v1_data()

fit_pair <- function(d, ...) {
  set.seed(1L)
  m <- three_step(d, v1_items, 3L, ...)
  set.seed(1L)
  f <- suppressWarnings(three_step(
    d, v1_items, 3L, Zp.names = "Zp", step1 = m, use.simple.cov = TRUE, ...
  ))
  list(m = m, f = f)
}

test_that("measurement-only posteriors are in data-row order", {
  cases <- list(
    binary = list(d = dl$cov_high),
    polytomous = list(d = dl$poly),
    fiml = list(d = dl$miss, incomplete = TRUE)
  )
  for (nm in names(cases)) {
    args <- cases[[nm]]
    p <- do.call(fit_pair, c(list(args$d), args[-1]))
    # Step-2 posteriors (computed independently in lca_step2) for the same rows
    expect_equal(posterior(p$m), posterior(p$f), tolerance = 1e-8,
                 ignore_attr = TRUE, info = nm)
    expect_equal(classes(p$m), classes(p$f), info = nm)
  }
})

test_that("measurement-only posteriors follow a rebased reference class", {
  set.seed(1L)
  m1 <- three_step(dl$cov_high, v1_items, 3L)
  set.seed(1L)
  m2 <- three_step(dl$cov_high, v1_items, 3L, rebase = "C2")
  expect_equal(posterior(m2), posterior(m1)[, c(2L, 1L, 3L)], ignore_attr = TRUE)
})

test_that("modal classes recover the true classes under high separation", {
  set.seed(1L)
  m <- three_step(dl$cov_high, v1_items, 3L)
  tab <- table(classes(m), dl$cov_high$X)
  # agreement up to label switching
  expect_gt(sum(apply(tab, 1L, max)) / nrow(dl$cov_high), 0.85)
})

test_that("polytomous measurement vcov and corrected Step-3 vcov are finite", {
  for (nm in c("poly", "sparse")) {
    set.seed(1L)
    m <- three_step(dl[[nm]], v1_items, 3L)
    expect_false(anyNA(vcov(m)), info = nm)
    set.seed(1L)
    f <- three_step(dl[[nm]], v1_items, 3L, Zp.names = "Zp", step1 = m)
    set.seed(1L)
    fs <- three_step(dl[[nm]], v1_items, 3L, Zp.names = "Zp", step1 = m,
                     use.simple.cov = TRUE)
    V <- vcov(f)
    expect_false(anyNA(V), info = nm)
    # correcting for Step-1 uncertainty cannot shrink the variances
    expect_true(all(diag(V) >= diag(vcov(fs)) - 1e-10), info = nm)
  }
})

test_that("the sparse data set has item probabilities at the boundary", {
  set.seed(1L)
  m <- three_step(dl$sparse, v1_items, 3L)
  expect_lt(min(item_probs(m)), 1e-3)
})

test_that("fit0$mU fallback decodes both multilevLCA polytomous codings", {
  d <- dl$poly
  d_miss <- d
  set.seed(9L)
  d_miss$Y1[sample(nrow(d), 40L)] <- NA
  for (case in list(list(d = d, inc = FALSE), list(d = d_miss, inc = TRUE))) {
    set.seed(1L)
    m <- suppressWarnings(three_step(case$d, v1_items, 3L, incomplete = case$inc))
    s1 <- m$measurement_model
    decoded <- step1_sample(list(fit0 = s1$fit0), s1$ivItemcat)
    # mU is sorted by response pattern: compare the multisets of rows
    key <- function(M) sort(apply(M, 1L, paste, collapse = ""))
    Yd <- decoded$Y.exp
    Ys <- s1$Y.exp
    if (!is.null(decoded$mDesign)) Yd[decoded$mDesign == 0] <- NA
    if (!is.null(s1$mDesign.exp)) Ys[s1$mDesign.exp == 0] <- NA
    expect_equal(key(Yd), key(Ys), info = paste("incomplete =", case$inc))
  }
})

test_that("a Step-1 model reused through `step1` keeps its own sample", {
  set.seed(1L)
  m <- three_step(dl$cov_high, v1_items, 3L)
  expect_equal(nrow(m$measurement_model$Y.exp), nrow(dl$cov_high))
  small <- dl$cov_mid[1:200, ]
  set.seed(1L)
  f <- three_step(small, v1_items, 3L, Zp.names = "Zp", step1 = m)
  expect_equal(nrow(f$measurement_model$Y.exp), nrow(dl$cov_high))
  expect_equal(nrow(posterior(f)), 200L)
  expect_false(anyNA(vcov(f)))
})

test_that("parameter counts behind logLik/AIC/BIC", {
  n_items <- 6L
  cfg <- function(nm) v1_fit(v1_configs[[nm]], dl)
  # measurement: (T-1) class sizes + T * sum(K-1) item parameters
  expect_equal(attr(logLik(cfg("meas_high")), "df"), 2 + 3 * n_items)
  expect_equal(attr(logLik(cfg("meas_poly")), "df"), 2 + 3 * 2 * n_items)
  # covariate: item parameters + Q * (T-1) logit coefficients (Q = 2)
  expect_equal(attr(logLik(cfg("cov_ml_modal")), "df"), 3 * n_items + 2 * 2)
  expect_equal(attr(logLik(cfg("cov_poly")), "df"), 3 * 2 * n_items + 2 * 2)
  # distal: class sizes + item parameters + T class means (+ sigma2 for
  # the gaussian family)
  expect_equal(attr(logLik(cfg("dis_gauss_ml")), "df"), 2 + 3 * n_items + 3 + 1)
  expect_equal(attr(logLik(cfg("dis_poisson")), "df"), 2 + 3 * n_items + 3)
  expect_equal(attr(logLik(cfg("dis_multinomial")), "df"), 2 + 3 * n_items + 3 * 2)
  # combined: covariate coefficients replace class sizes
  expect_equal(attr(logLik(cfg("both_ml_prop")), "df"), 2 * 2 + 3 * n_items + 3 + 1)
  fit <- cfg("cov_ml_modal")
  expect_equal(AIC(fit), -2 * fit$llik + 2 * 22)
  expect_equal(BIC(fit), -2 * fit$llik + 22 * log(nobs(fit)))
})

test_that("covariate models without an intercept can be fitted", {
  # 1.1.1 labeled the coefficient rows c("Intercept", Zp.names) regardless of
  # include.intercept and failed with a dimnames error.
  d <- generate_data(300L, "high", "covariate", seed = 1L)
  set.seed(1L)
  f <- three_step(d, v1_items, 3L, Zp.names = "Zp",
                  include.intercept = FALSE, use.simple.cov = TRUE)
  expect_equal(dimnames(coef(f, matrix = TRUE)), list("Zp", c("C2", "C3")))
  expect_equal(names(coef(f)), rownames(vcov(f)))
})

test_that("a legacy lca_step1() fit with two-step estimates can be rebased", {
  # Exercises normalize_fitZ_names() / permute_fitZ_classes(): the stored
  # two-step estimates are renamed and re-referenced to match the new
  # reference class.
  d <- dl$cov_high
  set.seed(1L)
  s1 <- lca_step1(d, v1_items, 3L, Zp.names = "Zp")
  expect_false(is.null(s1$fitZ))
  fit <- three_step(d, v1_items, 3L, Zp.names = "Zp", step1 = s1, rebase = "C2",
                    use.simple.cov = TRUE)
  ref <- three_step(d, v1_items, 3L, Zp.names = "Zp", step1 = s1, rebase = "C2",
                    use.simple.cov = TRUE, use.two.step = FALSE)
  expect_equal(colnames(coef(fit, step = "two_step", matrix = TRUE)), c("C1", "C3"))
  expect_equal(rownames(coef(fit, step = "two_step", matrix = TRUE)), c("(Intercept)", "Zp"))
  # the rebased two-step values are a valid start: same optimum as a cold start
  expect_equal(coef(fit), coef(ref), tolerance = 1e-4)
})

test_that("corrected SEs exist when a reference category is on the boundary", {
  # Class 1 never gives the first category of Y1, so P(Y1 = 0 | class 1) is on
  # the boundary; the corrected variance used to be all NA in this case.
  set.seed(7)
  n <- 1500
  Zp <- rnorm(n)
  cls <- apply(cbind(1, exp(0.5 * Zp), exp(-0.5 * Zp)), 1, function(p) sample(3, 1, prob = p))
  probs <- list(c(0, .5, .5), c(.8, .1, .1), c(.1, .1, .8))
  d <- data.frame(Zp = Zp)
  for (j in 1:5) {
    pj <- if (j == 1) probs else list(c(.1, .8, .1), c(.8, .1, .1), c(.1, .1, .8))
    d[[paste0("Y", j)]] <- vapply(cls, function(k) sample(0:2, 1, prob = pj[[k]]), integer(1))
  }
  set.seed(1)
  m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5) ~ 1, data = d, nclass = 3,
               control = tse_control(n_init = 5))
  bdry <- which(item_probs(m)[1, ] < 1e-3)
  expect_length(bdry, 1L)
  V1 <- vcov(m)
  expect_false(anyNA(V1))
  expect_equal(V1, t(V1))
  # the log-ratios of Y1 in that class: finite, moderate SEs (were ~3000),
  # and no variance along their common shift, log P(Y1 = 0 | class)
  blk <- grep(sprintf("^log\\(P\\(Y1=[12]\\|C%d\\)", bdry), rownames(V1))
  expect_length(blk, 2L)
  expect_true(all(sqrt(diag(V1)[blk]) < 1))
  expect_equal(sum(V1[blk, blk]), 0, tolerance = 1e-10)
  fc <- tse_covariate(tse_classify(m), ~ Zp)
  fr <- tse_covariate(tse_classify(m), ~ Zp, se = "robust")
  se_c <- sqrt(diag(vcov(fc)))
  expect_true(all(is.finite(se_c)))
  expect_true(all(se_c >= sqrt(diag(vcov(fr))) - 1e-8))
})

test_that("a singular Step-1 information matrix falls back to robust SEs, with a warning", {
  f <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
  set.seed(1)
  cl <- tse_classify(tse_lca(f, generate_data(500, "high", "covariate", seed = 3), 3))
  set.seed(1)
  cld <- tse_classify(tse_lca(f, generate_data(500, "high", "distal", seed = 3), 3))
  local_mocked_bindings(.step1_varmat = function(...) matrix(NA_real_, 2, 2))
  expect_warning(fc <- tse_covariate(cl, ~ Zp), "robust standard errors are reported")
  expect_equal(vcov(fc), vcov(tse_covariate(cl, ~ Zp, se = "robust")))
  expect_identical(fc$se, "robust")
  expect_warning(fd <- tse_distal(cld, Zo ~ 1, family = "gaussian"), "robust standard errors")
  expect_equal(vcov(fd), vcov(tse_distal(cld, Zo ~ 1, family = "gaussian", se = "robust")))
})

test_that("class labels do not depend on the reference class of the covariate model", {
  skip_on_cran()
  d <- generate_data(600, "high", "distal", seed = 2)
  set.seed(3)
  d$Zp <- rnorm(600)
  d$Zm <- factor(sample(c("a", "b", "c"), 600, TRUE))
  set.seed(1)
  m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
  cl <- tse_classify(m, assignment = "proportional")
  fd <- tse_distal(cl, Zo ~ 1)
  b1 <- tse_distal(tse_covariate(cl, ~ Zp), Zo ~ 1)
  bm1 <- tse_distal(tse_covariate(cl, ~ Zp), Zm ~ 1, family = "multinomial")
  for (r in 2:3) {
    fc <- tse_covariate(cl, ~ Zp, ref = r)
    expect_equal(posterior(fc), posterior(cl))
    expect_equal(classes(fc), classes(cl))
    expect_equal(item_probs(fc), item_probs(m))
    expect_equal(class_sizes(fc), class_sizes(m))

    # combined models: distal parameters of class t are those of class t
    br <- tse_distal(fc, Zo ~ 1)
    mu <- grep("^mu_", names(coef(br)), value = TRUE)
    expect_equal(coef(br)[mu], coef(b1)[mu], tolerance = 1e-5)
    expect_equal(vcov(br)[mu, mu], vcov(b1)[mu, mu], tolerance = 1e-3)
    expect_equal(coef(br)[mu], coef(fd)[mu], tolerance = 0.05)
    expect_equal(posterior(br), posterior(cl))
    bmr <- tse_distal(fc, Zm ~ 1, family = "multinomial")
    pk <- grep("^C[0-9]:", names(coef(bmr)), value = TRUE)
    expect_equal(coef(bmr)[pk], coef(bm1)[pk], tolerance = 1e-5)
    expect_equal(vcov(bmr)[pk, pk], vcov(bm1)[pk, pk], tolerance = 1e-4)
  }
})

test_that("the Step-2 Jacobian matches finite differences (modal and proportional)", {
  set.seed(1)
  d <- generate_data(400, "mid", "distal", seed = 5)
  m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
  fit0 <- m$measurement_model$fit0
  iT <- 3L
  Y <- as.matrix(d[paste0("Y", 1:6)])
  Y <- do.call(cbind, lapply(seq_len(ncol(Y)), function(j) cbind(1 - Y[, j], Y[, j])))
  ivI <- rep(2L, 6)
  D1 <- matrix(1L, nrow(Y), ncol(Y))
  # parameters in the order of the Step-1 variance: class-size log-ratios, then
  # class by class, item by item, logit P(Y = 1 | class)
  u0 <- c(log(fit0$vPi[-1] / fit0$vPi[1]), stats::qlogis(as.vector(fit0$mPhi)))
  for (modal in c(FALSE, TRUE)) {
    s2 <- lca_step2(Y, fit0, iT, modal, 1e-2, FALSE, ivI)
    w0 <- s2$w.is
    f <- function(u) {
      p <- exp(c(0, u[1:2]))
      post <- compute_posteriors(Y, D1, c((p / sum(p))[-1], stats::plogis(u[-(1:2)])), ivI, iT)
      w <- if (modal) w0 else post
      pj <- t(w) %*% post
      pwx <- sweep(pj, 2, colSums(pj), "/")
      g <- sweep(log(pwx), 2, log(diag(pwx)), "-")
      g[row(g) != col(g)]
    }
    h <- 1e-5
    J_num <- sapply(seq_along(u0), function(j) {
      e <- replace(numeric(length(u0)), j, h)
      (f(u0 + e) - f(u0 - e)) / (2 * h)
    })
    J <- s2$compute_J_unc(s2$p.xy, Y, D1, s2$theta1, ivI, iT)
    expect_equal(J, J_num, tolerance = 1e-5, ignore_attr = TRUE)
  }
})
