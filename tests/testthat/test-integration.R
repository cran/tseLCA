# tests/testthat/test-integration.R

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

# ---- Step 1: external starting values ---------------------------------------------------------------------------------------------------

test_that("lca_step1_startval fits from a user-supplied classification", {
  d <- generate_data(300L, "high", "distal", seed = 44L)
  s1 <- lca_step1_startval(
    d,
    Y.names = paste0("Y", 1:6),
    n_classes = 3L,
    startval = d$X
  )

  expect_type(s1, "list")
  expect_null(s1$fitZ)
  fit0 <- s1$fit0
  expect_equal(dim(fit0$mPhi), c(6L, 3L))
  expect_equal(sum(fit0$vPi), 1, tolerance = 1e-10)
  expect_true(all(fit0$mPhi > 0 & fit0$mPhi < 1))
})

test_that("lca_step1_startval validates the startval vector", {
  d <- generate_data(100L, "high", "distal", seed = 45L)
  expect_error(
    lca_step1_startval(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      startval = d$X[1:10]
    ),
    "length"
  )
  expect_error(
    lca_step1_startval(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      startval = rep(5L, nrow(d))
    ),
    "between 1 and n_classes"
  )
  expect_error(
    lca_step1_startval(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      startval = c(NA_integer_, d$X[-1])
    ),
    "NA"
  )
})

test_that("lca_step1(startval=) matches lca_step1_startval() and skips restarts", {
  d <- generate_data(300L, "high", "distal", seed = 46L)
  s1 <- lca_step1_startval(
    d,
    Y.names = paste0("Y", 1:6),
    n_classes = 3L,
    startval = d$X
  )
  s2 <- lca_step1(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    startval = d$X,
    verbose = FALSE
  )

  expect_equal(s1$fit0$vPi, s2$fit0$vPi)
  expect_equal(s1$fit0$mPhi, s2$fit0$mPhi)
})

test_that("three_step(startval=) reproduces the same measurement fit", {
  d <- generate_data(300L, "high", "covariate", seed = 47L)

  fit <- three_step(
    d,
    Y.names = paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    startval = d$X,
    use.simple.cov = TRUE
  )
  s1 <- lca_step1_startval(
    d,
    Y.names = paste0("Y", 1:6),
    n_classes = 3L,
    startval = d$X
  )

  expect_s3_class(fit, "tseLCA_covariate")
  expect_equal(fit$measurement_model$fit0$vPi, s1$fit0$vPi)
  expect_equal(dim(fit$three_step), c(2L, 2L))
})

test_that("three_step errors when both step1 and startval are supplied", {
  d <- generate_data(200L, "high", "distal", seed = 48L)
  s1 <- lca_step1_startval(d, paste0("Y", 1:6), n_classes = 3L, startval = d$X)

  expect_error(
    three_step(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      step1 = s1,
      startval = d$X
    ),
    "mutually exclusive"
  )
})

# ---- Step 1: item-response probability matrix as startval -------------------------------------------------------------------------------

phi_from_fit0 <- function(fit0, ivItemcat) {
  # Expand multilevLCA's compact mPhi (dichotomous items: 1 row = P(Y=1|C))
  # into the full one-row-per-category convention `startval` expects.
  blocks <- vector("list", length(ivItemcat))
  row <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    if (K_h == 2L) {
      blocks[[h]] <- rbind(1 - fit0$mPhi[row, ], fit0$mPhi[row, ])
      row <- row + 1L
    } else {
      blocks[[h]] <- fit0$mPhi[row:(row + K_h - 1L), , drop = FALSE]
      row <- row + K_h
    }
  }
  do.call(rbind, blocks)
}

test_that("lca_step1_startval derives a classification from a phi matrix", {
  d <- generate_data(300L, "high", "distal", seed = 51L)
  Y.names <- paste0("Y", 1:6)
  fit_ref <- lca_step1(d, Y.names, n_classes = 3L)$fit0
  phi <- phi_from_fit0(fit_ref, rep(2L, 6L))

  expect_equal(dim(phi), c(12L, 3L))
  expect_equal(colSums(phi[1:2, ]), c(1, 1, 1), ignore_attr = TRUE)

  s1_phi <- lca_step1_startval(d, Y.names, n_classes = 3L, startval = phi)
  expect_null(s1_phi$fitZ)
  expect_equal(dim(s1_phi$fit0$mPhi), c(6L, 3L))
  # Deriving startval from the reference fit's own phi should reproduce it.
  expect_equal(sort(s1_phi$fit0$vPi), sort(fit_ref$vPi), tolerance = 1e-4)
})

test_that("lca_step1(startval=phi matrix) matches lca_step1_startval()", {
  d <- generate_data(300L, "high", "distal", seed = 52L)
  Y.names <- paste0("Y", 1:6)
  fit_ref <- lca_step1(d, Y.names, n_classes = 3L)$fit0
  phi <- phi_from_fit0(fit_ref, rep(2L, 6L))

  s1 <- lca_step1_startval(d, Y.names, n_classes = 3L, startval = phi)
  s2 <- lca_step1(d, Y.names, n_classes = 3L, startval = phi)

  expect_equal(s1$fit0$vPi, s2$fit0$vPi)
  expect_equal(s1$fit0$mPhi, s2$fit0$mPhi)
})

test_that("startval phi matrix validates dimensions, probability range, and row sums", {
  d <- generate_data(150L, "high", "distal", seed = 53L)
  Y.names <- paste0("Y", 1:6)
  fit_ref <- lca_step1(d, Y.names, n_classes = 3L)$fit0
  phi <- phi_from_fit0(fit_ref, rep(2L, 6L))

  expect_error(
    lca_step1_startval(d, Y.names, n_classes = 3L, startval = phi[1:10, ]),
    "sum\\(category counts\\)"
  )
  expect_error(
    lca_step1_startval(d, Y.names, n_classes = 3L, startval = phi[, 1:2]),
    "n_classes"
  )
  bad_range <- phi
  bad_range[1L, 1L] <- 1.5
  expect_error(
    lca_step1_startval(d, Y.names, n_classes = 3L, startval = bad_range),
    "\\[0, 1\\]"
  )
  bad_sum <- phi
  bad_sum[1L, ] <- c(0.9, 0.9, 0.9)
  expect_error(
    lca_step1_startval(d, Y.names, n_classes = 3L, startval = bad_sum),
    "sum to ~1"
  )
})

test_that("three_step(startval=phi matrix) fits a covariate model", {
  d <- generate_data(300L, "high", "covariate", seed = 54L)
  Y.names <- paste0("Y", 1:6)
  fit_ref <- lca_step1(d, Y.names, n_classes = 3L)$fit0
  phi <- phi_from_fit0(fit_ref, rep(2L, 6L))

  fit <- three_step(
    d,
    Y.names,
    n_classes = 3L,
    Zp.names = "Zp",
    startval = phi,
    use.simple.cov = TRUE
  )

  expect_s3_class(fit, "tseLCA_covariate")
  expect_equal(dim(fit$three_step), c(2L, 2L))
})

# ---- Step 1: n_init random-classification restarts ---------------------------------------------------------------------------------------

test_that("lca_step1(n_init=) fits from independent random restarts", {
  d <- generate_data(300L, "high", "distal", seed = 55L)
  s1 <- lca_step1(d, paste0("Y", 1:6), n_classes = 3L, n_init = 5L, verbose = FALSE)

  expect_null(s1$fitZ)
  expect_equal(dim(s1$fit0$mPhi), c(6L, 3L))
  expect_equal(sum(s1$fit0$vPi), 1, tolerance = 1e-10)
})

test_that("lca_step1 errors when both startval and n_init are supplied", {
  d <- generate_data(150L, "high", "distal", seed = 56L)
  expect_error(
    lca_step1(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      startval = d$X,
      n_init = 5L
    ),
    "mutually exclusive"
  )
})

test_that("three_step(n_init=) fits measurement and covariate models", {
  d <- generate_data(300L, "high", "covariate", seed = 57L)

  fit_m <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    n_init = 5L,
    use.simple.cov = TRUE
  )
  expect_s3_class(fit_m, "tseLCA_measurement")

  fit_c <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    n_init = 5L,
    use.simple.cov = TRUE
  )
  expect_s3_class(fit_c, "tseLCA_covariate")
})

test_that("three_step errors when step1/startval/n_init are combined", {
  d <- generate_data(150L, "high", "distal", seed = 58L)
  s1 <- lca_step1(d, paste0("Y", 1:6), n_classes = 3L)

  expect_error(
    three_step(d, paste0("Y", 1:6), n_classes = 3L, step1 = s1, n_init = 5L),
    "mutually exclusive"
  )
  expect_error(
    three_step(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      startval = d$X,
      n_init = 5L
    ),
    "mutually exclusive"
  )
})

test_that("fitZ_from_multiLCA(n_init=) fits from independent random restarts", {
  d <- generate_data(300L, "high", "covariate", seed = 59L)
  fZ <- fitZ_from_multiLCA(
    data = d,
    Y.names = paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    maxIter.measurement = 5000L,
    measurement.tol = 1e-8,
    covariate.tol = 1e-6,
    iter.measurement = 10L,
    R2.threshold = 0.70,
    n_init = 5L
  )
  expect_equal(dim(fZ$mGamma), c(2L, 2L))
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

# ---- BCH weight orientation --------------------------------------------------------------------------------------------------------------

test_that("bch_weight_matrix rows sum to 1 and recovers posterior class totals", {
  # Regression test for the BCH weight-matrix orientation bug: pwx[s, t] =
  # P(W = s | X = t) is column-stochastic, and the correct BCH weight
  # matrix is w.is %*% t(pwx)^-1 (Mplus Web Note 21), not w.is %*% pwx^-1.
  # The two orientations only agree when pwx is symmetric, so this uses a
  # deliberately asymmetric classification-error matrix, built the same way
  # compute_pwx_adj() builds it, so the test fails loudly if the `t()` in
  # bch_weight_matrix() is ever removed.
  set.seed(321)
  N <- 3000L
  iT <- 4L
  true_class <- sample(
    seq_len(iT),
    N,
    replace = TRUE,
    prob = c(0.4, 0.3, 0.2, 0.1)
  )
  conc <- c(3, 8, 15, 25)[true_class] # unequal concentration -> asymmetric confusion
  post <- t(vapply(
    seq_len(N),
    function(i) {
      alpha <- rep(1, iT)
      alpha[true_class[i]] <- conc[i]
      p <- rgamma(iT, alpha)
      p / sum(p)
    },
    numeric(iT)
  ))

  w.is <- matrix(0, N, iT)
  w.is[cbind(seq_len(N), max.col(post))] <- 1

  p.wx_joint <- (t(w.is) %*% post) / N
  pwx <- sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")
  expect_equal(colSums(pwx), rep(1, iT), tolerance = 1e-10)
  # Confirm the scenario is genuinely asymmetric -- otherwise both
  # orientations would agree and this test wouldn't discriminate them.
  expect_true(max(abs(pwx - t(pwx))) > 0.02)

  w.it <- bch_weight_matrix(w.is, pwx)

  # Mplus Web Note 21: each case's BCH weights sum to 1.
  expect_equal(rowSums(w.it), rep(1, N))
  # Sample identity: weighted totals recover the posterior class sizes.
  expect_equal(colSums(w.it), colSums(post), tolerance = 1e-6)

  w.it_wrong <- w.is %*% qr.solve(pwx)
  expect_false(isTRUE(all.equal(rowSums(w.it_wrong), rep(1, N))))
  expect_false(
    isTRUE(all.equal(colSums(w.it_wrong), colSums(post), tolerance = 1e-6))
  )
})

test_that("three_step BCH covariate recovers true DGP slopes and intercepts", {
  d <- generate_data(2000L, "high", "covariate", seed = 42L)
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

  #True non-reference class params, sorted by slope ascending: (-1, 1)
  #Align estimated classes to true classes by Zp slope sign
  slopes <- fit$three_step["Zp", ]
  ord <- order(slopes)

  true_intercepts <- c(2.3446, -3.6554) # b0 for (C2, C3) in DGP ordering
  true_slopes <- c(-1, 1)

  ses <- sqrt(diag(fit$three_step_vcov))

  est_int <- fit$three_step["Intercept", ord]
  est_slope <- fit$three_step["Zp", ord]
  se_int <- ses[c(1L, 3L)][ord]
  se_slope <- ses[c(2L, 4L)][ord]

  for (j in 1:2) {
    expect_true(
      abs(est_slope[j] - true_slopes[j]) <= 2 * se_slope[j],
      label = sprintf(
        "BCH slope[%d]: est=%.3f, true=%.3f, 2SE=%.3f",
        j,
        est_slope[j],
        true_slopes[j],
        2 * se_slope[j]
      )
    )
    expect_true(
      abs(est_int[j] - true_intercepts[j]) <= 2 * se_int[j],
      label = sprintf(
        "BCH intercept[%d]: est=%.3f, true=%.3f, 2SE=%.3f",
        j,
        est_int[j],
        true_intercepts[j],
        2 * se_int[j]
      )
    )
  }
})

test_that("three_step BCH gaussian distal recovers true class means", {
  d <- generate_data(2000L, "high", "distal", seed = 88L)
  fit <- suppressWarnings(
    three_step(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      Zo.name = "Zo",
      family = "gaussian",
      use.bch = TRUE,
      use.simple.cov = TRUE,
      verbose = FALSE
    )
  )

  true_mu_sorted <- sort(c(-1, 0, 1))
  est_mu_sorted <- sort(fit$three_step)

  # BCH has higher sampling variance than ML, so this uses a looser 3 SE
  # bound (vs. 2 SE for the analogous ML test) to avoid single-seed
  # flakiness while still checking the estimates are unbiased.
  ses_sorted <- sort(sqrt(diag(fit$three_step_vcov)))
  for (j in seq_along(true_mu_sorted)) {
    expect_true(
      abs(est_mu_sorted[j] - true_mu_sorted[j]) <= 3 * ses_sorted[j],
      label = sprintf(
        "BCH mu[%d]: est=%.3f, true=%.3f, 3SE=%.3f",
        j,
        est_mu_sorted[j],
        true_mu_sorted[j],
        3 * ses_sorted[j]
      )
    )
  }
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
  #True mu = (-1, 0, 1) up to class labeling; range should span negatives and positives
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

# ---- Distal model: multinomial family ---------------------------------------------------------------------------------------------------

make_multinomial_distal_data <- function(n, seed) {
  d <- generate_data(n, "high", "distal", seed = seed)
  # True class-conditional category probabilities: each class has a
  # dominant category, giving a genuinely asymmetric confusion structure.
  pi_true <- matrix(
    c(
      0.70, 0.10, 0.10, 0.10,
      0.10, 0.70, 0.10, 0.10,
      0.10, 0.10, 0.10, 0.70
    ),
    nrow = 3,
    byrow = TRUE
  )
  d$Zcat <- factor(vapply(
    seq_len(nrow(d)),
    function(i) sample(c("a", "b", "c", "d"), 1, prob = pi_true[d$X[i], ]),
    character(1)
  ))
  d
}

test_that("three_step multinomial BCH returns a T x C probability matrix", {
  d <- make_multinomial_distal_data(1500L, seed = 55L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zcat",
    family = "multinomial",
    use.bch = TRUE,
    use.simple.cov = TRUE,
    verbose = FALSE
  )

  expect_s3_class(fit, "tseLCA_distal")
  pi_hat <- coef(fit)
  expect_equal(dim(pi_hat), c(3L, 4L))
  expect_equal(rowSums(pi_hat), setNames(rep(1, 3), rownames(pi_hat)))
  expect_true(all(pi_hat >= 0 & pi_hat <= 1))

  V <- vcov(fit)
  expect_equal(dim(V), c(12L, 12L))
  expect_true(all(is.finite(sqrt(diag(V)))))
})

test_that("print/summary don't error on a multinomial tseLCA_distal object", {
  d <- make_multinomial_distal_data(300L, seed = 55L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zcat",
    family = "multinomial",
    use.simple.cov = TRUE
  )
  expect_output(print(fit), "multinomial")
  expect_output(print(summary(fit)), "multinomial")
})

test_that("three_step multinomial ML (simple and full) agree and full SEs are >= simple SEs", {
  d <- make_multinomial_distal_data(1500L, seed = 55L)
  Y.names <- paste0("Y", 1:6)

  fit_simple <- three_step(
    d,
    Y.names,
    n_classes = 3L,
    Zo.name = "Zcat",
    family = "multinomial",
    use.bch = FALSE,
    use.simple.cov = TRUE
  )
  fit_full <- three_step(
    d,
    Y.names,
    n_classes = 3L,
    Zo.name = "Zcat",
    family = "multinomial",
    use.bch = FALSE,
    use.simple.cov = FALSE
  )

  expect_equal(coef(fit_simple), coef(fit_full), tolerance = 1e-6)
  se_simple <- sqrt(diag(vcov(fit_simple)))
  se_full <- sqrt(diag(vcov(fit_full)))
  # Full Step-1 propagation should add uncertainty on top of the robust
  # sandwich, not remove it.
  expect_true(all(se_full >= se_simple - 1e-8))
})

test_that("multinomial_ml_jacobian matches a numerical check of the estimating equation", {
  # multinomial_ml_jacobian() is a closed-form derivation (no numerical
  # differentiation in the package); this is a cheap regression check
  # against a numerical Jacobian computed inline, entirely independent of
  # multinomial_ml_jacobian() itself.
  d <- make_multinomial_distal_data(600L, seed = 60L)
  # Pre-encode as integer categories, matching what three_step() itself
  # does internally before calling clean_data() -- required here since we
  # call clean_data() directly below.
  d$Zcat <- as.integer(factor(d$Zcat))
  Y.names <- paste0("Y", 1:6)

  s1 <- lca_step1(d, Y.names, n_classes = 3L)
  fit0 <- s1$fit0
  cd <- clean_data(
    data = d, Y.names = Y.names, Zo.name = "Zcat",
    incomplete = FALSE, include.intercept = TRUE, verbose = FALSE
  )
  s2 <- lca_step2(
    cd$Y.obs, fit0, 3L, TRUE, 1e-2, FALSE,
    ivItemcat = cd$ivItemcat, mDesign = cd$mDesign
  )
  w.is <- s2$w.is[cd$keep_step3_Zo_in_Y, , drop = FALSE]
  pwx <- s2$p.wx_mat
  Y_cat <- cd$Zo_mat[, 1L]
  C <- length(unique(Y_cat))
  pi_adj <- matrix(fit0$vPi, nrow = length(Y_cat), ncol = 3L, byrow = TRUE)

  s3 <- lca_step3.distal.multinomial(
    Y_cat = Y_cat, C = C, iT = 3L, covariate.tol = 1e-8, use.bch = FALSE,
    w.is_cc = w.is, pwx = pwx, em.maxIter = 500L, vPi = fit0$vPi, pi_mat = pi_adj
  )
  theta_hat <- s3$res$par
  Psi <- function(theta) colSums(s3$three_step.score(theta))

  n <- length(theta_hat)
  Jac_numeric <- matrix(0, n, n)
  eps <- 1e-5
  for (j in seq_len(n)) {
    h <- eps * max(abs(theta_hat[j]), 1)
    xp <- theta_hat; xp[j] <- xp[j] + h
    xm <- theta_hat; xm[j] <- xm[j] - h
    Jac_numeric[, j] <- (Psi(xp) - Psi(xm)) / (2 * h)
  }

  pi_hat <- matrix(theta_hat, nrow = 3L, ncol = C)
  pzx_mat <- t(pi_hat[, Y_cat, drop = FALSE])
  ae <- w.is %*% pwx
  q_i <- rowSums(pi_adj * pzx_mat * ae)
  r_it <- pi_adj * pzx_mat * ae / q_i
  Jac_analytic <- multinomial_ml_jacobian(pi_hat, r_it, Y_cat)

  expect_equal(Jac_analytic, Jac_numeric, tolerance = 1e-4)
})

test_that("multinomial with C=2 matches the existing binomial family exactly", {
  d <- make_multinomial_distal_data(1200L, seed = 56L)
  Y.names <- paste0("Y", 1:6)
  d$Zbin <- rbinom(nrow(d), 1, c(0.2, 0.5, 0.8)[d$X])
  d$Zbin_cat <- factor(d$Zbin)

  fit_bin <- three_step(
    d,
    Y.names,
    n_classes = 3L,
    Zo.name = "Zbin",
    family = "binomial",
    use.bch = TRUE,
    use.simple.cov = TRUE
  )
  fit_multi <- three_step(
    d,
    Y.names,
    n_classes = 3L,
    Zo.name = "Zbin_cat",
    family = "multinomial",
    use.bch = TRUE,
    use.simple.cov = TRUE
  )

  mu_bin <- 1 / (1 + exp(-coef(fit_bin)))
  pi_multi <- coef(fit_multi)[, "1"]
  # Both estimators solve the same weighted-proportion closed form, but
  # binomial's BCH path gets there with Newton-Raphson (converged to within
  # covariate.tol on the parameter step, default 1e-6) while multinomial's
  # BCH is a direct closed form, so residual agreement is to ~1e-6-1e-7,
  # not machine precision.
  expect_equal(unname(mu_bin), unname(pi_multi), tolerance = 1e-5)
})

test_that("multinomial recovers the true category probability structure", {
  d <- make_multinomial_distal_data(3000L, seed = 57L)
  fit <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zcat",
    family = "multinomial",
    use.simple.cov = TRUE
  )
  pi_hat <- coef(fit)
  # Each class should have exactly one dominant (> 0.5) category.
  dominant <- apply(pi_hat, 1L, max)
  expect_true(all(dominant > 0.5))
  # The three dominant categories should be distinct (classes recovered
  # three different modes, matching the DGP's three distinct classes).
  expect_length(unique(apply(pi_hat, 1L, which.max)), 3L)
})

test_that("family = \"multinomial\" validates the category count", {
  d <- make_multinomial_distal_data(200L, seed = 58L)
  d$Zconst <- factor(rep("a", nrow(d)))

  expect_error(
    three_step(
      d,
      paste0("Y", 1:6),
      n_classes = 3L,
      Zo.name = "Zconst",
      family = "multinomial"
    ),
    "at least 2 distinct categories"
  )
})

test_that("combined Zp.names + family = \"multinomial\" works under full propagation", {
  d2 <- generate_data(200L, "high", "covariate", seed = 59L)
  d2$Zcat <- factor(sample(c("a", "b", "c"), nrow(d2), replace = TRUE))

  fit_full <- three_step(
    d2,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    Zo.name = "Zcat",
    family = "multinomial",
    use.simple.cov = FALSE
  )
  expect_s3_class(fit_full, "tseLCA_both")
  V_full <- vcov(fit_full, which = "distal")
  expect_true(all(is.finite(diag(V_full))))
  expect_true(all(diag(V_full) > 0))

  fit_simple <- three_step(
    d2,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    Zo.name = "Zcat",
    family = "multinomial",
    use.simple.cov = TRUE
  )
  expect_s3_class(fit_simple, "tseLCA_both")
  # Full Step-1/Step-2 propagation should add uncertainty on top of the
  # robust sandwich, not remove it (same check as the distal-only case).
  se_full <- sqrt(diag(V_full))
  se_simple <- sqrt(diag(vcov(fit_simple, which = "distal")))
  expect_true(all(se_full >= se_simple - 1e-8))
})

test_that("Step-2 covariate-uncertainty propagation is actually wired through for multinomial", {
  # The se_full >= se_simple check above would still pass even if
  # three_step() silently stopped passing Sigma.3/s3.par/p.xz.cov/Z_mat_cov
  # to lca_vcov_distal_multinomial() (it would just fall back to a
  # step1-only propagation, still >= the robust sandwich). This test
  # isolates the step-2 term's contribution directly through the internal
  # function, so a regression that drops the wiring shows up as "no
  # difference" rather than merely "still bigger than simple.cov".
  d2 <- generate_data(400L, "high", "covariate", seed = 59L)
  # Pre-encode as integer categories, matching what three_step() itself
  # does internally before calling clean_data() -- required here since we
  # call clean_data() directly below to reconstruct the internals.
  d2$Zcat <- as.integer(factor(sample(c("a", "b", "c"), nrow(d2), replace = TRUE)))
  Y.names <- paste0("Y", 1:6)

  fit <- three_step(
    d2, Y.names, n_classes = 3L, Zp.names = "Zp", Zo.name = "Zcat",
    family = "multinomial", use.bch = FALSE, use.simple.cov = FALSE
  )

  s1 <- lca_step1(d2, Y.names, n_classes = 3L)
  fit0 <- s1$fit0
  cd <- clean_data(
    data = d2, Y.names = Y.names, Zp.names = "Zp", Zo.name = "Zcat",
    incomplete = FALSE, include.intercept = TRUE, verbose = FALSE
  )
  s2 <- lca_step2(
    cd$Y.obs, fit0, 3L, TRUE, 1e-2, FALSE,
    ivItemcat = cd$ivItemcat, mDesign = cd$mDesign
  )
  J.2_dis <- s2$compute_J_unc(
    s2$p.xy[cd$keep_step3_Zo_in_Y, , drop = FALSE],
    cd$Y.obs[cd$keep_step3_Zo_in_Y, , drop = FALSE],
    matrix(1L, length(cd$keep_step3_Zo_in_Y), ncol(cd$Y.obs)),
    s2$theta1, cd$ivItemcat, 3L
  )
  s2_for_dis <- list(
    J.2 = J.2_dis, p.wx_mat = s2$p.wx_mat,
    w.is = s2$w.is[cd$keep_step3_Zo_in_Y, , drop = FALSE]
  )

  s3.par <- as.vector(fit$covariate$three_step)
  Z_full_raw <- cbind(1, as.matrix(d2[, "Zp", drop = FALSE]))
  Z_mat_dis <- Z_full_raw[cd$keep_step3_Zo, , drop = FALSE]
  p.xz_dis <- function(params) {
    eta_full <- cbind(0, Z_mat_dis %*% params)
    ex <- exp(eta_full - apply(eta_full, 1L, max))
    ex / rowSums(ex)
  }
  pi_adj <- p.xz_dis(matrix(s3.par, ncol = 2))

  res_adj <- compute_pwx_adj(
    cd$Y.obs[cd$keep_step3_Zo_in_Y, , drop = FALSE], fit0, cd$ivItemcat,
    NULL, TRUE, pi_adj = pi_adj
  )
  Y_cat_dis <- cd$Zo_mat[, 1L]
  C <- length(levels(factor(d2$Zcat)))

  s3.distal <- lca_step3.distal.multinomial(
    Y_cat = Y_cat_dis, C = C, iT = 3L, covariate.tol = 1e-8, use.bch = FALSE,
    w.is_cc = res_adj$w.is, pwx = res_adj$p.wx_mat, em.maxIter = 500L,
    vPi = fit0$vPi, pi_mat = pi_adj
  )
  Sigma.1 <- lca_indiv_varmat(
    cd$Y.obs, cd$mDesign, fit0, cd$ivItemcat, boundary.tol = 1e-2
  )$Varmat

  V_step1_only <- lca_vcov_distal_multinomial(
    theta_hat = s3.distal$res$par, three_step.score = s3.distal$three_step.score,
    pi_adj = pi_adj, w.is = res_adj$w.is, p.wx_mat = res_adj$p.wx_mat,
    Y_cat = Y_cat_dis, C = C, H.3.inv = s3.distal$H.3.inv, Sigma.1 = Sigma.1,
    s2 = s2_for_dis, iT = 3L, use.simple.cov = FALSE, use.bch = FALSE
    # Sigma.3/s3.par/p.xz.cov/Z_mat_cov omitted -> step1 term only
  )

  V_full <- vcov(fit, which = "distal")
  expect_false(isTRUE(all.equal(diag(V_step1_only), diag(V_full))))
  # Adding the step-2 term should increase (not decrease) the variance.
  expect_true(all(diag(V_full) >= diag(V_step1_only) - 1e-10))
})

# ---- Omnibus test -------------------------------------------------------------------------------------------------------------------------

test_that("omnibus_test degrees of freedom match theory for multinomial and scalar families", {
  d <- make_multinomial_distal_data(1500L, seed = 55L)
  fit_multi <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zcat",
    family = "multinomial",
    use.simple.cov = TRUE
  )
  ob_multi <- omnibus_test(fit_multi)
  expect_s3_class(ob_multi, "tseLCA_omnibus")
  # (T-1)*(C-1) = 2*3 = 6, the textbook df for a T x C homogeneity test.
  expect_equal(ob_multi$df, 6L)
  expect_true(ob_multi$p.value < 0.001)

  fit_gauss <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zo",
    use.simple.cov = TRUE
  )
  ob_gauss <- omnibus_test(fit_gauss)
  expect_equal(ob_gauss$df, 2L) # T - 1
})

test_that("omnibus_test does not reject when classes share the same distribution", {
  d <- make_multinomial_distal_data(1500L, seed = 55L)
  d$Zcat_null <- factor(sample(
    c("a", "b", "c", "d"),
    nrow(d),
    replace = TRUE,
    prob = c(0.4, 0.3, 0.2, 0.1)
  ))
  fit_null <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zo.name = "Zcat_null",
    family = "multinomial",
    use.simple.cov = TRUE
  )
  ob_null <- omnibus_test(fit_null)
  expect_true(ob_null$p.value > 0.10)
})

test_that("omnibus_test works on a tseLCA_both object's distal component", {
  d <- generate_data(1500L, "high", "covariate", seed = 66L)
  pi_true <- matrix(
    c(
      0.70, 0.10, 0.10, 0.10,
      0.10, 0.70, 0.10, 0.10,
      0.10, 0.10, 0.10, 0.70
    ),
    nrow = 3,
    byrow = TRUE
  )
  d$Zcat <- factor(vapply(
    seq_len(nrow(d)),
    function(i) sample(c("a", "b", "c", "d"), 1, prob = pi_true[d$X[i], ]),
    character(1)
  ))
  fit_both <- three_step(
    d,
    paste0("Y", 1:6),
    n_classes = 3L,
    Zp.names = "Zp",
    Zo.name = "Zcat",
    family = "multinomial",
    use.simple.cov = TRUE
  )
  ob <- omnibus_test(fit_both)
  expect_s3_class(ob, "tseLCA_omnibus")
  expect_equal(ob$df, 6L)
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
# Class labels from the estimator may differ from the DGP labeling.
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
