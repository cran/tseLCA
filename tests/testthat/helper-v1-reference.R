# tests/testthat/helper-v1-reference.R
#
# Shared definitions for the v1.1.1 numerical regression fixtures
# (tests/testthat/fixtures/v1_reference.rds), used by
# test-regression-v1.R and by fixtures/make_v1_reference.R.
#
# The fixtures pin the numbers produced by tseLCA 1.1.1, as amended by the
# documented 2.0 bug fixes (see the header of fixtures/make_v1_reference.R),
# so the 2.0 refactor can be checked to leave estimates unchanged. `v1_fit()` builds each fit and
# `v1_extract()` pulls out the numbers to compare. When the object structure
# changes during the refactor, update these two functions, not the fixtures.
# Fixtures are regenerated only for deliberate numerical changes (bug fixes),
# each documented in NEWS.md.

v1_items <- paste0("Y", 1:6)

# Synthetic polytomous data: 3 classes, 6 items with 3 categories each.
v1_poly_data <- function(n = 600L, seed = 11L) {
  set.seed(seed)
  X <- sample.int(3L, n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  probs <- list(
    c(0.70, 0.20, 0.10),
    c(0.15, 0.70, 0.15),
    c(0.10, 0.20, 0.70)
  )
  Y <- sapply(1:6, function(j) {
    vapply(X, function(x) sample(0:2, 1L, prob = probs[[x]]), integer(1))
  })
  colnames(Y) <- v1_items
  d <- data.frame(Y, X = X)
  d$Zp <- sample(1:5, n, replace = TRUE)
  d
}

# Synthetic polytomous data where some item-response probabilities sit at
# multilevLCA's 1e-5 floor (category 2 of Y1/Y2 is never chosen in class 1).
v1_sparse_poly_data <- function(n = 800L, seed = 3L) {
  set.seed(seed)
  X <- sample.int(3L, n, replace = TRUE)
  sparse <- list(c(.8, .2, 0), c(.2, .6, .2), c(.1, .2, .7))
  dense <- list(c(.6, .3, .1), c(.1, .6, .3), c(.3, .1, .6))
  Y <- sapply(1:6, function(j) {
    pr <- if (j <= 2L) sparse else dense
    vapply(X, function(x) sample(0:2, 1L, prob = pr[[x]]), integer(1))
  })
  colnames(Y) <- v1_items
  d <- data.frame(Y, X = X)
  d$Zp <- sample(1:5, n, replace = TRUE)
  d
}

# Data sets, all synthetic.
v1_data <- function() {
  cov_high <- generate_data(500L, "high", "covariate", seed = 101L)
  cov_mid <- generate_data(500L, "mid", "covariate", seed = 102L)
  dis <- generate_data(500L, "high", "distal", seed = 103L)
  set.seed(104L)
  dis$Zpois <- rpois(nrow(dis), lambda = c(1, 3, 6)[dis$X])
  dis$Zbin <- rbinom(nrow(dis), 1L, prob = c(0.2, 0.5, 0.8)[dis$X])
  cat_p <- rbind(c(.7, .15, .15), c(.15, .7, .15), c(.15, .15, .7))
  dis$Zcat <- factor(vapply(
    dis$X, function(x) sample(c("a", "b", "c"), 1L, prob = cat_p[x, ]),
    character(1)
  ))
  both <- generate_data(500L, "high", "covariate", seed = 105L)
  both$Zo <- draw_Zo(both$X, bk2018_params$distal_params)
  miss <- generate_data(500L, "high", "covariate", seed = 106L)
  set.seed(107L)
  M <- matrix(runif(nrow(miss) * 6) < 0.1, ncol = 6)
  for (j in 1:6) miss[[v1_items[j]]][M[, j]] <- NA
  list(
    cov_high = cov_high, cov_mid = cov_mid, dis = dis, both = both,
    miss = miss, poly = v1_poly_data(), sparse = v1_sparse_poly_data()
  )
}

# Configurations: data set + three_step() arguments.
v1_configs <- list(
  meas_high          = list(data = "cov_high", args = list()),
  meas_poly          = list(data = "poly", args = list()),
  cov_ml_modal       = list(data = "cov_high", args = list(Zp.names = "Zp")),
  cov_ml_prop        = list(data = "cov_mid", args = list(Zp.names = "Zp", use.modal.assignment = FALSE)),
  cov_bch_modal      = list(data = "cov_high", args = list(Zp.names = "Zp", use.bch = TRUE)),
  cov_bch_prop       = list(data = "cov_mid", args = list(Zp.names = "Zp", use.bch = TRUE, use.modal.assignment = FALSE)),
  cov_ml_simple      = list(data = "cov_high", args = list(Zp.names = "Zp", use.simple.cov = TRUE)),
  cov_rebase_C2      = list(data = "cov_high", args = list(Zp.names = "Zp", rebase = "C2")),
  cov_step1          = list(data = "cov_high", args = list(Zp.names = "Zp"), step1 = TRUE),
  cov_twostep_vcov   = list(data = "cov_high", args = list(Zp.names = "Zp", get.twostep.vcov = TRUE)),
  cov_fiml           = list(data = "miss", args = list(Zp.names = "Zp", incomplete = TRUE, use.two.step = FALSE)),
  cov_poly           = list(data = "poly", args = list(Zp.names = "Zp", use.simple.cov = TRUE)),
  cov_poly_corrected = list(data = "poly", args = list(Zp.names = "Zp")),
  cov_sparse         = list(data = "sparse", args = list(Zp.names = "Zp")),
  dis_gauss_ml       = list(data = "dis", args = list(Zo.name = "Zo", use.modal.assignment = FALSE)),
  dis_gauss_bch      = list(data = "dis", args = list(Zo.name = "Zo", use.bch = TRUE)),
  dis_poisson        = list(data = "dis", args = list(Zo.name = "Zpois", family = "poisson")),
  dis_binomial       = list(data = "dis", args = list(Zo.name = "Zbin", family = "binomial")),
  dis_multinomial    = list(data = "dis", args = list(Zo.name = "Zcat", family = "multinomial")),
  both_ml_prop       = list(data = "both", args = list(Zp.names = "Zp", Zo.name = "Zo", use.modal.assignment = FALSE))
)

# Fit one configuration. Seeds are reset before every fit because Step 1
# initialization is stochastic.
v1_fit <- function(cfg, data_list) {
  d <- data_list[[cfg$data]]
  args <- c(list(data = d, Y.names = v1_items, n_classes = 3L), cfg$args)
  if (isTRUE(cfg$step1)) {
    set.seed(1L)
    m <- three_step(d, Y.names = v1_items, n_classes = 3L)
    args$step1 <- m$measurement_model
  }
  set.seed(1L)
  suppressWarnings(suppressMessages(do.call(three_step, args)))
}

# Extract the numbers to compare.
v1_extract <- function(fit) {
  strip <- function(x) if (is.null(x)) NULL else unclass(x)
  if (inherits(fit, "tseLCA_measurement")) {
    f0 <- fit$measurement_model$fit0
    return(list(
      llik = fit$llik, AIC = fit$AIC, BIC = fit$BIC, R2entr = fit$R2entr,
      vPi = unname(f0$vPi), mPhi = unname(f0$mPhi),
      posteriors = unname(fit$posteriors)
    ))
  }
  structural <- function(x) {
    list(
      est = strip(x$three_step), vcov = strip(x$three_step_vcov),
      two_step = strip(x$two_step), two_step_vcov = strip(x$two_step_vcov),
      llik = x$llik, AIC = x$AIC, BIC = x$BIC, sigma2 = x$sigma2
    )
  }
  out <- if (inherits(fit, "tseLCA_both")) {
    list(covariate = structural(fit$covariate), distal = structural(fit$distal))
  } else {
    structural(fit)
  }
  out$posteriors <- unname(fit$posteriors)
  out
}
