# tests/testthat/test-distal-ml.R
#
# Three-step ML for distal outcomes (R/distal-ml.R):
#   * the log-likelihood is sum_i sum_s w_is log M_is over the expanded data,
#     so proportional-assignment weights stay outside the log (tseLCA 1.x
#     put them inside, biasing proportional-assignment estimates);
#   * the gaussian family estimates the within-class variance sigma2 jointly
#     with the means (fixed at 1 in 1.x).
# The analytic score, Hessian/Jacobian, and Step-1/Step-2 cross-derivatives
# are checked against numerical derivatives for every family, using
# proportional assignment weights so that the per-record form is exercised
# (with modal weights it reduces to the 1.x formulas).

# -- synthetic Step-3 inputs, independent of any Step-1 fit ---------------------

ml_inputs <- function(N = 300L, seed = 1L) {
  set.seed(seed)
  iT <- 3L
  X <- sample.int(iT, N, replace = TRUE, prob = c(.3, .3, .4))
  pwx <- matrix(c(.8, .1, .1, .15, .7, .15, .1, .2, .7), iT) # columns sum to 1
  # proportional assignment: noisy posterior rows centred on the true class
  w.is <- t(vapply(X, function(x) {
    p <- pwx[, x] + stats::rexp(iT, 8)
    p / sum(p)
  }, numeric(iT)))
  Zc <- cbind(1, rnorm(N))
  list(
    N = N, iT = iT, X = X, pwx = pwx, w.is = w.is, Zc = Zc,
    gamma = c(0.2, -0.3, 0.4, 0.5),
    z = list(
      gaussian = rnorm(N, c(-1, 1, 0)[X], 1.5),
      poisson = rpois(N, c(1, 3, 6)[X]),
      binomial = rbinom(N, 1L, c(.2, .5, .8)[X]),
      multinomial = vapply(X, function(x) {
        sample.int(3L, 1L, prob = rbind(c(.6, .3, .1), c(.2, .6, .2), c(.1, .3, .6))[x, ])
      }, integer(1))
    ),
    theta = list(
      gaussian = c(-0.8, 0.9, 0.1, 2.1),
      poisson = log(c(1.2, 2.7, 5.5)),
      binomial = qlogis(c(.25, .45, .75)),
      multinomial = as.vector(rbind(c(.5, .3, .2), c(.25, .5, .25), c(.15, .25, .6)))
    )
  )
}

p_xz <- function(gamma, Zc) {
  eta <- cbind(0, Zc %*% matrix(gamma, ncol = 2L))
  e <- exp(eta - apply(eta, 1L, max))
  e / rowSums(e)
}

# column-softmax parameterization of P(W = s | X = t): off-diagonal
# log-ratios log(pwx[s, t] / pwx[t, t]), ordered t0 then s0
theta2_of <- function(pwx) {
  iT <- ncol(pwx)
  unlist(lapply(seq_len(iT), function(t) log(pwx[-t, t] / pwx[t, t])))
}
pwx_of <- function(theta2, iT) {
  M <- matrix(0, iT, iT)
  k <- 0L
  for (t in seq_len(iT)) {
    M[-t, t] <- theta2[k + seq_len(iT - 1L)]
    k <- k + iT - 1L
  }
  E <- exp(M)
  sweep(E, 2L, colSums(E), "/")
}

num_jac <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(j) {
    h <- eps * max(abs(x[j]), 1)
    xp <- x
    xp[j] <- xp[j] + h
    xm <- x
    xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, numeric(length(f(x))))
}

# log f(z_i | X = t), written independently of the package
log_f <- function(family, theta, z, iT = 3L) {
  switch(
    family,
    gaussian = {
      s2 <- theta[iT + 1L]
      -0.5 * outer(z, theta[1:iT], "-")^2 / s2 - 0.5 * log(2 * pi * s2)
    },
    poisson = outer(z, theta, "*") - matrix(exp(theta), length(z), iT, byrow = TRUE) -
      lgamma(z + 1),
    binomial = outer(z, plogis(theta), function(y, p) y * log(p) + (1 - y) * log(1 - p)),
    multinomial = log(t(matrix(theta, nrow = iT)[, z, drop = FALSE]))
  )
}

# expanded-data log-likelihood, sum_i sum_s w_is log sum_t a_it f_it pwx[s, t]
loglik_expanded <- function(family, theta, z, a, pwx, w.is) {
  f <- exp(log_f(family, theta, z))
  M <- (a * f) %*% t(pwx) # N x S
  sum(w.is * log(M))
}

# estimating equation Psi(theta) (score for gaussian/poisson/binomial; the
# softmax-gradient form for multinomial)
Psi <- function(family, theta, g, a = NULL, pwx = g$pwx) {
  if (is.null(a)) a <- p_xz(g$gamma, g$Zc)
  z <- g$z[[family]]
  rec <- distal_records(log_f(family, theta, z), a, pwx)
  G <- distal_unit_derivs(theta, z, g$iT, family, C = 3L)$G
  colSums(distal_score(distal_lambda(rec$R, g$w.is), G))
}

# -- tests ----------------------------------------------------------------------

test_that("distal_records/distal_loglik give the expanded-data log-likelihood", {
  g <- ml_inputs()
  a <- p_xz(g$gamma, g$Zc)
  for (fam in c("gaussian", "poisson", "binomial", "multinomial")) {
    z <- g$z[[fam]]
    rec <- distal_records(log_f(fam, g$theta[[fam]], z), a, g$pwx)
    expect_equal(distal_loglik(rec$logM, g$w.is),
                 loglik_expanded(fam, g$theta[[fam]], z, a, g$pwx, g$w.is),
                 info = fam)
    # E-step weights are proper: each row sums to sum_s w_is = 1
    expect_equal(rowSums(distal_lambda(rec$R, g$w.is)), rep(1, g$N), info = fam)
  }
})

test_that("score equals the gradient of the log-likelihood", {
  g <- ml_inputs()
  a <- p_xz(g$gamma, g$Zc)
  for (fam in c("gaussian", "poisson", "binomial")) {
    z <- g$z[[fam]]
    grad <- num_jac(function(th) loglik_expanded(fam, th, z, a, g$pwx, g$w.is), g$theta[[fam]])
    expect_equal(Psi(fam, g$theta[[fam]], g), as.vector(grad), tolerance = 1e-6, info = fam)
  }
})

test_that("analytic Hessian / Jacobian match numerical derivatives", {
  g <- ml_inputs()
  a <- p_xz(g$gamma, g$Zc)
  for (fam in c("gaussian", "poisson", "binomial")) {
    th <- g$theta[[fam]]
    z <- g$z[[fam]]
    rec <- distal_records(log_f(fam, th, z), a, g$pwx)
    H <- distal_neg_hessian(rec, g$w.is, distal_unit_derivs(th, z, g$iT, fam))
    H_num <- -num_jac(function(t) Psi(fam, t, g), th)
    expect_equal(H, H_num, tolerance = 1e-6, info = fam)
  }
  th <- g$theta$multinomial
  z <- g$z$multinomial
  rec <- distal_records(log_f("multinomial", th, z), a, g$pwx)
  Jac <- distal_multinomial_jacobian(matrix(th, nrow = 3L), rec, g$w.is, z)
  expect_equal(Jac, num_jac(function(t) Psi("multinomial", t, g), th), tolerance = 1e-6)
})

test_that("Step-1 and Step-2 cross-derivatives match numerical derivatives", {
  g <- ml_inputs()
  a <- p_xz(g$gamma, g$Zc)
  for (fam in c("gaussian", "poisson", "binomial", "multinomial")) {
    th <- g$theta[[fam]]
    z <- g$z[[fam]]
    rec <- distal_records(log_f(fam, th, z), a, g$pwx)
    G <- distal_unit_derivs(th, z, g$iT, fam, C = 3L)$G
    cross <- distal_cross_derivs(rec, g$w.is, g$pwx, G, Z_mat = g$Zc)
    C1_num <- num_jac(function(t2) Psi(fam, th, g, pwx = pwx_of(t2, g$iT)), theta2_of(g$pwx))
    C_num <- num_jac(function(gm) Psi(fam, th, g, a = p_xz(gm, g$Zc)), g$gamma)
    expect_equal(cross$C1, C1_num, tolerance = 1e-5, info = fam)
    expect_equal(cross$C_mat, C_num, tolerance = 1e-5, info = fam)
  }
})

test_that("lca_vcov_distal assembles the propagated sandwich (gaussian)", {
  g <- ml_inputs()
  iT <- g$iT
  a <- p_xz(g$gamma, g$Zc)
  z <- g$z$gaussian
  th <- g$theta$gaussian
  p.zx <- function(params) log_f("gaussian", params, z)
  score_fn <- function(params) {
    rec <- distal_records(p.zx(params), a, g$pwx)
    distal_score(distal_lambda(rec$R, g$w.is), distal_unit_derivs(params, z, iT, "gaussian")$G)
  }
  H.3.inv <- solve(distal_neg_hessian(
    distal_records(p.zx(th), a, g$pwx), g$w.is, distal_unit_derivs(th, z, iT, "gaussian")
  ))
  set.seed(2L)
  A <- matrix(rnorm(36), 6)
  Sigma.1 <- crossprod(A) / 500
  B <- matrix(rnorm(16), 4)
  Sigma.3 <- crossprod(B) / 500

  C1_num <- num_jac(function(t2) Psi("gaussian", th, g, pwx = pwx_of(t2, iT)), theta2_of(g$pwx))
  C_num <- num_jac(function(gm) Psi("gaussian", th, g, a = p_xz(gm, g$Zc)), g$gamma)
  V_num <- H.3.inv %*%
    (crossprod(score_fn(th)) + C1_num %*% Sigma.1 %*% t(C1_num) +
       C_num %*% Sigma.3 %*% t(C_num)) %*%
    H.3.inv

  V <- lca_vcov_distal(
    mu_hat = th, three_step.score = score_fn, pi_adj = a, w.is = g$w.is,
    p.wx_mat = g$pwx, p.zx = p.zx, family = "gaussian", H.3.inv = H.3.inv,
    Sigma.1 = Sigma.1, s2 = list(J.2 = diag(6)), Sigma.3 = Sigma.3,
    s3.par = g$gamma, p.xz.cov = function(params) p_xz(as.vector(params), g$Zc),
    Z_mat_cov = g$Zc, iT = iT, use.simple.cov = FALSE, use.bch = FALSE,
    unit_scores = function(params) distal_unit_derivs(params, z, iT, "gaussian")$G
  )
  expect_equal(unname(V), unname(V_num), tolerance = 1e-5)
  expect_equal(rownames(V), c("mu_C1", "mu_C2", "mu_C3", "sigma2"))
})

test_that("modal assignment reduces the per-record form to one record per person", {
  g <- ml_inputs()
  a <- p_xz(g$gamma, g$Zc)
  w_modal <- diag(3)[max.col(g$w.is), ]
  z <- g$z$gaussian
  th <- g$theta$gaussian
  rec <- distal_records(log_f("gaussian", th, z), a, g$pwx)
  ae <- w_modal %*% g$pwx
  joint <- a * exp(log_f("gaussian", th, z)) * ae
  expect_equal(distal_lambda(rec$R, w_modal), joint / rowSums(joint))
  expect_equal(distal_loglik(rec$logM, w_modal), sum(log(rowSums(joint))))
})

test_that("gaussian distal estimates are scale equivariant (ML and BCH)", {
  d <- generate_data(600L, "mid", "distal", seed = 7L)
  d$Zo10 <- 10 * d$Zo
  set.seed(1L)
  m <- three_step(d, paste0("Y", 1:6), 3L)
  for (bch in c(FALSE, TRUE)) {
    for (modal in c(TRUE, FALSE)) {
      f1 <- three_step(d, paste0("Y", 1:6), 3L, Zo.name = "Zo", step1 = m,
                       use.bch = bch, use.modal.assignment = modal)
      f10 <- three_step(d, paste0("Y", 1:6), 3L, Zo.name = "Zo10", step1 = m,
                        use.bch = bch, use.modal.assignment = modal)
      info <- paste("bch =", bch, "modal =", modal)
      expect_equal(coef(f10), 10 * coef(f1), tolerance = 1e-4, info = info)
      expect_equal(sqrt(diag(vcov(f10))), 10 * sqrt(diag(vcov(f1))),
                   tolerance = 1e-4, info = info)
      expect_equal(f10$sigma2[["estimate"]], 100 * f1$sigma2[["estimate"]],
                   tolerance = 1e-4, info = info)
    }
  }
})

test_that("ML recovers means and variance of a wide outcome (modal and proportional)", {
  d <- generate_data(3000L, "high", "distal", seed = 11L)
  mu <- c(-1, 1, 0)[d$X]
  d$Zo3 <- mu + 3 * (d$Zo - mu) # same class means, residual SD 3
  set.seed(1L)
  m <- three_step(d, paste0("Y", 1:6), 3L)
  perm <- apply(table(classes(m), d$X), 2L, which.max)
  for (modal in c(TRUE, FALSE)) {
    f <- three_step(d, paste0("Y", 1:6), 3L, Zo.name = "Zo3", step1 = m,
                    use.modal.assignment = modal)
    expect_equal(unname(coef(f)[perm]), c(-1, 1, 0), tolerance = 0.25, info = modal)
    expect_equal(f$sigma2[["estimate"]], 9, tolerance = 0.1, info = modal)
    expect_true(is.finite(f$sigma2[["se"]]))
  }
})
