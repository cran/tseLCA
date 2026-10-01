# tseLCA/R/step3-distal.R
#
# Step 3 for distal outcomes: class-specific outcome parameters (BCH or ML).
# The ML likelihood machinery shared by all families is in R/distal-ml.R.

#' Step 3 (distal, multinomial): estimate class-conditional category
#' probabilities for a nominal categorical distal outcome
#'
#' Estimates the T x C matrix `pi_hat[t, c] = P(Zo = c | X = t)` for a
#' saturated (nominal) multinomial distal outcome with C categories, using the
#' closed-form weighted-proportion estimator
#' `pi_hat[t, c] = sum_i w_it * 1(y_i = c) / sum_i w_it`, using either the
#' BCH weight matrix (fixed given Step 1/2) or, for ML, EM-updated posterior
#' responsibilities (both give the same closed-form M-step; only the E-step
#' weights differ). Returns the same contract as `lca_step3.distal()`
#' (`res$par`, `H.3.inv`, `three_step.score`), with `res$par` the
#' column-major-flattened `pi_hat` (`matrix(par, nrow = iT, ncol = C)`
#' recovers it) so downstream sandwich-variance code is unchanged.
#'
#' For BCH, the score is the moment/estimating-equation form the weighted
#' proportion solves, `s_i,(t,c) = w_it * (1(y_i=c) - pi_hat[t,c])`; its
#' bread is exactly `diag(1 / colSums(w_it))` (repeated across categories)
#' since `w_it` does not depend on `pi_hat`. For ML, the same estimating
#' equation is used with `w_it` replaced by the (theta-dependent) E-step
#' weights `lambda_it = sum_s w_is R_ist` of the expanded-data likelihood
#' (R/distal-ml.R); its bread is `solve(-Jacobian(Psi))`, with the Jacobian
#' given in closed form by `distal_multinomial_jacobian()`.
#' @noRd
lca_step3.distal.multinomial <- function(
  Y_cat,
  C,
  iT,
  covariate.tol,
  use.bch = FALSE,
  w.is_cc = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  vPi = NULL,
  pi_mat = NULL,
  verbose = FALSE
) {
  N <- length(Y_cat)
  onehot_y <- matrix(0, N, C)
  onehot_y[cbind(seq_len(N), Y_cat)] <- 1

  pi_s <- if (!is.null(pi_mat)) {
    pi_mat
  } else {
    matrix(vPi, ncol = iT, nrow = N, byrow = TRUE)
  }

  # T x C closed-form weighted-proportion M-step for any N x T weight matrix.
  weighted_props <- function(w) {
    sweep(t(w) %*% onehot_y, 1, colSums(w), "/")
  }

  # N x (iT*C) estimating-equation matrix, column-major in (t, c): columns
  # 1:iT are category 1 for classes 1..iT, columns (iT+1):(2*iT) are
  # category 2, etc. -- matches matrix(theta, nrow = iT, ncol = C).
  score_matrix <- function(w, pi_hat) {
    J <- matrix(0, N, iT * C)
    for (c in seq_len(C)) {
      idx <- ((c - 1L) * iT + 1L):(c * iT)
      J[, idx] <- w * onehot_y[, c] - sweep(w, 2, pi_hat[, c], "*")
    }
    J
  }

  if (use.bch) {
    w.it <- bch_weight_matrix(w.is_cc, pwx)
    w_colsums <- colSums(w.it)
    if (any(w_colsums <= 0)) {
      stop(
        "BCH weights have non-positive column sums for at least one class. ",
        "The variance-covariance matrix will not be positive semi-definite. ",
        "Consider use.bch = FALSE.",
        call. = FALSE
      )
    }
    pi_hat <- weighted_props(w.it)
    log_pzx_bch <- log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300))

    three_step.score <- function(params) {
      score_matrix(w.it, matrix(params, nrow = iT, ncol = C))
    }

    H.3.inv <- diag(rep(1 / w_colsums, times = C))
    res <- list(
      par = as.vector(pi_hat),
      value = -sum(w.it * log_pzx_bch),
      convergence = 0L
    )
  } else {
    # EM: E-step (posterior responsibilities) / M-step (weighted
    # proportions, closed form). Initialize from a smoothed crosstab of the
    # modal Step-2 assignment against the observed category.
    init_class <- factor(max.col(w.is_cc), levels = seq_len(iT))
    init_cat <- factor(Y_cat, levels = seq_len(C))
    pi_hat <- unclass(table(init_class, init_cat)) + 0.5
    pi_hat <- pi_hat / rowSums(pi_hat)
    storage.mode(pi_hat) <- "double"

    # Per-record posteriors (see R/distal-ml.R) at a T x C probability matrix
    records_at <- function(pi_hat) {
      distal_records(
        log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300)), # N x T
        pi_s,
        pwx
      )
    }

    ll_prev <- -Inf
    for (iter in seq_len(em.maxIter)) {
      rec <- records_at(pi_hat)
      pi_hat_new <- weighted_props(distal_lambda(rec$R, w.is_cc))
      pi_hat_new <- pmax(pmin(pi_hat_new, 1 - 1e-10), 1e-10)
      pi_hat_new <- pi_hat_new / rowSums(pi_hat_new)

      ll_new <- distal_loglik(rec$logM, w.is_cc)
      if (iter > 1L && abs(ll_new - ll_prev) < covariate.tol) {
        pi_hat <- pi_hat_new
        if (verbose) {
          message(sprintf(
            "Multinomial ML EM converged in %d iterations.",
            iter
          ))
        }
        break
      }
      pi_hat <- pi_hat_new
      ll_prev <- ll_new
      if (iter == em.maxIter) {
        warning("Multinomial ML EM reached maximum iterations.")
      }
    }

    three_step.score <- function(params) {
      pi_hat_p <- matrix(params, nrow = iT, ncol = C)
      score_matrix(distal_lambda(records_at(pi_hat_p)$R, w.is_cc), pi_hat_p)
    }

    theta_hat <- as.vector(pi_hat)
    rec <- records_at(pi_hat)
    Jac <- distal_multinomial_jacobian(pi_hat, rec, w.is_cc, Y_cat)

    H.3.inv <- tryCatch(
      qr.solve(-Jac),
      error = function(e) {
        warning("Hessian inversion failed. SEs will be NA.")
        matrix(NA_real_, iT * C, iT * C)
      }
    )
    res <- list(
      par = theta_hat,
      value = -distal_loglik(rec$logM, w.is_cc),
      convergence = 0L
    )
  }

  list(
    res = res,
    H.3.inv = H.3.inv,
    three_step.score = three_step.score
  )
}

#' Step 3 (distal): estimate class-specific distal outcome parameters
#'
#' Estimates mu = (mu_1, ..., mu_T) for Gaussian (means), Poisson (log-rates),
#' or Binomial (logits) distal outcomes with either BCH (closed-form or Newton-Rhapson) or
#' ML EM. Returns the parameter estimates, the inverted Hessian H.3.inv, and
#' the case-wise score function three_step.score for sandwich variance
#' propagation in lca_vcov_distal.
#' @noRd
lca_step3.distal <- function(
  neg.ll,
  beta_init,
  iT,
  covariate.tol,
  use.bch = FALSE,
  Zo_cc = NULL,
  w.is_cc = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  family = "gaussian",
  p.zx = NULL,
  vPi = NULL,
  pi_mat = NULL,
  verbose = FALSE
) {
  N <- length(Zo_cc)
  #use covariate-adjusted pi if provided, otherwise flat vPi
  pi_s <- if (!is.null(pi_mat)) {
    pi_mat
  } else {
    matrix(vPi, ncol = iT, nrow = N, byrow = TRUE)
  }

  if (use.bch) {
    w.it <- bch_weight_matrix(w.is_cc, pwx) # N x T

    score_nt_bch <- function(mu) {
      if (family == "gaussian") {
        resid <- outer(Zo_cc, mu, "-")
        w.it * resid
      } else if (family == "poisson") {
        mu_val <- exp(mu)
        w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
      } else if (family == "binomial") {
        mu_val <- 1 / (1 + exp(-mu))
        w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
      }
    }

    w_colsums <- colSums(w.it)

    if (any(w_colsums < 0)) {
      stop(
        "BCH weights have negative column sums for at least one class. ",
        "The variance-covariance matrix will not be positive semi-definite. ",
        "Consider use.bch = FALSE."
      )
    }

    if (family == "gaussian") {
      beta <- colSums(w.it * Zo_cc) / w_colsums # closed-form weighted mean
      resid <- outer(Zo_cc, beta, "-")
      sigma2 <- sum(w.it * resid^2) / sum(w.it)

      three_step.score <- function(params) {
        mu <- params[1:iT]
        resid <- outer(Zo_cc, mu, "-")
        w.it * resid / sigma2
      }

      H.3.inv <- diag(sigma2 / w_colsums)

      res <- list(
        par = beta,
        value = neg.ll(beta),
        convergence = 0L,
        sigma2 = sigma2
      )
    } else {
      beta <- beta_init
      for (nr in seq_len(em.maxIter)) {
        grad_vec <- colSums(score_nt_bch(beta))

        if (family == "gaussian") {
          H_diag <- w_colsums
        } else if (family == "poisson") {
          H_diag <- exp(beta) * w_colsums
        } else if (family == "binomial") {
          mu_val <- 1 / (1 + exp(-beta))
          H_diag <- mu_val * (1 - mu_val) * w_colsums
        }

        direction <- grad_vec / H_diag

        step <- 1.0
        ll_cur <- -neg.ll(beta)
        for (ls in seq_len(20L)) {
          beta_new <- beta + step * direction
          ll_new <- tryCatch(-neg.ll(beta_new), error = function(e) -Inf)
          if (is.finite(ll_new) && ll_new > ll_cur) {
            break
          }
          step <- step * 0.5
        }

        delta <- step * direction
        beta <- beta + delta

        if (max(abs(delta)) < covariate.tol) {
          if (verbose) {
            message(sprintf("BCH NR converged in %d iterations.", nr))
          }
          break
        }
        if (nr == em.maxIter) warning("BCH NR reached maximum iterations.")
      }

      resid <- outer(Zo_cc, beta, "-")
      sigma2 <- sum(w.it * resid^2) / sum(w.it)

      three_step.score <- function(params) {
        mu <- params[1:iT]
        if (family == "gaussian") {
          resid <- outer(Zo_cc, mu, "-")
          w.it * resid
        } else if (family == "poisson") {
          mu_val <- exp(mu)
          w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
        } else if (family == "binomial") {
          mu_val <- 1 / (1 + exp(-mu))
          w.it * (outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val))
        }
      }

      H.3.inv <- tryCatch(
        diag(sigma2 / w_colsums),
        error = function(e) {
          warning("Hessian inversion failed. SEs will be NA.")
          matrix(NA_real_, iT, iT)
        }
      )
      res <- list(
        par = beta,
        value = neg.ll(beta),
        convergence = 0L,
        sigma2 = sigma2
      )
    }
  } else {
    if (is.null(p.zx)) {
      stop("p.zx must be provided for ML distal outcome estimation.")
    }

    # ML score helper
    # Three-step ML over the expanded data (R/distal-ml.R): the E-step
    # weights lambda_it = sum_s w_is R_ist average the per-record posteriors,
    # so proportional-assignment weights stay outside the log.
    #
    # For the gaussian family the common within-class variance sigma2 is
    # estimated jointly with the class means (Bakk, Tekle & Vermunt 2013:
    # normal distal outcome with constant error variance). Unlike ordinary
    # regression, sigma2 does not factor out of the mean estimates here: it
    # enters the posterior weights P(X = t | W_i, Zo_i), so fixing it would
    # bias the means whenever the true variance differs.
    gaussian <- family == "gaussian"
    theta_of <- function(mu, sigma2) if (gaussian) c(mu, sigma2) else mu
    records_at <- function(theta) distal_records(p.zx(theta), pi_s, pwx)
    derivs_at <- function(theta) distal_unit_derivs(theta, Zo_cc, iT, family)

    beta <- beta_init
    sigma2 <- if (gaussian) stats::var(Zo_cc) else NULL
    Z_long <- rep(Zo_cc, iT)
    X_long <- factor(rep(seq_len(iT), each = N))

    for (iter in seq_len(em.maxIter)) {
      theta <- theta_of(beta, sigma2)
      lambda <- distal_lambda(records_at(theta)$R, w.is_cc)

      # quasibinomial: same estimates as binomial, without glm's warning about
      # the fractional E-step weights
      fit <- glm(
        Z_long ~ X_long - 1,
        family = if (family == "binomial") stats::quasibinomial() else family,
        weights = as.vector(lambda)
      )
      beta_new <- coef(fit)
      sigma2_new <- if (gaussian) {
        sum(lambda * outer(Zo_cc, beta_new, "-")^2) / sum(lambda)
      } else {
        NULL
      }

      converged <- abs(neg.ll(theta_of(beta_new, sigma2_new)) - neg.ll(theta)) <
        covariate.tol
      beta <- beta_new
      sigma2 <- sigma2_new
      if (converged) {
        break
      }
      if (iter == em.maxIter) {
        warning("ML distal EM reached maximum iterations.")
      }
    }
    theta_hat <- theta_of(beta, sigma2)
    P <- length(theta_hat)

    three_step.score <- function(params) {
      distal_score(distal_lambda(records_at(params)$R, w.is_cc), derivs_at(params)$G)
    }

    H.3.inv <- tryCatch(
      qr.solve(distal_neg_hessian(records_at(theta_hat), w.is_cc, derivs_at(theta_hat))),
      error = function(e) {
        warning("Hessian inversion failed. SEs will be NA.")
        matrix(NA_real_, P, P)
      }
    )
    res <- list(
      par = beta,
      theta = theta_hat,
      value = neg.ll(theta_hat),
      convergence = 0L,
      sigma2 = sigma2
    )
  }

  return(list(
    res = res,
    H.3.inv = H.3.inv,
    three_step.score = three_step.score
  ))
}

#' Step 3 for three_step(): distal outcome model
#'
#' Class-specific distal outcome parameters (BCH or ML) with their variance
#' (Step-1 uncertainty propagated through `Sigma.1`, NULL for robust-only
#' standard errors) and joint-model fit statistics. When covariates are also
#' modeled (`cov`, the .fit_covariate() result), the class prior is the
#' fitted P(X = t | Zp) and the covariate-model uncertainty is propagated too.
#'
#' @return The distal component: a list with `three_step`,
#'   `three_step_vcov`, `three_step.llik`, `llik`, `AIC`, `BIC`, `npar`,
#'   `nobs`, and `sigma2` (gaussian) or `zo_levels` (multinomial).
#' @noRd
.fit_distal <- function(dat, s1, s2, Sigma.1, cov, n_classes, family, opts) {
  if (!(family %in% c("gaussian", "poisson", "binomial", "multinomial"))) {
    message(
      'Provided family is not one of "gaussian", "poisson", "binomial", nor "multinomial". Defaulting to family="gaussain".'
    )
  }
  prior <- .distal_prior(dat, s1$fit0, s2$dis, cov, n_classes, opts)
  cov_terms <- if (!is.null(cov)) {
    list(
      Sigma.3 = cov$Sigma.3,
      s3.par = cov$par,
      p.xz.cov = prior$p.xz.cov,
      Z_mat_cov = prior$Z_mat_cov
    )
  } else {
    list(Sigma.3 = NULL, s3.par = NULL, p.xz.cov = NULL, Z_mat_cov = NULL)
  }
  fit_fn <- if (family == "multinomial") .fit_distal_multinomial else .fit_distal_glm
  fit_fn(dat, s1$fit0, s2$dis, prior, Sigma.1, cov_terms, n_classes, family, opts)
}

#' Class prior and assignment weights for the distal model
#'
#' Without covariates: the Step-1 class sizes and the Step-2 assignments. With
#' covariates: the fitted P(X = t | Zp_i) as the prior, and assignment weights
#' and classification errors recomputed with that covariate-adjusted prior.
#'
#' @return list(pi_adj = n x T prior, res_adj = list(w.is, p.wx_mat),
#'   p.xz.cov = prior as a function of the covariate coefficients, Z_mat_cov =
#'   covariate design on the distal rows); the last two are NULL without
#'   covariates.
#' @noRd
.distal_prior <- function(dat, fit0, s2_for_dis, cov, n_classes, opts) {
  iT <- n_classes
  rows <- dat$keep_step3_Zo_in_Y
  flat_prior <- matrix(fit0$vPi, nrow = length(rows), ncol = iT, byrow = TRUE)

  if (is.null(cov)) {
    return(list(
      pi_adj = flat_prior,
      res_adj = list(w.is = s2_for_dis$w.is, p.wx_mat = s2_for_dis$p.wx_mat),
      p.xz.cov = NULL,
      Z_mat_cov = NULL
    ))
  }

  # Covariate design on the distal rows (which also have complete covariates)
  Z_mat_dis <- if (!is.null(dat$Z_mat) && length(dat$keep_step3_Zo) > 0L) {
    dat$Z_mat[dat$keep_step3_Zo_in_Z, , drop = FALSE]
  } else {
    NULL
  }

  p.xz_dis <- NULL
  if (!is.null(Z_mat_dis)) {
    p.xz_dis <- function(params) {
      eta_full <- cbind(0, Z_mat_dis %*% params)
      row_max <- apply(eta_full, 1L, max)
      exp_eta <- exp(eta_full - row_max)
      exp_eta / rowSums(exp_eta)
    }
    pi_adj <- p.xz_dis(matrix(cov$par, ncol = iT - 1))
  } else {
    pi_adj <- flat_prior
  }

  res_adj <- compute_pwx_adj(
    dat$Y.obs[rows, , drop = FALSE],
    fit0,
    dat$ivItemcat,
    if (!is.null(dat$mDesign)) dat$mDesign[rows, , drop = FALSE] else NULL,
    opts$use.modal.assignment,
    pi_adj = pi_adj
  )
  if (isTRUE(opts$uncorrected)) {
    res_adj$p.wx_mat <- diag(iT)
  }
  list(
    pi_adj = pi_adj,
    res_adj = res_adj,
    p.xz.cov = if (!is.null(p.xz_dis)) p.xz_dis else cov$p.xz,
    Z_mat_cov = if (!is.null(Z_mat_dis)) Z_mat_dis else dat$Z_mat
  )
}

#' Joint-model fit statistics of a distal outcome model
#'
#' Log-likelihood sum_i log sum_t P(X = t | Zp_i) P(Zo_i | X = t) P(Y_i | X = t)
#' with Step-1 parameters fixed. Free parameters: class sizes (or the
#' covariate logit coefficients, when the prior depends on covariates), item
#' parameters, and the `n_distal_params` outcome parameters.
#' @noRd
.distal_fit_stats <- function(log_pZo_t, pi_adj, dat, fit0, n_classes, n_distal_params) {
  iT <- n_classes
  rows <- dat$keep_step3_Zo_in_Y
  llik <- joint_log_lik_distal(
    Y = dat$Y.obs[rows, , drop = FALSE],
    mPhi = expand_Phi(fit0$mPhi, dat$ivItemcat),
    log_pZo_t = log_pZo_t,
    pi_mat = pi_adj,
    mDesign = if (!is.null(dat$mDesign)) dat$mDesign[rows, , drop = FALSE] else NULL
  )
  n_meas_params <- (if (!is.null(dat$Z_mat)) ncol(dat$Z_mat) * (iT - 1L) else iT - 1L) +
    sum(dat$ivItemcat - 1L) * iT
  k <- n_meas_params + n_distal_params
  n <- length(dat$keep_step3_Zo)
  list(llik = llik, AIC = -2 * llik + 2 * k, BIC = -2 * llik + k * log(n), npar = k, nobs = n)
}

#' Distal model for a nominal outcome (family = "multinomial")
#' @noRd
.fit_distal_multinomial <- function(
  dat,
  fit0,
  s2_for_dis,
  prior,
  Sigma.1,
  cov_terms,
  n_classes,
  family,
  opts
) {
  iT <- n_classes
  Y_cat <- dat$Zo_mat[, 1L] # already 1..C integer-coded
  zo_levels <- dat$zo_levels
  C <- length(zo_levels)

  s3.distal <- lca_step3.distal.multinomial(
    Y_cat = Y_cat,
    C = C,
    iT = iT,
    covariate.tol = opts$covariate.tol,
    use.bch = opts$use.bch,
    w.is_cc = prior$res_adj$w.is,
    pwx = prior$res_adj$p.wx_mat,
    em.maxIter = opts$em.maxIter,
    vPi = fit0$vPi,
    pi_mat = prior$pi_adj,
    verbose = opts$verbose
  )

  V <- lca_vcov_distal_multinomial(
    theta_hat = s3.distal$res$par,
    three_step.score = s3.distal$three_step.score,
    pi_adj = prior$pi_adj,
    w.is = prior$res_adj$w.is,
    p.wx_mat = prior$res_adj$p.wx_mat,
    Y_cat = Y_cat,
    C = C,
    H.3.inv = s3.distal$H.3.inv,
    Sigma.1 = Sigma.1,
    s2 = s2_for_dis,
    iT = iT,
    use.simple.cov = opts$use.simple.cov,
    use.bch = opts$use.bch,
    Sigma.3 = cov_terms$Sigma.3,
    s3.par = cov_terms$s3.par,
    p.xz.cov = cov_terms$p.xz.cov,
    Z_mat_cov = cov_terms$Z_mat_cov
  )

  class_labels <- paste0("C", seq_len(iT))
  pi_hat <- matrix(s3.distal$res$par, nrow = iT, ncol = C)
  dimnames(pi_hat) <- list(class_labels, zo_levels)
  param_labels <- as.vector(outer(class_labels, zo_levels, paste, sep = ":"))
  dimnames(V) <- list(param_labels, param_labels)

  stats <- .distal_fit_stats(
    log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300)),
    prior$pi_adj,
    dat,
    fit0,
    iT,
    iT * (C - 1L) # T x (C-1) free simplex parameters
  )
  c(
    list(
      three_step = pi_hat,
      three_step_vcov = V,
      three_step.llik = -s3.distal$res$value
    ),
    stats,
    list(zo_levels = zo_levels)
  )
}

#' Distal model for a gaussian, poisson, or binomial outcome
#' @noRd
.fit_distal_glm <- function(
  dat,
  fit0,
  s2_for_dis,
  prior,
  Sigma.1,
  cov_terms,
  n_classes,
  family,
  opts
) {
  iT <- n_classes
  z <- dat$Zo_mat[, 1L]
  pi_adj <- prior$pi_adj
  res_adj <- prior$res_adj

  # log f(z_i | X = t) as a function of the class parameters, and starting
  # values from a GLM of Zo on the modal assignment
  modal <- as.factor(max.col(res_adj$w.is))
  if (family == "poisson") {
    p.zx <- function(params) {
      log_mu <- params[1:iT]
      outer(z, log_mu, "*") - # N x T: z_i * log(mu_t)
        outer(rep(1, length(z)), exp(log_mu), "*") - # N x T: mu_t
        lgamma(z + 1L) # N x 1, recycled
    }
    beta_init <- coef(glm(z ~ -1 + modal, family = poisson()))
  } else if (family == "binomial") {
    # params: logit(mu_t)
    p.zx <- function(params) {
      mu <- 1 / (1 + exp(-params[1:iT]))
      outer(z, log(mu), "*") + outer(1 - z, log(1 - mu), "*")
    }
    beta_init <- coef(glm(z ~ -1 + modal, family = binomial()))
  } else {
    # gaussian. params: class means, optionally followed by the common
    # within-class variance sigma2 (defaults to 1 only when omitted, e.g. for
    # the closed-form BCH means, which do not depend on it)
    p.zx <- function(params) {
      mu <- params[1:iT]
      sigma2 <- if (length(params) > iT) params[iT + 1L] else 1
      resid <- outer(z, mu, "-")
      -0.5 * resid^2 / sigma2 - 0.5 * log(2 * pi * sigma2)
    }
    beta_init <- coef(lm(z ~ -1 + modal))
  }

  if (opts$use.bch) {
    w.it <- bch_weight_matrix(res_adj$w.is, res_adj$p.wx_mat)
    neg.ll <- function(params) -sum(w.it * p.zx(params))
  } else {
    # expanded-data log-likelihood, sum_i sum_s w_is log M_is (R/distal-ml.R)
    neg.ll <- function(params) {
      rec <- distal_records(p.zx(params), pi_adj, res_adj$p.wx_mat)
      -distal_loglik(rec$logM, res_adj$w.is)
    }
  }

  s3.distal <- lca_step3.distal(
    neg.ll = neg.ll,
    em.maxIter = opts$em.maxIter,
    pwx = res_adj$p.wx_mat,
    w.is_cc = res_adj$w.is,
    Zo_cc = z,
    use.bch = opts$use.bch,
    covariate.tol = opts$covariate.tol,
    iT = iT,
    beta_init = beta_init,
    family = family,
    p.zx = p.zx,
    vPi = fit0$vPi,
    pi_mat = pi_adj
  )

  # Full Step-3 parameter vector: class parameters, plus sigma2 for the
  # gaussian family (estimated jointly under ML; the BCH means do not depend
  # on it and it is estimated from the weighted residuals).
  gaussian_ml <- family == "gaussian" && !opts$use.bch
  theta_zx <- if (family == "gaussian") {
    c(s3.distal$res$par, s3.distal$res$sigma2)
  } else {
    s3.distal$res$par
  }

  V <- lca_vcov_distal(
    mu_hat = if (gaussian_ml) theta_zx else s3.distal$res$par,
    three_step.score = s3.distal$three_step.score,
    pi_adj = pi_adj,
    w.is = res_adj$w.is,
    p.wx_mat = res_adj$p.wx_mat,
    p.zx = p.zx,
    family = family,
    H.3.inv = s3.distal$H.3.inv,
    Sigma.1 = Sigma.1,
    s2 = s2_for_dis,
    Sigma.3 = cov_terms$Sigma.3,
    s3.par = cov_terms$s3.par,
    p.xz.cov = cov_terms$p.xz.cov,
    Z_mat_cov = cov_terms$Z_mat_cov,
    iT = iT,
    use.simple.cov = opts$use.simple.cov,
    use.bch = opts$use.bch,
    unit_scores = function(theta) distal_unit_derivs(theta, z, iT, family)$G
  )

  # Report the class parameters; keep sigma2 (and its SE under ML) separately.
  sigma2_hat <- NULL
  if (family == "gaussian") {
    sigma2_hat <- c(
      estimate = s3.distal$res$sigma2,
      se = if (gaussian_ml) sqrt(V["sigma2", "sigma2"]) else NA_real_
    )
    V <- V[seq_len(iT), seq_len(iT), drop = FALSE]
  }

  par <- s3.distal$res$par
  names(par) <- paste0("mu_C", seq_len(iT))

  stats <- .distal_fit_stats(
    p.zx(theta_zx),
    pi_adj,
    dat,
    fit0,
    iT,
    iT + as.integer(family == "gaussian") # + sigma2
  )
  c(
    list(
      three_step = par,
      three_step_vcov = V,
      three_step.llik = -neg.ll(theta_zx)
    ),
    stats,
    list(sigma2 = sigma2_hat)
  )
}
