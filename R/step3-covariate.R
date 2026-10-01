# tseLCA/R/step3-covariate.R
#
# Step 3 for explanatory covariates: bias-adjusted multinomial logit of
# latent class on covariates (BCH or ML).

#' Step 3 (covariate): estimate multinomial logit gamma with either BCH or ML EM
#'
#' Optimizes the (Q+1) x (T-1) coefficient matrix gamma for P(X=t|Z_i) with
#' Newton-Raphson (BCH) or EM with an inner NR M-step (ML). Returns the
#' parameter vector, the inverted Hessian H.3.inv (or NA matrix on failure),
#' used by lca_vcov for sandwich variance propagation.
#' @noRd
lca_step3 <- function(
  neg.ll,
  gamma_init,
  Q,
  iT,
  covariate.tol,
  use.bch = FALSE,
  gradient = NULL,
  Z_mat_cc = NULL,
  w.is_cc = NULL,
  p.xz = NULL,
  pwx = NULL,
  em.maxIter = 200L,
  verbose = FALSE,
  correct.spec = FALSE
) {
  N <- nrow(Z_mat_cc)
  beta <- matrix(gamma_init, nrow = Q, ncol = iT - 1)
  ll_prev <- -neg.ll(c(beta))
  H <- NULL
  # print(ll_prev)
  if (use.bch) {
    w.it <- bch_weight_matrix(w.is_cc, pwx) # N x T
    w.it_plus <- rowSums(w.it)

    for (nr in seq_len(em.maxIter)) {
      # print(beta)
      # print(-neg.ll(c(beta)))
      grad_vec <- -gradient(c(beta)) # gradient of pos. ll: Q*(T-1) vector

      H <- matrix(0, Q * (iT - 1), Q * (iT - 1))
      pi_ <- p.xz(beta)
      for (k in seq_len(iT - 1)) {
        for (l in k:(iT - 1)) {
          w_kl <- w.it_plus * pi_[, k + 1L] * ((k == l) - pi_[, l + 1L])
          idx_k <- ((k - 1) * Q + 1):(k * Q)
          idx_l <- ((l - 1) * Q + 1):(l * Q)
          block <- -t(Z_mat_cc) %*% (w_kl * Z_mat_cc)
          H[idx_k, idx_l] <- block
          if (k != l) H[idx_l, idx_k] <- t(block) # Clairaut: H symmetric
        }
      }

      if (
        nr >= max(em.maxIter / 5, 1) && #Wait a little bit before testing for PSD
          inherits(
            tryCatch(chol(H), error = function(e) e),
            "error"
          )
      ) {
        stop(
          sprintf(
            "BCH Newton-Raphson failed after %d iterations: Hessian is not positive semi-definite. ",
            nr
          ),
          "This typically occurs under low class separation. ",
          "Try the ML estimator (use.bch = FALSE), or increase em.maxIter.",
          call. = FALSE
        )
      }

      direction <- tryCatch(
        qr.solve(-H, grad_vec),
        error = function(e) rep(0, Q * (iT - 1))
      )

      alpha <- 1
      current_ll <- -neg.ll(c(beta))

      beta_vec <- c(beta)

      while (alpha > (covariate.tol / 2)) {
        trial <- beta_vec + alpha * direction
        trial_ll <- tryCatch(-neg.ll(trial), error = function(e) -Inf)
        if (is.finite(trial_ll) && trial_ll > current_ll) {
          break
        }
        alpha <- alpha / 2
      }

      delta <- alpha * direction
      #print(max(abs(delta)))
      beta <- beta + matrix(delta, nrow = Q, ncol = iT - 1)
      if (max(abs(delta)) < covariate.tol) {
        if (verbose) {
          message(sprintf("BCH NR converged in %d iterations.", nr))
        }
        break
      }
      if (nr == em.maxIter) stop("BCH NR reached maximum iterations.")
    }
    #print(H)
  } else {
    for (iter in seq_len(em.maxIter)) {
      #print(beta)
      #print(ll_prev)

      # E-step ######################################################################
      p <- p.xz(beta)

      q <- p %*% t(pwx)

      gamma <- matrix(0, nrow = N, ncol = iT)
      for (t in seq_len(iT)) {
        gamma[, t] <- p[, t] * rowSums(w.is_cc * outer(rep(1, N), pwx[, t]) / q)
      }
      ###############################################################################

      #M-step ########################################################################
      gamma_plus <- rowSums(gamma)
      Gamma_nr <- gamma[, -1, drop = FALSE]

      for (nr in seq_len(10L)) {
        p_nr <- p.xz(beta)
        p_nr1 <- p_nr[, -1, drop = FALSE]

        grad <- t(Z_mat_cc) %*% (Gamma_nr - p_nr1 * gamma_plus)

        H <- matrix(0, Q * (iT - 1), Q * (iT - 1))
        for (k in seq_len(iT - 1)) {
          for (l in k:(iT - 1)) {
            w_kl <- gamma_plus * p_nr1[, k] * ((k == l) - p_nr1[, l])
            idx_k <- ((k - 1) * Q + 1):(k * Q)
            idx_l <- ((l - 1) * Q + 1):(l * Q)
            block <- -t(Z_mat_cc) %*% (w_kl * Z_mat_cc)
            H[idx_k, idx_l] <- block
            # Clairaut: H symmetric
            if (k != l) H[idx_l, idx_k] <- t(block)
          }
        }

        delta <- tryCatch(
          qr.solve(-H, as.vector(grad)),
          error = function(e) rep(0, Q * (iT - 1))
        )
        beta <- beta + matrix(delta, nrow = Q, ncol = iT - 1)
        if (max(abs(delta)) < covariate.tol) break
      }
      if (nr == 10L && max(abs(delta)) >= covariate.tol) {
        warning(sprintf(
          "M-step in EM algorithm did not converge at iteration %d",
          iter
        ))
      }

      #Alternatively, fit a multinomial logistic regression model ##########################
      # class_exp <- rep(seq_len(iT), each = N)
      # Z_exp <- Z_mat_cc[rep(seq_len(N), iT), , drop = FALSE]
      # w_exp <- as.vector(gamma)

      # fit <- nnet::multinom(
      #   class_exp ~ Z_exp - 1,
      #   weights = w_exp,
      #   trace = FALSE,
      #   maxit = 500L
      # )

      # beta <- t(coef(fit))
      ######################################################################################

      ll_curr <- -neg.ll(c(beta))
      if (abs(ll_curr - ll_prev) < covariate.tol && iter > 1L) {
        if (verbose) {
          message(sprintf("EM converged in %d iterations.", iter))
        }
        break
      }
      ll_prev <- ll_curr

      if (iter == em.maxIter) {
        warning("EM reached maximum iterations without converging.")
      }
    }
  }
  res <- list(par = c(beta), value = neg.ll(c(beta)), convergence = 0L)

  H.3.inv <- tryCatch(
    {
      if (!use.bch) {
        if (correct.spec) {
          matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
        } else {
          # -- Analytic observed-data Hessian of neg.ll (checked with sympy)--------------------------
          # neg.ll = -sum_i sum_s w_{is} * log(q_{is})
          # q_{is} = sum_t pi_{it}(beta) * pwx[s,t]
          # r_{is} = w_{is} / q_{is}
          #
          # H_{(q,k),(p,l)} = sum_i z_{iq}*z_{ip} * [
          #   pi_{i,k+1}*(I(k==l)-pi_{i,l+1}) * F_k
          # - pi_{i,k+1} * pi_{i,l+1} * G_{kl}
          # ]
          # where:
          #   F_k  = sum_s w_{is}*pwx[s,k+1]/q_{is} - sum_s w_{is}
          #   G_{kl} = sum_s w_{is}*pwx[s,k+1]*(pwx[s,l+1]-q_{is})/q_{is}^2

          p_ <- p.xz(beta) # N x T
          q_mat <- p_ %*% t(pwx) # N x T: q[i,s]
          r_mat <- w.is_cc / q_mat # N x T: r[i,s]

          # F_k for each non-reference class k (N x (T-1))
          F_mat <- matrix(0, N, iT - 1L)
          for (k in seq_len(iT - 1L)) {
            F_mat[, k] <- rowSums(r_mat * pwx[, k + 1L][col(r_mat)]) -
              rowSums(w.is_cc)
          }

          # G_{kl} for each pair (k,l)
          H_obs <- matrix(0, Q * (iT - 1L), Q * (iT - 1L))
          for (k in seq_len(iT - 1L)) {
            for (l in k:(iT - 1L)) {
              idx_k <- ((k - 1L) * Q + 1L):(k * Q)
              idx_l <- ((l - 1L) * Q + 1L):(l * Q)
              tA <- p_[, k + 1L] * ((k == l) - p_[, l + 1L]) * F_mat[, k]
              G_kl <- rowSums(
                w.is_cc *
                  pwx[, k + 1L][col(w.is_cc)] *
                  (pwx[, l + 1L][col(w.is_cc)] - q_mat) /
                  q_mat^2
              )
              tB <- -p_[, k + 1L] * p_[, l + 1L] * G_kl

              block <- -t(Z_mat_cc) %*% ((tA + tB) * Z_mat_cc) # Q x Q
              H_obs[idx_k, idx_l] <- block
              if (k != l) {
                H_obs[idx_l, idx_k] <- t(block)
              } # Clairaut: H is symmetric
            }
          }
          qr.solve(H_obs)
        }
      } else {
        qr.solve(-H)
      }
    },
    # the caller falls back to the outer product of the scores and warns
    error = function(e) matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
  )
  return(list(res = res, H.3.inv = H.3.inv))
}


# -- Variance estimation (Bakk et al., 2014) ----------------------------------

#' Step 3 for three_step(): covariate model
#'
#' Bias-adjusted multinomial logit of class on the covariates (BCH or ML),
#' its variance with Step-1 uncertainty propagation (`Sigma.1`, NULL for
#' robust-only standard errors), joint-model fit statistics, the optional
#' two-step variance, and the covariate-adjusted entropy R^2.
#'
#' @return list(fit = the tseLCA_covariate object, s1 = the measurement
#'   model (its `$fitZ` may be filled in by `get.twostep.vcov`), par = the
#'   Step-3 coefficient vector, Sigma.3 = its variance, p.xz = the class-prior
#'   function of the coefficients).
#' @noRd
.fit_covariate <- function(dat, s1, s2, Sigma.1, n_classes, opts) {
  iT <- n_classes
  Z_mat <- dat$Z_mat
  Q <- ncol(Z_mat)
  s2_for_cov <- s2$cov
  fit0 <- s1$fit0
  fitZ <- s1$fitZ

  p.xz <- function(params) {
    eta_full <- cbind(0, Z_mat %*% params)
    row_max <- apply(eta_full, 1, max)
    exp_eta <- exp(eta_full - row_max)
    exp_eta / rowSums(exp_eta)
  }

  if (opts$use.bch) {
    w.it <- bch_weight_matrix(s2_for_cov$w.is, s2_for_cov$p.wx_mat)

    neg.ll <- function(params) {
      beta.cur <- matrix(params, ncol = iT - 1)
      -sum(w.it * log(pmax(p.xz(beta.cur), 1e-6)))
    }

    # (pwx is unused by BCH; kept so both estimators share one signature)
    three_step.grad <- function(params, pwx = s2_for_cov$p.wx_mat) {
      beta.cur <- matrix(params, ncol = iT - 1)
      pi_ <- p.xz(beta.cur)
      resid <- w.it[, -1L, drop = FALSE] -
        pi_[, -1L, drop = FALSE] * rowSums(w.it)
      -as.vector(t(Z_mat) %*% resid)
    }

    three_step.score <- function(params, pwx = s2_for_cov$p.wx_mat) {
      beta.cur <- matrix(params, ncol = iT - 1)
      pi_ <- p.xz(beta.cur)
      resid <- w.it[, -1L, drop = FALSE] -
        pi_[, -1L, drop = FALSE] * rowSums(w.it)
      resid[, rep(seq_len(iT - 1L), each = Q)] *
        Z_mat[, rep(seq_len(Q), iT - 1L)]
    }
  } else {
    three_step.ll <- function(params, pwx = s2_for_cov$p.wx_mat) {
      probs <- p.xz(matrix(params, ncol = iT - 1))
      rowSums(s2_for_cov$w.is * log(probs %*% t(pwx)))
    }

    three_step.grad <- function(params, pwx = s2_for_cov$p.wx_mat) {
      beta <- matrix(params, ncol = iT - 1)
      p <- p.xz(beta)
      q <- p %*% t(pwx)
      r <- s2_for_cov$w.is / q
      grad <- matrix(0, nrow = iT - 1, ncol = Q)
      for (k in seq_len(iT - 1L)) {
        score_i <- p[, k + 1L] *
          (r %*% pwx[, k + 1L] - rowSums(s2_for_cov$w.is))
        grad[k, ] <- t(Z_mat) %*% score_i
      }
      as.vector(t(grad))
    }

    three_step.score <- function(params, pwx = s2_for_cov$p.wx_mat) {
      beta <- matrix(params, ncol = iT - 1)
      p <- p.xz(beta)
      q <- p %*% t(pwx)
      r <- s2_for_cov$w.is / q
      score_ik <- matrix(0, nrow = nrow(Z_mat), ncol = iT - 1)
      for (k in seq_len(iT - 1L)) {
        score_ik[, k] <- p[, k + 1L] *
          (r %*% pwx[, k + 1L] - rowSums(s2_for_cov$w.is))
      }
      score_ik[, rep(seq_len(iT - 1L), each = Q)] *
        Z_mat[, rep(seq_len(Q), iT - 1L)]
    }

    neg.ll <- function(params) -sum(three_step.ll(params))
  }

  # exact matching: `opts$start` would partially match `opts$startval`
  gamma_start <- opts[["gamma_start", exact = TRUE]]
  gamma_init <- if (!is.null(gamma_start)) {
    c(gamma_start)
  } else if (!is.null(fitZ$mGamma) && opts$use.two.step) {
    c(fitZ$mGamma)
  } else {
    rep(0, Q * (iT - 1))
  }

  s3 <- lca_step3(
    neg.ll,
    gamma_init,
    Q,
    iT,
    opts$covariate.tol,
    gradient = three_step.grad,
    use.bch = opts$use.bch,
    Z_mat_cc = Z_mat,
    w.is_cc = s2_for_cov$w.is,
    p.xz = p.xz,
    pwx = s2_for_cov$p.wx_mat,
    em.maxIter = opts$em.maxIter,
    verbose = opts$verbose,
    correct.spec = opts$correct.spec
  )
  use_opg <- opts$correct.spec && !opts$use.bch
  if (!use_opg && (is.null(s3$H.3.inv) || !all(is.finite(s3$H.3.inv)))) {
    warning(
      "The Step-3 Hessian could not be inverted; the information matrix is ",
      "estimated by the outer product of the case-wise scores, which assumes ",
      "a correctly specified Step-3 model.",
      call. = FALSE
    )
    use_opg <- TRUE
  }
  if (use_opg) {
    s3$H.3.inv <- qr.solve(crossprod(three_step.score(s3$res$par)))
  }

  coefs <- matrix(s3$res$par, ncol = iT - 1)
  ref_idx <- parse_rebase(opts$rebase, iT)
  colnames(coefs) <- paste0("C", seq_len(iT)[-ref_idx])
  rownames(coefs) <- colnames(Z_mat)

  # -- Variance ---------------------------------------------------------------
  Sigma.3 <- lca_vcov(
    coefs = coefs,
    three_step.score = three_step.score,
    H.3.inv = s3$H.3.inv,
    Sigma.1 = Sigma.1,
    J.2 = s2_for_cov$J.2,
    p.wx_mat = s2_for_cov$p.wx_mat,
    w.is = s2_for_cov$w.is,
    Z_mat = Z_mat,
    n_classes = n_classes,
    p.xz = p.xz,
    s2 = s2_for_cov,
    use.simple.cov = opts$use.simple.cov || opts$use.bch
  )

  # -- Model fit --------------------------------------------------------------
  rows <- dat$keep_step3_Z_in_Y
  Y_cc <- dat$Y.obs[rows, , drop = FALSE]
  mDes_cc <- if (!is.null(dat$mDesign)) dat$mDesign[rows, , drop = FALSE] else NULL
  total.llik <- joint_log_lik(
    Y_cc,
    Z_mat,
    expand_Phi(fit0$mPhi, dat$ivItemcat),
    coefs,
    mDes_cc
  )
  # Free parameters of the joint model: item-response log-ratios plus the
  # multinomial-logit coefficients (whose intercepts replace class sizes).
  total.k <- iT * sum(dat$ivItemcat - 1L) + Q * (iT - 1L)

  # -- Optional two-step variance ---------------------------------------------
  tsv <- .two_step_vcov(fitZ, dat, n_classes, opts)
  if (is.null(fitZ) && !is.null(tsv$fitZ)) {
    fitZ <- tsv$fitZ
    s1$fitZ <- tsv$fitZ
  }

  # -- Covariate-adjusted entropy R^2 -----------------------------------------
  entropy.R2 <- .covariate_entropy_R2(
    p.xz(matrix(s3$res$par, ncol = iT - 1L)),
    Y_cc,
    fit0,
    dat$ivItemcat,
    mDes_cc
  )

  fit <- structure(
    list(
      measurement_model = s1,
      two_step = if (!is.null(fitZ)) fitZ$mGamma else NULL,
      two_step_vcov = tsv$vcov,
      three_step = coefs,
      three_step_vcov = Sigma.3,
      three_step.llik = -s3$res$value,
      neg.ll = neg.ll,
      llik = total.llik,
      AIC = -2 * total.llik + 2 * total.k,
      BIC = -2 * total.llik + total.k * log(nrow(Y_cc)),
      npar = total.k,
      nobs = nrow(Y_cc),
      n_classes = iT,
      estimator = .estimator_label(opts),
      entropy.R2 = entropy.R2,
      posteriors = s2$all$p.xy,
      classifications = max.col(s2$all$p.xy)
    ),
    class = c("tseLCA_covariate", "tseLCA_structural", "tseLCA")
  )
  list(fit = fit, s1 = s1, par = s3$res$par, Sigma.3 = Sigma.3, p.xz = p.xz)
}

#' Two-step variance (multilevLCA's bias-corrected Varmat) for the covariate
#' model
#'
#' A Varmat_cor already attached to `fitZ` is always returned: at
#' fitZ$Varmat_cor (plain multiLCA output) or fitZ$raw_fit$Varmat_cor
#' (fitZ_from_multiLCA() output). Otherwise, with `get.twostep.vcov = TRUE`,
#' the two-step model is re-estimated with multiLCA(); fitZ_from_fit0() output
#' carries no Varmat_cor.
#'
#' @return list(vcov = named variance matrix or NULL, fitZ = the
#'   re-estimated two-step fit, if any).
#' @noRd
.two_step_vcov <- function(fitZ, dat, n_classes, opts) {
  extract_varmat <- function(fZ) {
    if (is.null(fZ)) {
      return(NULL)
    }
    if (!is.null(fZ$Varmat_cor)) {
      return(fZ$Varmat_cor)
    }
    if (!is.null(fZ$raw_fit$Varmat_cor)) {
      return(fZ$raw_fit$Varmat_cor)
    }
    if (!is.null(fZ$raw_fit$SEs_cor_gamma)) {
      return(diag(as.vector(fZ$raw_fit$SEs_cor_gamma)^2))
    }
    NULL
  }
  name_varmat <- function(V, fZ) {
    if (is.null(V) || is.null(fZ$mGamma)) {
      return(V)
    }
    nms <- as.vector(outer(rownames(fZ$mGamma), colnames(fZ$mGamma), paste, sep = ":"))
    dimnames(V) <- list(nms, nms)
    V
  }

  existing <- extract_varmat(fitZ)
  if (!is.null(existing)) {
    return(list(vcov = name_varmat(existing, fitZ), fitZ = NULL))
  }
  if (!opts$get.twostep.vcov) {
    return(list(vcov = NULL, fitZ = NULL))
  }
  fZ_ml <- fitZ_from_multiLCA(
    data = dat$data,
    Y.names = dat$Y.names,
    n_classes = n_classes,
    Zp.names = dat$Zp.names,
    maxIter.measurement = opts$maxIter.measurement,
    measurement.tol = opts$measurement.tol,
    covariate.tol = opts$covariate.tol,
    iter.measurement = opts$iter.measurement,
    R2.threshold = opts$R2.threshold,
    incomplete = opts$incomplete,
    rebase = opts$rebase,
    startval = opts$startval,
    n_init = opts$n_init,
    verbose = opts$verbose
  )
  raw <- extract_varmat(fZ_ml)
  if (is.null(raw)) {
    warning("get.twostep.vcov: neither Varmat_cor nor SEs_cor_gamma found.")
  }
  list(vcov = name_varmat(raw, fZ_ml), fitZ = fZ_ml)
}

#' Covariate-adjusted entropy R^2
#'
#' How much the indicators reduce classification uncertainty beyond what the
#' covariates already explain: (H(X|Z) - H(X|Y,Z)) / H(X|Z), with H(X|Z) the
#' average entropy of the fitted class priors P(X|Z_i) and H(X|Y,Z) that of
#' the covariate-adjusted posteriors (soft assignment).
#' @noRd
.covariate_entropy_R2 <- function(pi_adj, Y_cc, fit0, ivItemcat, mDes_cc) {
  h <- function(p) {
    p <- p[p > sqrt(.Machine$double.eps)]
    -sum(p * log(p))
  }
  error_prior <- mean(apply(pi_adj, 1L, h))
  adj_res <- compute_pwx_adj(
    Y.obs = Y_cc,
    fit0 = fit0,
    ivItemcat = ivItemcat,
    mDesign = mDes_cc,
    use.modal.assignment = FALSE,
    pi_adj = pi_adj
  )
  error_post <- mean(apply(adj_res$post, 1L, h))
  if (error_prior > 1e-8) {
    (error_prior - error_post) / error_prior
  } else {
    1.0 # covariates already explain all class membership
  }
}
