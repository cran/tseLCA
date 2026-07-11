#' One-hot expand an integer response matrix
#'
#' Converts an N x H matrix of 0-based integer category values into an
#' N x sum(ivItemcat) binary indicator matrix, one column per category per item.
#' @noRd
expand_Y <- function(mY_int, ivItemcat) {
  # mY_int: N x H matrix of integer category values (0-based)
  # ivItemcat: length-H vector of number of categories per item
  N <- nrow(mY_int)
  H <- ncol(mY_int)
  out <- matrix(0, N, sum(ivItemcat))
  col_start <- 1L
  for (h in seq_len(H)) {
    K_h <- ivItemcat[h]
    for (k in seq_len(K_h)) {
      out[, col_start + k - 1L] <- as.integer(mY_int[, h] == (k - 1L))
    }
    col_start <- col_start + K_h
  }
  out
}

#' Expand a compact mPhi to a full item-probability matrix
#'
#' Converts multilevLCA's storage convention (one row per dichotomous item,
#' K rows per polytomous item) into a sum(ivItemcat) x T expanded matrix where
#' each item block contains all K category probabilities including the reference.
#' @noRd
expand_Phi <- function(phi_mat, ivItemcat) {
  dichotomous <- ivItemcat == 2L
  result <- vector("list", length(ivItemcat))
  h_phi <- 1L
  for (h in seq_along(ivItemcat)) {
    if (dichotomous[h]) {
      result[[h]] <- rbind(1 - phi_mat[h_phi, ], phi_mat[h_phi, ])
      h_phi <- h_phi + 1L
    } else {
      K_h <- ivItemcat[h]
      result[[h]] <- phi_mat[h_phi:(h_phi + K_h - 1L), , drop = FALSE]
      h_phi <- h_phi + K_h
    }
  }
  do.call(rbind, result)
}

#' Expand a free-parameter phi matrix to full category probabilities
#'
#' Inverse of the simplex constraint: given (K-1) free rows per polytomous item
#' and 1 row per dichotomous item, prepends the reference P(Y=0|C) row for each
#' polytomous item so that the result aligns with expand_Y output.
#' @noRd
expand_Phi_free <- function(phi_free, ivItemcat) {
  result <- vector("list", length(ivItemcat))
  h_phi <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    if (K_h == 2L) {
      result[[h]] <- phi_free[h_phi, , drop = FALSE]
      h_phi <- h_phi + 1L
    } else {
      rows_h <- phi_free[h_phi:(h_phi + K_h - 2L), , drop = FALSE]
      result[[h]] <- rbind(1 - colSums(rows_h), rows_h)
      h_phi <- h_phi + K_h - 1L
    }
  }
  do.call(rbind, result)
}

#' Per-observation class log-likelihood matrix
#'
#' Returns an N x T matrix where entry `[i, t]` is the conditional log-likelihood
#' log P(Y_i | X=t) under the expanded item-probability matrix mPhi.
#' mDesign masks missing indicators (0 = missing, 1 = observed).
#' @noRd
log_lik_matrix <- function(Y, mPhi, mDesign = NULL) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }
  (mDesign * Y) %*% log(mPhi)
}

#' Joint observed-data log-likelihood with covariates
#'
#' Computes sum_i log P(Y_i, Z_i) = sum_i log sum_t P(Y_i|X=t) P(X=t|Z_i)
#' under a multinomial logit structural model with coefficient matrix gamma.coefs
#' (Q x (T-1), reference class absorbed into the intercept column of Z).
#' @noRd
joint_log_lik <- function(Y, Z, mPhi, gamma.coefs, mDesign = NULL) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }

  log_P_Y_given_X <- log_lik_matrix(Y, mPhi, mDesign)

  eta <- Z %*% gamma.coefs
  eta_full <- cbind(0, eta)

  row_maxes_eta <- apply(eta_full, 1, max)
  log_denom_eta <- row_maxes_eta + log(rowSums(exp(eta_full - row_maxes_eta)))
  log_P_X_given_Z <- eta_full - log_denom_eta

  log_joint_prob <- log_P_Y_given_X + log_P_X_given_Z

  row_maxes_joint <- apply(log_joint_prob, 1, max)
  log_marg_prob <- row_maxes_joint +
    log(rowSums(exp(log_joint_prob - row_maxes_joint)))

  sum(log_marg_prob)
}

#' Joint log-likelihood for the distal outcome model
#'
#' Computes sum_i log( sum_t P(X=t|Zp_i) * P(Zo_i|X=t) * P(Y_i|X=t) )
#' When Zp is absent, P(X=t|Zp) = vPi (flat prevalences).
#'
#' @param Y       N x K_total expanded one-hot response matrix.
#' @param Zo      Length-N distal outcome vector.
#' @param mPhi    K_total x T expanded item-response probability matrix.
#' @param p.zx    N x T matrix of log-densities log P(Zo_i|X=t).
#' @param pi_mat  N x T matrix of class priors P(X=t|Zp_i). If NULL, uses
#'   flat prevalences from the row means of p.zx (not used; vPi supplied).
#' @param vPi     Length-T flat prevalences, used when pi_mat is NULL.
#' @param mDesign N x K_total design matrix (NULL for complete data).
#' @noRd
joint_log_lik_distal <- function(
  Y,
  mPhi,
  log_pZo_t,
  pi_mat = NULL,
  vPi = NULL,
  mDesign = NULL
) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y), ncol(Y))
  }

  # log P(Y_i | X=t): N x T
  log_P_Y_t <- log_lik_matrix(Y, mPhi, mDesign)

  # log P(X=t | Zp_i): N x T
  if (!is.null(pi_mat)) {
    log_P_X_t <- log(pmax(pi_mat, 1e-300))
  } else {
    # flat prevalences: broadcast vPi across rows
    log_P_X_t <- matrix(
      log(pmax(vPi, 1e-300)),
      nrow(Y),
      length(vPi),
      byrow = TRUE
    )
  }

  # log P(Zo_i | X=t): N x T  (passed in as log_pZo_t)
  log_joint <- log_P_X_t + log_pZo_t + log_P_Y_t # N x T

  row_max <- apply(log_joint, 1L, max)
  log_marg <- row_max + log(rowSums(exp(log_joint - row_max)))
  sum(log_marg)
}

#' Compute posterior class probabilities from the unconstrained theta1 vector
#'
#' Reconstructs vPi and phi from the stacked parameter vector theta1 =
#' `c(vPi[-1], phi_free)` used internally by lca_step2, then returns the
#' N x T soft posterior matrix P(X=t|Y_i).
#' @noRd
compute_posteriors <- function(Y, mDesign, theta1, ivItemcat, iT) {
  vPi_free <- theta1[1:(iT - 1L)]
  vPi <- c(1 - sum(vPi_free), vPi_free)
  n_free <- sum(ifelse(ivItemcat == 2L, 1L, ivItemcat - 1L))
  phi_free <- matrix(
    theta1[iT:(iT + n_free * iT - 1L)],
    nrow = n_free,
    ncol = iT
  )

  mPhi <- expand_Phi(expand_Phi_free(phi_free, ivItemcat), ivItemcat)

  log_joint <- sweep(log_lik_matrix(Y, mPhi, mDesign), 2, log(vPi), "+")
  row_max <- apply(log_joint, 1, max)
  log_denom <- row_max + log(rowSums(exp(log_joint - row_max)))
  exp(log_joint - log_denom)
}

#' Compute classification-error matrix with optional covariate-adjusted prior
#'
#' Returns posteriors, modal/soft assignments (w.is), and the T x T
#' classification-error probability matrix p.wx_mat = P(W=s|X=t).
#' When pi_adj (N x T) is supplied, uses person-specific class priors from the
#' covariate model; otherwise falls back to the flat vPi from fit0.
#' @noRd
compute_pwx_adj <- function(
  Y.obs,
  fit0,
  ivItemcat,
  mDesign = NULL,
  use.modal.assignment = TRUE,
  pi_adj = NULL # N x T covariate-adjusted class probs, or NULL for flat vPi
) {
  N <- nrow(Y.obs)
  iT <- ncol(fit0$mPhi)

  mPhi_exp <- expand_Phi(fit0$mPhi, ivItemcat)
  phi_clamped <- pmax(pmin(mPhi_exp, 1 - 1e-10), 1e-10)

  log_p_it <- if (is.null(mDesign)) {
    Y.obs %*% log(phi_clamped)
  } else {
    (mDesign * Y.obs) %*% log(phi_clamped)
  }

  # use adjusted or flat priors
  if (!is.null(pi_adj)) {
    log_prior <- log(pi_adj) # N x T, person-specific
  } else {
    log_prior <- matrix(log(fit0$vPi), nrow = N, ncol = iT, byrow = TRUE)
  }

  log_joint <- log_p_it + log_prior # N x T
  row_max <- apply(log_joint, 1, max)
  post <- exp(log_joint - row_max - log(rowSums(exp(log_joint - row_max)))) # N x T posteriors

  w.is <- if (use.modal.assignment) {
    w <- matrix(0L, N, iT)
    w[cbind(seq_len(N), max.col(post))] <- 1L
    w
  } else {
    post
  }

  p.wx_joint <- (t(w.is) %*% post) / N
  p.wx_mat <- sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")

  list(
    post = post,
    w.is = w.is,
    p.wx_mat = p.wx_mat
  )
}

# -- Step 2: Posteriors and classification-error matrix -----------------------

#' Step 2: posteriors, classification-error matrix, and Jacobian closure
#'
#' Computes all Step-2 quantities needed for Step 3 and variance propagation:
#' theta1 and theta2 (constrained and unconstrained parameterizations of the
#' classification-error matrix), w.is (modal or soft assignments), p.wx_mat,
#' and optionally a closure compute_J_unc for the analytic Jacobian
#' d theta2 / d u used in the measurement-uncertainty correction.
#' Returns NULL for compute_J_unc when use.simple.cov = TRUE.
#' @noRd
lca_step2 <- function(
  Y.obs,
  fit0,
  n_classes,
  use.modal.assignment,
  boundary.tol,
  use.simple.cov,
  ivItemcat,
  mDesign = NULL
) {
  if (is.null(mDesign)) {
    mDesign <- matrix(1L, nrow(Y.obs), ncol(Y.obs))
  }

  N <- nrow(Y.obs)
  iT <- n_classes
  K <- ncol(Y.obs) #sum(K_h)

  # Number of free rows per item in mPhi
  n_free <- ivItemcat - 1L

  starts <- c(
    1L,
    cumsum(ifelse(ivItemcat == 2L, 1L, ivItemcat))[-length(ivItemcat)] + 1L
  )

  free_idx <- unlist(mapply(
    \(s, K_h) if (K_h == 2L) s else (s + 1L):(s + K_h - 1L),
    starts,
    ivItemcat,
    SIMPLIFY = FALSE
  ))

  phi_free <- fit0$mPhi[free_idx, ]

  theta1 <- c(
    fit0$vPi[2:iT],
    phi_free
  )

  p.xy <- compute_posteriors(Y.obs, mDesign, theta1, ivItemcat, iT)
  assignment <- max.col(p.xy)

  make_w <- function(posteriors) {
    if (use.modal.assignment) {
      w <- matrix(0L, nrow = N, ncol = iT)
      w[cbind(seq_len(N), max.col(posteriors))] <- 1L
      w
    } else {
      posteriors
    }
  }

  w.is <- make_w(p.xy)

  compute_pwx <- function(t1) {
    post <- compute_posteriors(Y.obs, mDesign, t1, ivItemcat, iT)
    w_local <- make_w(post)
    p.wx_joint <- (t(w_local) %*% post) / N
    sweep(p.wx_joint, 2, colSums(p.wx_joint), "/")
  }

  theta2_from_theta1 <- function(th1) {
    # rho <- c(1 - sum(th1[1:(iT - 1)]), th1[1:(iT - 1)])
    # phi <- matrix(th1[iT:length(th1)], nrow = K, ncol = iT)
    p.wx_mat <- compute_pwx(th1)
    log_ref <- log(diag(p.wx_mat))
    gamma_mat <- sweep(log(p.wx_mat), 2, log_ref, "-")
    gamma_mat[row(gamma_mat) != col(gamma_mat)]
  }

  gamma_vec_to_pwx <- function(gamma_vec) {
    gamma_mat <- matrix(0, nrow = iT, ncol = iT)
    gamma_mat[row(gamma_mat) != col(gamma_mat)] <- gamma_vec
    exp_mat <- exp(gamma_mat)
    sweep(exp_mat, 2, colSums(exp_mat), "/")
  }

  theta2 <- theta2_from_theta1(theta1)
  p.wx_mat <- gamma_vec_to_pwx(theta2)

  if (!use.simple.cov) {
    compute_J_unc_analytical <- function(
      p_ik,
      Y_obs,
      mDes,
      th1,
      ivItemcat,
      T_classes
    ) {
      N <- nrow(p_ik)

      A <- t(p_ik) %*% p_ik
      A[A < 1e-12] <- 1e-12

      #Extract item probabilities to match the free parameter structure
      phi_mat <- matrix(
        th1[T_classes:length(th1)],
        nrow = sum(ivItemcat - 1L),
        ncol = T_classes
      )

      L_rho <- T_classes - 1L
      L_phi <- sum((ivItemcat - 1L) * T_classes)
      L <- L_rho + L_phi

      J <- matrix(0, nrow = T_classes * (T_classes - 1L), ncol = L)

      #Offsets for locating items and categories in the expanded matrices
      starts_Y <- c(1L, cumsum(ivItemcat)[-length(ivItemcat)] + 1L)
      item_offsets <- c(
        0L,
        cumsum((ivItemcat - 1L) * T_classes)[-length(ivItemcat)]
      )

      # Map (s, t) to the correct row in J (matching gamma_mat[row != col] column-major)
      st_idx <- 1L
      st_map <- matrix(0L, nrow = T_classes, ncol = T_classes)
      for (t in seq_len(T_classes)) {
        for (s in seq_len(T_classes)) {
          if (s != t) {
            st_map[s, t] <- st_idx
            st_idx <- st_idx + 1L
          }
        }
      }

      for (t in seq_len(T_classes)) {
        P_tt <- (p_ik[, t]^2) / A[t, t]

        for (s in seq_len(T_classes)) {
          if (s == t) {
            next
          }
          row_J <- st_map[s, t]
          P_st <- (p_ik[, s] * p_ik[, t]) / A[s, t]

          for (c_prime in seq_len(T_classes)) {
            I_s <- if (s == c_prime) 1.0 else 0.0
            I_t <- if (t == c_prime) 1.0 else 0.0

            #shared derivative component for class c_prime
            Q_stc <- P_st *
              (I_s + I_t - 2 * p_ik[, c_prime]) -
              2 * P_tt * (I_t - p_ik[, c_prime])

            sum_Q <- sum(Q_stc)

            #Derivative for class prevalence u^rho (only c' >= 2)
            if (c_prime >= 2L) {
              col_rho <- c_prime - 1L
              J[row_J, col_rho] <- sum_Q
            }

            #Derivative for item response u^phi
            for (h in seq_along(ivItemcat)) {
              K_h <- ivItemcat[h]
              n_free <- K_h - 1L
              Y_cols <- (starts_Y[h] + 1L):(starts_Y[h] + K_h - 1L)

              phi_row_start <- if (h == 1L) {
                1L
              } else {
                sum(ivItemcat[1:(h - 1L)] - 1L) + 1L
              }
              phi_vals <- phi_mat[
                phi_row_start:(phi_row_start + n_free - 1L),
                c_prime
              ]

              col_J_start <- L_rho + item_offsets[h] + (c_prime - 1L) * n_free

              for (k in seq_len(n_free)) {
                Y_col <- Y_cols[k]

                val <- sum(Q_stc * Y_obs[, Y_col]) -
                  phi_vals[k] * sum(Q_stc * mDes[, Y_col])
                J[row_J, col_J_start + k] <- val
              }
            }
          }
        }
      }
      J
    }
  }

  list(
    theta1 = theta1,
    theta2 = theta2,
    w.is = w.is,
    p.wx_mat = p.wx_mat,
    gamma_vec_to_pwx = gamma_vec_to_pwx,
    theta2_from_theta1 = theta2_from_theta1,
    p.xy = p.xy,
    compute_J_unc = if (!use.simple.cov) compute_J_unc_analytical else NULL
  )
}

#' Analytic Hessian of the ML distal outcome negative log-likelihood
#'
#' Computes the T x T observed-data Hessian matrix at the current parameter
#' vector beta = (mu_1, ..., mu_T) for Gaussian, Poisson, or Binomial families,
#' accounting for classification error via the p.wx_mat correction matrix.
#' Used by lca_step3.distal to invert for the naive SE estimate.
#' @noRd
ml_hessian_distal <- function(
  beta,
  p.zx,
  pi_mat,
  pwx,
  Zo_cc,
  w.is_cc,
  iT,
  family,
  sigma2 = 1
) {
  N <- length(Zo_cc)
  pzx <- exp(pmax(p.zx(beta), -500))
  assignment_errors <- w.is_cc %*% pwx # N x T: P(W=s_i|X=t)
  q_i <- rowSums(pi_mat * pzx * assignment_errors) # N x 1
  r_it <- pi_mat * pzx * assignment_errors / q_i # N x T

  if (family == "gaussian") {
    mu <- beta
    g_it <- outer(Zo_cc, mu, "-") / sigma2 # N x T
    h_t <- rep(-1 / sigma2, iT)
  } else if (family == "poisson") {
    mu <- exp(beta)
    g_it <- outer(Zo_cc, rep(1, iT)) -
      outer(rep(1, N), mu) # N x T: z_i - mu_t
    h_t <- -mu
  } else if (family == "binomial") {
    mu <- 1 / (1 + exp(-beta))
    g_it <- outer(Zo_cc, rep(1, iT)) -
      outer(rep(1, N), mu) # N x T: z_i - mu_t
    h_t <- -mu * (1 - mu)
  }
  diag_term <- colSums(r_it) * h_t + colSums(r_it * g_it^2 * (1 - r_it))

  rg <- r_it * g_it # N x T
  off_diag <- -crossprod(rg) # T x T

  H_pos <- off_diag
  diag(H_pos) <- diag_term

  -H_pos # Hessian of neg.ll
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
    D <- qr.solve(pwx)
    w.it <- w.is_cc %*% D # N x T

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
    score_nt_ml <- function(mu, pzx, sigma2 = 1) {
      assignment_errors <- w.is_cc %*% pwx # N x T: P(W=s_i|X=t)
      q_i <- rowSums(pi_s * pzx * assignment_errors) # N x 1

      if (family == "gaussian") {
        resid <- outer(Zo_cc, mu, "-")
        pi_s * pzx * assignment_errors * resid / (q_i * sigma2)
      } else if (family == "poisson") {
        mu_val <- exp(mu)
        score <- outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val)
        pi_s * pzx * assignment_errors * score / q_i
      } else if (family == "binomial") {
        mu_val <- 1 / (1 + exp(-mu))
        score <- outer(Zo_cc, rep(1, iT)) - outer(rep(1, N), mu_val)
        pi_s * pzx * assignment_errors * score / q_i
      }
    }

    beta <- beta_init
    Z_long <- rep(Zo_cc, iT)
    X_long <- factor(rep(seq_len(iT), each = N))

    for (iter in seq_len(em.maxIter)) {
      pzx <- exp(pmax(p.zx(beta), -500))
      joint <- pi_s * pzx * (w.is_cc %*% pwx)
      w_tilde <- joint / rowSums(joint)
      w_long <- as.vector(w_tilde)

      fit <- glm(Z_long ~ X_long - 1, family = family, weights = w_long)
      beta_new <- coef(fit)

      if (abs(-neg.ll(beta_new) - (-neg.ll(beta))) < covariate.tol) {
        beta <- beta_new
        break
      }
      beta <- beta_new
      if (iter == em.maxIter) {
        warning("ML distal EM reached maximum iterations.")
      }
    }

    # estimate sigma2 after EM convergence (just for diagnostics)
    sigma2 <- if (family == "gaussian") {
      pzx <- exp(pmax(p.zx(beta), -500))
      joint <- pi_s * pzx * (w.is_cc %*% pwx)
      w_tilde <- joint / rowSums(joint)
      resid <- outer(Zo_cc, beta, "-")
      sum(w_tilde * resid^2) / sum(w_tilde)
    } else {
      NULL
    }

    three_step.score <- function(params) {
      mu <- params[1:iT]
      pzx <- exp(pmax(p.zx(params), -500))
      score_nt_ml(mu, pzx, sigma2 = 1)
    }

    H.3.inv <- tryCatch(
      qr.solve(ml_hessian_distal(
        beta,
        p.zx,
        pi_s,
        pwx,
        Zo_cc,
        w.is_cc,
        iT,
        family
      )),
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

  return(list(
    res = res,
    H.3.inv = H.3.inv,
    three_step.score = three_step.score
  ))
}

#' Step 3 (covariate): estimate multinomial logit gamma with either BCH or ML EM
#'
#' Optimizes the Q x (T-1) coefficient matrix gamma for P(X=t|Z_i) via
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
    D <- qr.solve(pwx)
    w.it <- w.is_cc %*% D # N x T
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
          "Try use.bch = FALSE to use the ML estimator instead. Or you can try increasing em.maxIter...",
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
    error = function(e) {
      if (!use.bch) {
        warning(
          "Exact Hessian inversion failed. Falling back to Hessian estimated from the cross-product of the case-wise log-likelihood gradients. This assumes a correct third-stage model specification."
        )
        matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
      } else {
        warning(
          "Hessian inversion failed after BCH. Try again with use.bch = FALSE."
        )
        matrix(NA_real_, Q * (iT - 1), Q * (iT - 1))
      }
    }
  )
  return(list(res = res, H.3.inv = H.3.inv))
}


# -- Variance estimation (Bakk et al., 2014) ----------------------------------

#' Individual-level BHHH variance matrix for binary and polytomous LCA
#'
#' Computes the outer-product (BHHH) information matrix and variance-covariance
#' matrix for LCA measurement model parameters in the unconstrained
#' (logit/log-ratio) space, matching \pkg{multilevLCA}'s \code{$Varmat}.
#'
#' The score in unconstrained space is
#' \eqn{s_{it} = u_{it}(y_i - d_i \circ p_{it})},
#' where \eqn{d_i} is the missing-data design indicator matrix.
#'
#' Assumes \code{fit0$mPhi} follows the \pkg{multilevLCA} storage convention:
#' \itemize{
#'   \item Dichotomous item h (\code{ivItemcat[h] == 2}): 1 row =
#'     \eqn{P(Y=1|C)}; the base level \eqn{P(Y=0|C)} is excluded.
#'   \item Polytomous item h (\code{ivItemcat[h] > 2}): \code{K_h} rows =
#'     \eqn{P(Y=0|C), \ldots, P(Y=K_h-1|C)}; the base level is included.
#' }
#' \code{expand_Y} produces one-hot columns in the same order so that
#' \code{expand_Phi(fit0$mPhi, ivItemcat)} aligns column-wise with
#' \code{expand_Y(mY, ivItemcat)}.  Free (estimable) parameters per item are
#' the single \eqn{P(Y=1|C)} row for dichotomous items, and rows 2 through
#' \eqn{K_h} for polytomous items (row 1, \eqn{P(Y=0|C)}, is the reference).
#' Boundary parameters (within \code{boundary.tol} of 0 or 1) are treated as
#' fixed: their score columns are zeroed and they do not contribute to the
#' information matrix.
#'
#' @param Y.exp       N x sum(K_h) expanded one-hot indicator matrix.
#' @param mDesign.exp Expanded design matrix (same dimensions as \code{Y.exp}),
#'   or \code{NULL} for complete data.
#' @param fit0        Step-1 fit object with \code{$vPi} and \code{$mPhi}.
#' @param ivItemcat   Integer vector of category counts per item.
#' @param boundary.tol Scalar tolerance for boundary detection. Default
#'   \code{1e-2}.
#' @param use.freq    Logical. Collapse duplicate score rows before computing
#'   the cross-product, weighting by frequency. Default \code{TRUE}.
#' @param u_post      Optional N x T matrix of posterior class probabilities.
#'   When supplied (e.g. extracted from \code{fit0$mU} via
#'   \code{extract_Y_from_mU}), \code{compute_posteriors} is skipped.
#'   Default \code{NULL}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`Infomat`}{Square BHHH information matrix of dimension p x p,
#'       where p = (iT-1) + sum(ivItemcat - 1) * iT is the total number of free
#'       parameters. Boundary parameters have zero rows and columns.}
#'     \item{`Varmat`}{Inverse of \code{Infomat} divided by N, giving the
#'       asymptotic variance-covariance matrix on the same scale as
#'       \pkg{multilevLCA}'s \code{$Varmat}. Boundary parameters have zero
#'       rows and columns.}
#'     \item{`SEs`}{Numeric vector of length p. Square root of the diagonal of
#'       \code{Varmat}; zero for boundary parameters.}
#'     \item{`mScore`}{N x p matrix of individual score contributions in the
#'       unconstrained parameterization, used for sandwich variance propagation
#'       in \code{lca_vcov} and \code{lca_vcov_distal}.}
#'   }
#' @keywords internal
lca_indiv_varmat <- function(
  Y.exp,
  mDesign.exp,
  fit0,
  ivItemcat,
  boundary.tol = 1e-2,
  use.freq = TRUE,
  u_post = NULL
) {
  pi_ <- fit0$vPi
  phi <- fit0$mPhi
  iT <- length(pi_)
  N <- nrow(Y.exp)

  if (is.null(mDesign.exp)) {
    mDesign.exp <- matrix(1L, N, ncol(Y.exp))
  }

  #Protect against parameter estimates on the boundary of the support (zero out their score contributions)
  pi_bdry <- pi_ <= boundary.tol | pi_ >= (1 - boundary.tol)
  phi_bdry <- phi <= boundary.tol | phi >= (1 - boundary.tol)

  #Clamp boundary parameters before computing posteriors
  pi_[pi_bdry] <- pmax(pmin(pi_[pi_bdry], 1 - 1e-6), 1e-6)
  phi[phi_bdry] <- pmax(pmin(phi[phi_bdry], 1 - 1e-6), 1e-6)

  # ---- Build theta1 and compute posteriors -----------------------------------
  starts <- c(
    1L,
    cumsum(ifelse(ivItemcat == 2L, 1L, ivItemcat))[-length(ivItemcat)] + 1L
  )
  free_idx <- unlist(mapply(
    \(s, K_h) if (K_h == 2L) s else (s + 1L):(s + K_h - 1L),
    starts,
    ivItemcat,
    SIMPLIFY = FALSE
  ))
  phi_free <- phi[free_idx, , drop = FALSE]

  # Use pre-computed posteriors when available (e.g. extracted from fit0$mU),
  # otherwise compute them from theta1.
  if (is.null(u_post)) {
    theta1 <- c(pi_[-1L], phi_free)
    u_post <- compute_posteriors(Y.exp, mDesign.exp, theta1, ivItemcat, iT)
  }

  # ---- Expand phi for residual computation -----------------------------------
  phi_exp <- expand_Phi(phi, ivItemcat) # K_total x T

  # ---- Build free_cols (columns of phi_exp for free parameters) -------------
  # free_cols: indices into phi_exp rows for the free categories.
  # phi_exp has K_total rows; free categories are the non-reference rows
  # within each item block. The mapping from free_idx (mPhi rows) to
  # phi_exp rows differs between dichotomous and polytomous items:
  #
  #   Dichotomous (K=2): phi_exp block = [P(Y=0), P(Y=1)], 2 rows.
  #     free_idx points to the single mPhi row = P(Y=1) = phi_exp row 2
  #     within the block (col_start + 1).
  #
  #   Polytomous (K>2): phi_exp block = [P(Y=0)..P(Y=K-1)], K rows.
  #     free_idx points to mPhi rows 2..K within the block (P(Y=1)..P(Y=K-1))
  #     = phi_exp rows col_start+1 .. col_start+K-1.
  #
  # In both cases the free phi_exp rows are exactly col_start+1..col_start+K-1
  # (dropping col_start = reference P(Y=0) for binary and polytomous alike).

  free_cols <- integer(0L)
  col_start <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    free_cols <- c(free_cols, (col_start + 1L):(col_start + K_h - 1L))
    col_start <- col_start + K_h
  }
  n_free_phi <- length(free_cols) # = sum(ivItemcat - 1) = nrow(phi_free) for all items

  # phi_bdry_free: n_free_phi x T, boundary flags for free phi parameters.
  # free_idx already selects the free mPhi rows in item order, so:
  phi_bdry_free <- phi_bdry[free_idx, , drop = FALSE]

  # ---- Pi scores (T-1 columns, one per non-reference class t=2..T) -----------
  s_u_pi <- sweep(u_post[, -1L, drop = FALSE], 2L, pi_[-1L], "-")
  pi_free_bdry <- pi_bdry[-1L] # drop reference class t=1
  if (any(pi_free_bdry)) {
    s_u_pi[, pi_free_bdry] <- 0
  }

  # ---- Phi scores (n_free_phi * T columns, class-major) ----------------------
  Y_free <- Y.exp[, free_cols, drop = FALSE]
  D_free <- mDesign.exp[, free_cols, drop = FALSE]

  s_u_phi <- matrix(0, N, n_free_phi * iT)

  for (t in seq_len(iT)) {
    idx <- ((t - 1L) * n_free_phi + 1L):(t * n_free_phi)
    resid <- Y_free -
      D_free *
        matrix(
          phi_exp[free_cols, t],
          nrow = N,
          ncol = n_free_phi,
          byrow = TRUE
        )
    s_col <- u_post[, t] * resid

    #Zero boundary free parameters for this class
    bdry_t <- phi_bdry_free[, t]
    if (any(bdry_t)) {
      s_col[, bdry_t] <- 0
    }

    s_u_phi[, idx] <- s_col
  }

  S <- cbind(s_u_pi, s_u_phi) # N x p

  # ---- Identify active (non-boundary) columns --------------------------------
  # Boundary parameters have all-zero score columns to avoid rank deficiency, then restore zero rows/cols after
  active <- which(colSums(S != 0) > 0L)
  p_full <- ncol(S)

  # ---- BHHH information matrix (on active columns only) ----------------------
  S_active <- S[, active, drop = FALSE]

  if (use.freq) {
    S_char <- apply(S_active, 1L, paste, collapse = "\r")
    uniq <- !duplicated(S_char)
    freq <- tabulate(match(S_char, S_char[uniq]))
    Infomat_active <- crossprod(S_active[uniq, , drop = FALSE] * sqrt(freq)) / N
  } else {
    Infomat_active <- crossprod(S_active) / N
  }

  Varmat_active <- tryCatch(
    qr.solve(Infomat_active) / N,
    error = function(e) {
      warning(
        "lca_indiv_varmat: Infomat is singular even after removing boundary ",
        "parameters; returning NA matrix. Check for near-empty classes.",
        call. = FALSE
      )
      matrix(NA_real_, length(active), length(active))
    }
  )

  # ---- Restore full-size Infomat and Varmat ----------------------------------
  Infomat <- matrix(0, p_full, p_full)
  Infomat[active, active] <- Infomat_active

  Varmat <- matrix(0, p_full, p_full)
  Varmat[active, active] <- Varmat_active

  list(
    Infomat = Infomat,
    Varmat = Varmat,
    SEs = sqrt(diag(Varmat)),
    mScore = S
  )
}

#' Covariate model variance-covariance with measurement-uncertainty correction
#'
#' Assembles the sandwich variance matrix for the Step-3 gamma estimates,
#' optionally propagating Step-1 measurement uncertainty through the analytic
#' `C_mat = d/d(theta2) [sum_i score_3_i]` and Jacobian `J.2 = d(theta2)/d(u)`.
#' When use.simple.cov = TRUE, returns the plain robust sandwich H^{-1} S'S H^{-1}.
#' @noRd
lca_vcov <- function(
  coefs,
  three_step.score,
  H.3.inv,
  Sigma.1,
  theta2,
  J.2,
  p.wx_mat,
  w.is,
  Z_mat,
  n_classes,
  p.xz,
  s2,
  use.simple.cov
) {
  J.3 <- three_step.score(c(coefs))

  Sigma.3.robust <- H.3.inv %*% crossprod(J.3) %*% H.3.inv
  if (!use.simple.cov) {
    # -- Analytic C_mat = d/d theta2 [colSums(score_3)] (checked with sympy) ----------------------------

    T_ <- n_classes
    pwx <- p.wx_mat
    pi_ <- p.xz(matrix(coefs, ncol = T_ - 1L)) # N x T
    N <- nrow(Z_mat)
    Q_ <- ncol(Z_mat)

    q <- pi_ %*% t(pwx) # N x T: q[i,s]
    V <- s2$w.is / q # N x T: w[i,s]/q[i,s]
    U <- s2$w.is / q^2 # N x T: w[i,s]/q[i,s]^2

    # VP[i, t0] = sum_s V[i,s] * pwx[s, t0]
    VP <- V %*% pwx

    n_theta2 <- T_ * (T_ - 1L)
    n_coef <- Q_ * (T_ - 1L)
    C_mat <- matrix(0, nrow = n_coef, ncol = n_theta2)

    # Map (s0, t0) mapping to the correct column index in C_mat
    theta2_idx <- matrix(0, T_, T_)
    theta2_idx[row(theta2_idx) != col(theta2_idx)] <- seq_len(n_theta2)

    for (t0 in seq_len(T_)) {
      for (s0 in seq_len(T_)) {
        if (s0 == t0) {
          next
        }
        idx_theta2 <- theta2_idx[s0, t0]

        # Store derivatives for all classes k for this specific (s0, t0)
        d_score_mat <- matrix(0, nrow = N, ncol = T_ - 1L)

        for (k in seq_len(T_ - 1L)) {
          # UP[i] = sum_s U[i,s] * pwx[s, k+1] * pwx[s, t0]
          UP_k_t0 <- rowSums(
            U *
              matrix(
                pwx[, k + 1] * pwx[, t0],
                nrow = N,
                ncol = T_,
                byrow = TRUE
              )
          )

          term1 <- 0
          if (k + 1L == t0) {
            term1 <- pi_[, k + 1] * (V[, s0] - VP[, t0])
          }

          term2 <- pi_[, k + 1] *
            pi_[, t0] *
            (U[, s0] * pwx[s0, k + 1] - UP_k_t0)

          d_score_mat[, k] <- pwx[s0, t0] * (term1 - term2)
        }

        C_mat[, idx_theta2] <- as.vector(crossprod(Z_mat, d_score_mat))
      }
    }

    step1.uncertainty <- C_mat %*% J.2 %*% Sigma.1 %*% t(J.2) %*% t(C_mat)
    Sigma.3 <- H.3.inv %*% (crossprod(J.3) + step1.uncertainty) %*% H.3.inv
  } else {
    Sigma.3 <- Sigma.3.robust
    step1.uncertainty <- NULL
  }

  param_names <- as.vector(outer(
    rownames(coefs),
    colnames(coefs),
    paste,
    sep = ":"
  ))
  rownames(Sigma.3) <- param_names
  colnames(Sigma.3) <- param_names

  Sigma.3
}

#' Distal outcome variance-covariance with full uncertainty propagation
#'
#' Assembles the T x T sandwich variance matrix for the distal outcome mu
#' estimates, propagating Step-1 measurement uncertainty
#' and, when both a covariate and distal model are fitted, Step-3 covariate
#' uncertainty. Skips steps 1/2 uncertainty corrections when
#' use.bch = TRUE or use.simple.cov = TRUE.
#' @noRd
lca_vcov_distal <- function(
  mu_hat,
  three_step.score,
  pi_adj,
  w.is,
  p.wx_mat,
  p.zx,
  family,
  H.3.inv,
  Sigma.1,
  s2,
  Sigma.3 = NULL,
  s3.par = NULL,
  p.xz.cov = NULL,
  Z_mat_cov = NULL,
  iT,
  use.simple.cov,
  use.bch
) {
  J.3 <- three_step.score(mu_hat)
  meat <- crossprod(J.3)

  if (use.bch || use.simple.cov) {
    result <- H.3.inv %*% meat %*% H.3.inv
    class_names <- paste0("mu_C", seq_len(iT))
    rownames(result) <- class_names
    colnames(result) <- class_names
    return(result)
  }

  N <- nrow(w.is)
  pwx <- p.wx_mat # T x T

  # -- Quantities at converged mu_hat --------------------------------------------
  pzx_mat <- exp(pmax(p.zx(mu_hat), -500)) # N x T: f(z_i|mu_t)
  ae <- w.is %*% pwx # N x T: sum_s w_{is}*pwx[s,t]
  q_i <- rowSums(pi_adj * pzx_mat * ae) # N: marginal density
  r_it <- pi_adj * pzx_mat * ae / q_i # N x T: posterior r_{it}

  # g_it = unit score of log f(z_i|mu_t) wrt mu_t
  # three_step.score[i,t] = r_{it} * g_{it}  =>  g_it = score / r_it
  g_it <- J.3 / pmax(r_it, 1e-300) # N x T

  # -- C1_mat: d/d theta2 [colSums(score_distal)]  (T x T*(T-1)) (checked with sympy) ----------------
  # theta2 = off-diagonal log-ratios of pwx (column-softmax parameterzation).
  # d ae_{i,t0}/d g_{s0,t0} = pwx[s0,t0] * (w_{i,s0} - ae_{i,t0})
  # c_i = dae / ae_{i,t0}
  # t==t0: dr_{i,t0} = r_{i,t0} * c_i * (1 - r_{i,t0})
  # t!=t0: dr_{i,t}  = -r_{i,t} * r_{i,t0} * c_i
  # C1[t, col] = sum_i dr_{it} * g_{it}

  n_theta2 <- iT * (iT - 1L)
  C1_mat <- matrix(0, iT, n_theta2)
  col_idx <- 0L

  for (t0 in seq_len(iT)) {
    for (s0 in seq_len(iT)) {
      if (s0 == t0) {
        next
      }
      col_idx <- col_idx + 1L
      v <- pwx[s0, t0]
      dae <- v * (w.is[, s0] - ae[, t0]) # N
      c_i <- dae / pmax(ae[, t0], 1e-300) # N
      for (t in seq_len(iT)) {
        dr <- if (t == t0) {
          r_it[, t0] * c_i * (1 - r_it[, t0])
        } else {
          -r_it[, t] * r_it[, t0] * c_i
        }
        C1_mat[t, col_idx] <- sum(dr * g_it[, t])
      }
    }
  }

  step1.uncertainty <- C1_mat %*%
    s2$J.2 %*%
    Sigma.1 %*%
    t(s2$J.2) %*%
    t(C1_mat)

  # -- C_mat: d/d gamma [colSums(score_distal)]  (T x Q*(T-1)) (checked with sympy) -----------------
  # gamma enters through pi_adj = p.xz(gamma); pwx_adj treated as fixed.
  # m_{it} = pzx_{it}*ae_i[t] / q_i   (proportional to r_{it}/pi_{it})
  # A_{i,l} = pi_{i,l+1} * (m_{i,l+1} - sum_t' m_{it'}*pi_{it'})
  # d r_{it}/d gamma_{q,l} = z_{iq} * [m_{it}*pi_{it}*(I(t==l+1)-pi_{i,l+1}) - r_{it}*A_{il}]
  # C_mat[t, l*Q+q] = sum_i score_mu[it] * z_{iq} * (above)

  step2.uncertainty <- matrix(0, iT, iT)

  if (!is.null(s3.par) && !is.null(p.xz.cov) && !is.null(Z_mat_cov)) {
    Q_cov <- ncol(Z_mat_cov)
    C_mat <- matrix(0, iT, (iT - 1L) * Q_cov)

    m_it <- pzx_mat * ae / q_i # N x T
    pi_cov <- p.xz.cov(matrix(s3.par, ncol = iT - 1L)) # N x T
    m_pi_sum <- rowSums(m_it * pi_cov) # N: sum_t m_{it}*pi_{it}

    for (l in seq_len(iT - 1L)) {
      A_il <- pi_cov[, l + 1L] * (m_it[, l + 1L] - m_pi_sum) # N
      for (t in seq_len(iT)) {
        delta_tl <- as.integer(t == l + 1L)
        inner_t <- m_it[, t] *
          pi_cov[, t] *
          (delta_tl - pi_cov[, l + 1L]) -
          r_it[, t] * A_il # N
        idx_l <- ((l - 1L) * Q_cov + 1L):(l * Q_cov)
        C_mat[t, idx_l] <- colSums(g_it[, t] * inner_t * Z_mat_cov)
      }
    }

    step2.uncertainty <- C_mat %*% Sigma.3 %*% t(C_mat)
  }

  result <- H.3.inv %*%
    (meat + step1.uncertainty + step2.uncertainty) %*%
    H.3.inv
  class_names <- paste0("mu_C", seq_len(iT))
  rownames(result) <- class_names
  colnames(result) <- class_names
  result
}


#' Three-step LCA estimation with covariates and/or distal outcomes
#'
#' Fits a three-step latent class model through the following steps:
#' \enumerate{
#'   \item \strong{Measurement model}: estimates latent class parameters
#'     (\eqn{\pi}, \eqn{\phi}) using \pkg{multilevLCA}
#'     (Lyrvall et al., 2025).
#'   \item \strong{Classification-error matrix}: computes posterior class
#'     probabilities and the T x T misclassification probability matrix
#'     \eqn{P(W = s \mid X = t)}, with standard errors corrected for
#'     classification-error propagation (Bakk, Oberski & Vermunt, 2014).
#'   \item \strong{Structural model}: estimates covariate effects using
#'     two-step starting values (Bakk & Kuha, 2018) and/or distal outcome
#'     means following Bakk, Tekle & Vermunt (2013), with the ML correction (Vermunt, 2010) or BCH correction
#'     (Bolck, Croon & Hagenaars, 2004). See
#'     \code{vignette("tseLCA", package = "tseLCA")} for a worked example.
#' }
#'
#' @param data A data.frame containing all columns referenced by \code{Y.names},
#'   \code{Zp.names}, and \code{Zo.name}.
#' @param Y.names Character vector of indicator column names. Need to be coded as consecutive integers with base level starting at `0`.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names, or \code{NULL}
#'   for a measurement-only fit. Default \code{NULL}.
#' @param Zo.name Single character name of the distal outcome column, or
#'   \code{NULL}. Default \code{NULL}.
#' @param step1 Pre-fitted Step-1 object (output of [tseLCA::lca_step1()] or a
#'   prior \code{three_step()} call), or \code{NULL} to run Step 1 internally.
#'   Default \code{NULL}.
#' @param use.two.step Logical. Initialize Step-3 from two-step estimates.
#'   Default \code{TRUE}.
#' @param use.modal.assignment Logical. Use modal (hard) class assignments in
#'   Step 2 and 3. \code{FALSE} uses soft posterior weights. Default \code{TRUE}.
#' @param include.intercept Logical. Prepend an intercept column to the
#'   covariate design matrix. Default \code{TRUE}.
#' @param use.simple.cov Logical. Skip the Step-1 measurement-uncertainty
#'   correction and return only the robust sandwich variance. Faster but
#'   underestimates standard errors when class separation is low. Default
#'   \code{FALSE}.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}.
#'   Default \code{FALSE}.
#' @param boundary.tol Scalar. Parameters within this tolerance of 0 or 1 are
#'   treated as fixed when computing the Step-1 variance matrix for numerical stability. Default
#'   \code{1e-2}.
#' @param maxIter.measurement Integer. Maximum EM iterations for Step 1.
#'   Default \code{5000L}.
#' @param measurement.tol Scalar. Convergence tolerance for the Step-1 EM
#'   algorithm. Default \code{1e-8}.
#' @param covariate.tol Scalar. Convergence tolerance for the Step-3
#'   Newton-Raphson or EM algorithm. Default \code{1e-6}.
#' @param iter.measurement Integer. Number of random restarts triggered when
#'   the Step-1 entropy R\eqn{^2} falls below \code{R2.threshold}. Default
#'   \code{10L}.
#' @param R2.threshold Scalar. Entropy R\eqn{^2} threshold below which Step-1
#'   random restarts are triggered. Default \code{0.70}.
#' @param use.bch Logical. Use BCH-corrected weights instead of the ML
#'   estimator in Step 3. May error if BCH weights induce a non-positive semi-definite Hessian in the third step (common in cases of low separation). Default \code{FALSE}.
#' @param em.maxIter Integer. Maximum EM iterations for the Step-3 covariate
#'   or distal outcome model. Default \code{200L}.
#' @param get.twostep.vcov Logical. If \code{TRUE}, obtain \pkg{multilevLCA}'s
#'   bias-corrected variance-covariance matrix for the two-step gamma estimates
#'   and store it in \code{$two_step_vcov}. If the \code{fitZ} object passed
#'   via \code{step1} already contains a \code{Varmat_cor} (from a prior
#'   [fitZ_from_multiLCA()] or plain \code{multiLCA} call), it is attached
#'   automatically even when \code{get.twostep.vcov = FALSE}. Default
#'   \code{FALSE}.
#' @param rebase Character (e.g. \code{"C1"}, \code{"C2"}) or integer
#'   specifying which latent class to use as the reference category in the
#'   multinomial logit. The measurement model is permuted so this class becomes
#'   column 1 before any structural estimation. Default \code{"C1"}.
#' @param family Character. Distal outcome family: one of \code{"gaussian"}
#'   (class means), \code{"poisson"} (log-rates), or \code{"binomial"}
#'   (logits). Default \code{"gaussian"}.
#' @param correct.spec Logical. Use the model-robust outer-product Hessian for
#'   Step-3 standard errors rather than the observed-data Hessian. Not appropriate
#'   when the Step-3 model may be misspecified. Default \code{FALSE}.
#' @param verbose Logical. Print convergence messages. Default \code{FALSE}.
#'
#' @return An S3 object of class \code{tseLCA}. The subclass depends on which
#'   models were estimated:
#'   \describe{
#'     \item{`tseLCA_measurement`}{Returned when neither \code{Zp.names} nor
#'       \code{Zo.name} is supplied. Contains the following elements:
#'       \describe{
#'         \item{`measurement_model`}{Step-1 output list from [tseLCA::lca_step1()].}
#'         \item{`llik`}{Final Step-1 log-likelihood.}
#'         \item{`AIC`, `BIC`}{Information criteria from the measurement model.}
#'         \item{`R2entr`}{Entropy R\eqn{^2} of the measurement model.}
#'         \item{`n_classes`}{Number of latent classes.}
#'         \item{`posteriors`}{N x T matrix of soft posterior class probabilities.}
#'         \item{`classifications`}{Length-N integer vector of modal class assignments.}
#'       }
#'     }
#'     \item{`tseLCA_covariate`}{Returned when \code{Zp.names} is supplied and
#'       \code{Zo.name} is \code{NULL}. Contains all elements of
#'       \code{tseLCA_measurement} plus:
#'       \describe{
#'         \item{`three_step`}{Q x (T-1) matrix of Step-3 gamma coefficients.}
#'         \item{`three_step_vcov`}{Q(T-1) x Q(T-1) variance-covariance matrix
#'           for \code{three_step}, with measurement-uncertainty correction
#'           unless \code{use.simple.cov = TRUE}.}
#'         \item{`two_step`}{Q x (T-1) matrix of two-step starting values, or
#'           \code{NULL} if \code{use.two.step = FALSE}.}
#'         \item{`two_step_vcov`}{\pkg{multilevLCA} bias-corrected vcov for the
#'           two-step estimates, or \code{NULL}.}
#'         \item{`estimator`}{Character: \code{"ML"} or \code{"BCH"}.}
#'         \item{`entropy.R2`}{Covariate-adjusted entropy R\eqn{^2}.}
#'         \item{`llik`}{Profile log-likelihood
#'           \eqn{\sum_i \log \sum_t P(X=t|Z_{p,i};\hat{\gamma}) P(Y_i|X=t;\hat{\phi})},
#'           with Step-1 parameters \eqn{\hat{\phi}} held fixed. By construction
#'           smaller than the equivalent one-step MLE likelihood.}
#'       }
#'     }
#'     \item{`tseLCA_distal`}{Returned when \code{Zo.name} is supplied and
#'       \code{Zp.names} is \code{NULL}. Contains:
#'       \describe{
#'         \item{`three_step`}{Named length-T vector of Step-3 distal outcome
#'           parameters (means, log-rates, or logits depending on \code{family}).}
#'         \item{`three_step_vcov`}{T x T variance-covariance matrix for
#'           \code{three_step}, named \code{mu_C1} through \code{mu_CT}.}
#'         \item{`three_step.llik`}{Step-3 distal log-likelihood
#'           \eqn{\log P(Z_o|X=t)} at converged estimates.}
#'         \item{`llik`}{Profile log-likelihood
#'           \eqn{\sum_i \log \sum_t P(X=t|\hat{\pi}) P(Z_{o,i}|X=t;\hat{\mu}) P(Y_i|X=t;\hat{\phi})},
#'           with Step-1 parameters \eqn{\hat{\pi}, \hat{\phi}} held fixed.
#'           By construction smaller than the equivalent one-step MLE likelihood.}
#'         \item{`AIC`}{Akaike information criterion based on \code{llik}.}
#'         \item{`BIC`}{Bayesian information criterion based on \code{llik},
#'           using the number of distal-complete observations.}
#'         \item{`family`}{Character. The distal outcome family used.}
#'         \item{`estimator`}{Character: \code{"ML"} or \code{"BCH"}.}
#'         \item{`posteriors`}{N x T soft posterior matrix.}
#'         \item{`classifications`}{Length-N modal class assignment vector.}
#'       }
#'     }
#'     \item{`tseLCA_both`}{Returned when both \code{Zp.names} and
#'       \code{Zo.name} are supplied. Contains:
#'       \describe{
#'         \item{`covariate`}{A \code{tseLCA_covariate}-structured sub-list
#'           (see above), including \code{llik}, \code{AIC}, \code{BIC},
#'           \code{entropy.R2}.}
#'         \item{`distal`}{A \code{tseLCA_distal}-structured sub-list
#'           (see above), including \code{llik}, \code{AIC}, \code{BIC},
#'           \code{three_step.llik}.}
#'         \item{`family`, `n_classes`, `estimator`}{Shared top-level fields.}
#'         \item{`posteriors`, `classifications`}{Shared N x T posterior
#'           matrix and length-N modal class vector.}
#'       }
#'     }
#'   }
#'
#' @references
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#'   between latent class membership and external variables using bias-adjusted
#'   three-step approaches. \emph{Sociological Methodology}, 43(1), 272--311.
#'   \doi{10.1177/0081175012470644}
#'
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#'   classes and external variables. \emph{Psychometrika}, 83(4), 871--892.
#'   \doi{10.1007/s11336-017-9592-7}
#'
#' Bakk, Z., Pohle, M. J., & Kuha, J. (2025). Bias-adjusted three-step
#'   estimation of structural models for latent classes. \emph{Multivariate
#'   Behavioral Research}. \doi{10.1080/00273171.2025.2473935}
#'
#' @seealso \code{vignette("tseLCA", package = "tseLCA")} for a full worked
#'   example; [tseLCA::lca_step1()] for standalone Step-1 estimation;
#'   [fitZ_from_fit0()] and [fitZ_from_multiLCA()] for two-step covariate
#'   estimation.
#'
#' @examples
#' d <- generate_data(n = 200, separation = "high",
#'                    scenario = "covariate", seed = 1)
#'
#' # Measurement model only
#' fit_m <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' summary(fit_m)
#'
#' # ML three-step with simple SEs (fast)
#' fit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#'
#' # Full measurement-uncertainty correction (see vignette for interpretation)
#' fit_cor <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.simple.cov = FALSE,
#'                       use.modal.assignment = FALSE)
#' summary(fit_cor)
#'
#' # BCH estimator
#' fit_bch <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.bch = TRUE,
#'                       use.simple.cov = TRUE)
#' summary(fit_bch)
#'
#' # Change reference class
#' fit_c2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                      Zp.names = "Zp", use.simple.cov = TRUE,
#'                      rebase = "C2")
#' summary(fit_c2)
#'
#' # Gaussian distal outcome
#' d2 <- generate_data(200, "high", "distal", seed = 2)
#' fit_dis <- three_step(d2, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zo.name = "Zo", family = "gaussian",
#'                       use.simple.cov = TRUE)
#' summary(fit_dis)
#'
#' # Pass a pre-fitted measurement model to skip Step 1
#' fit_step1 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' fit2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                    Zp.names = "Zp", step1 = fit_step1,
#'                    use.simple.cov = TRUE)
#' summary(fit2)
#'
#' # Plot item-response profiles from the measurement model
#' plot(fit)
#'
#' @export
three_step <- function(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  Zo.name = NULL,
  step1 = NULL,
  use.two.step = TRUE,
  use.modal.assignment = TRUE,
  include.intercept = TRUE,
  use.simple.cov = FALSE,
  incomplete = FALSE,
  boundary.tol = 1e-2,
  maxIter.measurement = 5000,
  measurement.tol = 1e-8,
  covariate.tol = 1e-6,
  iter.measurement = 10L,
  R2.threshold = 0.70,
  use.bch = FALSE,
  em.maxIter = 200L,
  get.twostep.vcov = FALSE,
  rebase = "C1",
  family = "gaussian",
  correct.spec = FALSE,
  verbose = FALSE
) {
  # -- Step 1: measurement model -----------------------------------------------
  ref_idx <- parse_rebase(rebase, n_classes)

  if (!is.null(step1)) {
    # Normalize: accept raw lca_step1() list or any tseLCA object
    s1 <- if (inherits(step1, "tseLCA")) step1$measurement_model else step1
    # Apply rebase permutation so the desired reference class is column 1.
    s1$fit0 <- permute_fit0_classes(s1$fit0, ref_idx)
    if (!is.null(s1$fitZ)) {
      s1$fitZ <- normalize_fitZ_names(
        s1$fitZ,
        n_classes = n_classes
      )
      s1$fitZ <- permute_fitZ_classes(s1$fitZ, ref_idx)
    }
  } else {
    s1 <- lca_step1(
      data,
      Y.names,
      n_classes,
      Zp.names,
      maxIter.measurement,
      measurement.tol,
      covariate.tol,
      iter.measurement,
      R2.threshold,
      get.twostep.vcov,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  }
  fit0 <- s1$fit0

  # Compute two-step coefficients with EM algorithm holding step1 fixed.
  # This runs when use.two.step = TRUE and no fitZ was already computed
  # (e.g. when step1 was passed in from a measurement-only fit).
  if (use.two.step && is.null(s1$fitZ) && !is.null(Zp.names)) {
    s1$fitZ <- fitZ_from_fit0(
      fit0 = s1$fit0,
      data = data,
      Y.names = Y.names,
      Zp.names = Zp.names,
      tol = covariate.tol,
      maxIter = em.maxIter,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  }
  fitZ <- s1$fitZ

  # -- Data preparation --------------------------------------------------------
  cd <- clean_data(
    data = data,
    Y.names = Y.names,
    Zp.names = Zp.names,
    Zo.name = Zo.name,
    incomplete = incomplete,
    include.intercept = include.intercept,
    verbose = verbose
  )
  Y.obs <- cd$Y.obs
  mDesign <- cd$mDesign
  ivItemcat <- cd$ivItemcat
  keep_Y <- cd$keep_Y
  Z_mat <- cd$Z_mat
  keep_step3_Z_in_Y <- cd$keep_step3_Z_in_Y
  Zo_mat <- cd$Zo_mat
  keep_step3_Zo_in_Y <- cd$keep_step3_Zo_in_Y
  keep_step3_Zo <- cd$keep_step3_Zo

  s1$Y.names <- Y.names
  s1$ivItemcat <- ivItemcat
  s1$ref_idx <- ref_idx

  # -- Early return if no covariates -------------------------------------------
  if (is.null(Zp.names) && is.null(Zo.name)) {
    # Extract posteriors and modal classifications from fit0$mU
    mU_posts <- if (!is.null(s1$fit0$mU)) {
      K_total <- sum(ivItemcat)
      n_mU_Y <- sum(ifelse(ivItemcat == 2L, 1L, ivItemcat))
      posts <- s1$fit0$mU[, (n_mU_Y + 1L):(n_mU_Y + n_classes), drop = FALSE]
      mode(posts) <- "double"
      posts
    } else {
      NULL
    }
    return(structure(
      list(
        measurement_model = s1,
        llik = s1$fit0$LLKSeries[nrow(s1$fit0$LLKSeries)],
        AIC = s1$fit0$AIC,
        BIC = s1$fit0$BIC,
        R2entr = s1$fit0$R2entr,
        n_classes = n_classes,
        posteriors = mU_posts,
        classifications = if (!is.null(mU_posts)) max.col(mU_posts) else NULL
      ),
      class = c("tseLCA_measurement", "tseLCA")
    ))
  }

  # -- Step 2 ------------------------------------------------------------------
  s2 <- lca_step2(
    Y.obs,
    fit0,
    n_classes,
    use.modal.assignment,
    boundary.tol,
    use.simple.cov || use.bch,
    ivItemcat = ivItemcat,
    mDesign = mDesign
  )

  Q <- if (!is.null(Z_mat)) ncol(Z_mat) else 0L
  iT <- n_classes

  # Subset s2 outputs to the rows used in each Step 3 model.
  # s2 was estimated on all keep_Y rows; Step 3 models only use complete-Z rows.
  s2_for_cov <- if (!is.null(Z_mat)) {
    J.2_cov <- if (!is.null(s2$compute_J_unc)) {
      Y_cov <- Y.obs[keep_step3_Z_in_Y, , drop = FALSE]
      mDes_cov <- if (!is.null(mDesign)) {
        mDesign[keep_step3_Z_in_Y, , drop = FALSE]
      } else {
        matrix(1L, nrow(Y_cov), ncol(Y_cov))
      }
      p_xy_cov <- s2$p.xy[keep_step3_Z_in_Y, , drop = FALSE]
      s2$compute_J_unc(
        p_xy_cov,
        Y_cov,
        mDes_cov,
        s2$theta1,
        ivItemcat,
        iT
      )
    } else {
      NULL
    }

    list(
      theta1 = s2$theta1,
      theta2 = s2$theta2,
      p.wx_mat = s2$p.wx_mat,
      gamma_vec_to_pwx = s2$gamma_vec_to_pwx,
      theta2_from_theta1 = s2$theta2_from_theta1,
      J.2 = J.2_cov,
      w.is = s2$w.is[keep_step3_Z_in_Y, , drop = FALSE],
      post = if (!is.null(s2$post)) {
        s2$post[keep_step3_Z_in_Y, , drop = FALSE]
      } else {
        NULL
      }
    )
  } else {
    NULL
  }

  s2_for_dis <- if (!is.null(Zo_mat)) {
    J.2_dis <- if (!is.null(s2$compute_J_unc)) {
      Y_dis <- Y.obs[keep_step3_Zo_in_Y, , drop = FALSE]
      mDes_dis <- if (!is.null(mDesign)) {
        mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
      } else {
        matrix(1L, nrow(Y_dis), ncol(Y_dis))
      }
      p_xy_dis <- s2$p.xy[keep_step3_Zo_in_Y, , drop = FALSE]
      s2$compute_J_unc(
        p_xy_dis,
        Y_dis,
        mDes_dis,
        s2$theta1,
        ivItemcat,
        iT
      )
    } else {
      NULL
    }

    list(
      theta1 = s2$theta1,
      theta2 = s2$theta2,
      p.wx_mat = s2$p.wx_mat,
      gamma_vec_to_pwx = s2$gamma_vec_to_pwx,
      theta2_from_theta1 = s2$theta2_from_theta1,
      J.2 = J.2_dis,
      w.is = s2$w.is[keep_step3_Zo_in_Y, , drop = FALSE],
      post = if (!is.null(s2$post)) {
        s2$post[keep_step3_Zo_in_Y, , drop = FALSE]
      } else {
        NULL
      }
    )
  } else {
    NULL
  }

  # Extract Step-1 sample inputs for Sigma.1. Uses fit0$mU when available.
  step1_Y <- if (!is.null(fit0$mU)) {
    raw <- extract_Y_from_mU(fit0, ivItemcat)
    if (ref_idx != 1L) {
      T_ <- n_classes
      ord <- c(ref_idx, seq_len(T_)[-ref_idx])
      raw$u_post <- raw$u_post[, ord, drop = FALSE]
    }
    raw
  } else {
    list(Y.exp = Y.obs, mDesign = mDesign, ivItemcat = ivItemcat, u_post = NULL)
  }

  #For covariate estimation
  if (!is.null(Z_mat)) {
    p.xz <- function(params) {
      eta_full <- cbind(0, Z_mat %*% params)
      row_max <- apply(eta_full, 1, max)
      exp_eta <- exp(eta_full - row_max)
      exp_eta / rowSums(exp_eta)
    }

    if (use.bch) {
      D <- qr.solve(s2_for_cov$p.wx_mat)
      w.it <- s2_for_cov$w.is %*% D

      .ll_bch <- function(params, pwx = NULL) {
        beta.cur <- matrix(params, ncol = iT - 1)
        rowSums(w.it * log(p.xz(beta.cur)))
      }

      neg.ll <- function(params) {
        beta.cur <- matrix(params, ncol = iT - 1)
        -sum(w.it * log(pmax(p.xz(beta.cur), 1e-6)))
      }

      .grad_bch <- function(params) {
        beta.cur <- matrix(params, ncol = iT - 1)
        pi_ <- p.xz(beta.cur)
        resid <- w.it[, -1L, drop = FALSE] -
          pi_[, -1L, drop = FALSE] * rowSums(w.it)
        -as.vector(t(Z_mat) %*% resid)
      }

      .score_bch <- function(params) {
        beta.cur <- matrix(params, ncol = iT - 1)
        pi_ <- p.xz(beta.cur)
        resid <- w.it[, -1L, drop = FALSE] -
          pi_[, -1L, drop = FALSE] * rowSums(w.it)
        resid[, rep(seq_len(iT - 1L), each = Q)] *
          Z_mat[, rep(seq_len(Q), iT - 1L)]
      }

      three_step.ll <- .ll_bch
      three_step.grad <- .grad_bch
      three_step.score <- .score_bch
    } else {
      .ll_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
        probs <- p.xz(matrix(params, ncol = iT - 1))
        rowSums(s2_for_cov$w.is * log(probs %*% t(pwx)))
      }

      .grad_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
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

      .score_ml <- function(params, pwx = s2_for_cov$p.wx_mat) {
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

      three_step.ll <- .ll_ml
      three_step.grad <- .grad_ml
      three_step.score <- .score_ml
      neg.ll <- function(params) -sum(three_step.ll(params))
    }

    gamma_init <- if (!is.null(fitZ$mGamma) && use.two.step) {
      c(fitZ$mGamma)
    } else {
      rep(0, Q * (iT - 1))
    }

    # -- Step 3 ------------------------------------------------------------------
    s3 <- lca_step3(
      neg.ll,
      gamma_init,
      Q,
      iT,
      covariate.tol,
      gradient = three_step.grad,
      use.bch = use.bch,
      Z_mat_cc = Z_mat,
      w.is_cc = s2_for_cov$w.is,
      p.xz = p.xz,
      pwx = s2_for_cov$p.wx_mat,
      em.maxIter = em.maxIter,
      verbose = verbose,
      correct.spec = correct.spec
    )
    if (
      (correct.spec && !use.bch) ||
        is.null(s3$H.3.inv) ||
        any(is.na(s3$H.3.inv))
    ) {
      s3$H.3.inv <- qr.solve(crossprod(three_step.score(s3$res$par)))
    }

    coefs <- matrix(s3$res$par, ncol = iT - 1)
    ref_idx <- parse_rebase(rebase, iT)
    non_ref_classes <- seq_len(iT)[-ref_idx]
    colnames(coefs) <- paste0("C", non_ref_classes)
    rownames(coefs) <- c("Intercept", Zp.names)

    # -- Variance -----------------------------------------------------------------
    Sigma.3 <- lca_vcov(
      coefs = coefs,
      three_step.score = three_step.score,
      H.3.inv = s3$H.3.inv,
      Sigma.1 = if (use.simple.cov || use.bch) {
        NULL
      } else {
        lca_indiv_varmat(
          step1_Y$Y.exp,
          step1_Y$mDesign,
          fit0,
          step1_Y$ivItemcat,
          boundary.tol = boundary.tol,
          u_post = step1_Y$u_post
        )$Varmat
      },
      J.2 = s2_for_cov$J.2,
      p.wx_mat = s2_for_cov$p.wx_mat,
      w.is = s2_for_cov$w.is,
      Z_mat = Z_mat,
      n_classes = n_classes,
      p.xz = p.xz,
      s2 = s2_for_cov,
      use.simple.cov = use.simple.cov || use.bch
    )
    Sigma.3.covariate <- Sigma.3

    # -- Model fit ----------------------------------------------------------------
    Y_cc <- Y.obs[keep_step3_Z_in_Y, , drop = FALSE]
    mDes_cc <- if (!is.null(mDesign)) {
      mDesign[keep_step3_Z_in_Y, , drop = FALSE]
    } else {
      NULL
    }
    total.llik <- joint_log_lik(
      Y_cc,
      Z_mat,
      expand_Phi(fit0$mPhi, ivItemcat),
      coefs,
      mDes_cc
    )
    total.k <- (iT * ncol(Y.obs)) + (Q * (iT - 1))
    total.AIC <- -2 * total.llik + 2 * total.k
    total.BIC <- -2 * total.llik + total.k * log(nrow(Y_cc))

    # -- Optional two-step vcov from multiLCA ------------------------------------
    # get.twostep.vcov = TRUE requests multilevLCA's bias-corrected SEs.
    # We skip re-estimation if a Varmat_cor is already attached to fitZ --
    # this handles three cases:
    #   (a) fitZ_from_multiLCA output: Varmat_cor lives at fitZ$raw_fit$Varmat_cor
    #   (b) Plain multiLCA output passed directly: fitZ$Varmat_cor
    #   (c) fitZ_from_fit0 output: no Varmat_cor anywhere -> must re-estimate
    #
    # If a Varmat_cor is already present on fitZ (regardless of get.twostep.vcov),
    # we always attach it to the output -- the user shouldn't lose it just because
    # they didn't set get.twostep.vcov = TRUE.

    .extract_varmat <- function(fZ) {
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

    .name_varmat <- function(V, fZ) {
      if (is.null(V) || is.null(fZ$mGamma)) {
        return(V)
      }
      param_names <- as.vector(outer(
        rownames(fZ$mGamma),
        colnames(fZ$mGamma),
        paste,
        sep = ":"
      ))
      rownames(V) <- param_names
      colnames(V) <- param_names
      V
    }

    existing_varmat <- .extract_varmat(fitZ)

    two_step_vcov <- if (!is.null(existing_varmat)) {
      # Already have it
      .name_varmat(existing_varmat, fitZ)
    } else if (get.twostep.vcov) {
      # No existing vcov
      fZ_ml <- fitZ_from_multiLCA(
        data = data,
        Y.names = Y.names,
        n_classes = n_classes,
        Zp.names = Zp.names,
        maxIter.measurement = maxIter.measurement,
        measurement.tol = measurement.tol,
        covariate.tol = covariate.tol,
        iter.measurement = iter.measurement,
        R2.threshold = R2.threshold,
        incomplete = incomplete,
        rebase = rebase,
        verbose = verbose
      )
      if (is.null(fitZ)) {
        fitZ <- fZ_ml
        s1$fitZ <- fZ_ml
      }
      raw_varmat <- .extract_varmat(fZ_ml)
      if (!is.null(raw_varmat)) {
        .name_varmat(raw_varmat, fZ_ml)
      } else {
        warning(
          "get.twostep.vcov: neither Varmat_cor nor SEs_cor_gamma found."
        )
        NULL
      }
    } else {
      NULL
    }

    # -- Covariate-adjusted entropy R^2 ------------------------------------------
    # Measures how much the items reduce classification uncertainty *beyond*
    # what the covariates already explain.
    #   error_prior = H(X|Z):   average entropy of P(X|Z_i) under fitted gamma
    #   error_post  = H(X|Y,Z): average entropy of P(X|Y_i,Z_i), recomputed
    #                            with the covariate-adjusted prior via
    #                            compute_pwx_adj (soft posterior assignment)
    #   R^2 = (H(X|Z) - H(X|Y,Z)) / H(X|Z)
    .h <- function(p) {
      p <- p[p > sqrt(.Machine$double.eps)]
      -sum(p * log(p))
    }
    pi_adj_cov <- p.xz(matrix(s3$res$par, ncol = iT - 1L)) # N x T

    error_prior <- mean(apply(pi_adj_cov, 1L, .h))

    adj_res <- compute_pwx_adj(
      Y.obs = Y_cc,
      fit0 = fit0,
      ivItemcat = ivItemcat,
      mDesign = mDes_cc,
      use.modal.assignment = FALSE,
      pi_adj = pi_adj_cov
    )
    error_post <- mean(apply(adj_res$post, 1L, .h))

    entropy.R2 <- if (error_prior > 1e-8) {
      (error_prior - error_post) / error_prior
    } else {
      1.0 # covariates already explain all class membership
    }

    s3.covariate <- structure(
      list(
        measurement_model = s1,
        two_step = if (!is.null(fitZ)) fitZ$mGamma else NULL,
        two_step_vcov = two_step_vcov,
        three_step = coefs,
        three_step_vcov = Sigma.3,
        three_step.llik = -s3$res$value,
        neg.ll = neg.ll,
        llik = total.llik,
        AIC = total.AIC,
        BIC = total.BIC,
        n_classes = iT,
        estimator = if (use.bch) "BCH" else "ML",
        entropy.R2 = entropy.R2,
        posteriors = s2$p.xy,
        classifications = max.col(s2$p.xy)
      ),
      class = c("tseLCA_covariate", "tseLCA")
    )
  }

  if (!is.null(Zo_mat)) {
    if (!(family %in% c("gaussian", "poisson", "binomial"))) {
      message(
        'Provided family is not one of "gaussian", "poisson", nor "binomial". Defaulting to family="gaussain".'
      )
    }

    if (!is.null(Zp.names)) {
      Z_mat_dis <- if (!is.null(Z_mat) && length(keep_step3_Zo) > 0L) {
        Z_full_raw <- if (include.intercept) {
          m <- cbind(1, as.matrix(data[, Zp.names, drop = FALSE]))
          colnames(m) <- c("Intercept", Zp.names)
          m
        } else {
          as.matrix(data[, Zp.names, drop = FALSE])
        }
        Z_full_raw[keep_step3_Zo, , drop = FALSE]
      } else {
        NULL
      }

      if (!is.null(Z_mat_dis)) {
        pi_adj_full <- cbind(1, p.xz(matrix(s3$res$par, ncol = iT - 1)))
        p.xz_dis <- function(params) {
          eta_full <- cbind(0, Z_mat_dis %*% params)
          row_max <- apply(eta_full, 1L, max)
          exp_eta <- exp(eta_full - row_max)
          exp_eta / rowSums(exp_eta)
        }
        pi_adj <- p.xz_dis(matrix(s3$res$par, ncol = iT - 1))
      } else {
        pi_adj <- matrix(
          fit0$vPi,
          nrow = length(keep_step3_Zo_in_Y),
          ncol = iT,
          byrow = TRUE
        )
      }

      res_adj <- compute_pwx_adj(
        Y.obs[keep_step3_Zo_in_Y, , drop = FALSE],
        fit0,
        ivItemcat,
        if (!is.null(mDesign)) {
          mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
        } else {
          NULL
        },
        use.modal.assignment,
        pi_adj = pi_adj
      )
    } else {
      pi_adj <- matrix(
        fit0$vPi,
        nrow = length(keep_step3_Zo_in_Y),
        ncol = iT,
        byrow = TRUE
      )
      res_adj <- list(
        w.is = s2_for_dis$w.is,
        p.wx_mat = s2_for_dis$p.wx_mat
      )
    }

    w.is_dis <- res_adj$w.is

    # Create p.zx : The function that maps latent indicators to oberved distal outcome Zo (need a choice of likelihood)

    if (family == "poisson") {
      p.zx <- function(params) {
        log_mu <- params[1:iT]
        mu <- exp(log_mu)
        z <- Zo_mat[, 1L]
        outer(z, log_mu, "*") - # N x T: z_i * log(mu_t)
          outer(rep(1, nrow(Zo_mat)), mu, "*") - # N x T: mu_t
          lgamma(z + 1L) # N x 1, recycled
      }
      starting.lm <- glm(
        Zo_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)),
        family = poisson()
      )
      beta_init <- coef(starting.lm)
    } else if (family == "binomial") {
      # params: logit(mu_t)
      p.zx <- function(params) {
        logit_mu <- params[1:iT]
        mu <- 1 / (1 + exp(-logit_mu))
        z <- Zo_mat[, 1L]
        outer(z, log(mu), "*") + # N x T
          outer(1 - z, log(1 - mu), "*")
      }
      starting.lm <- glm(
        Zo_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)),
        family = binomial()
      )
      beta_init <- coef(starting.lm) # already on logit scale
    } else {
      #(family == "gaussian")
      p.zx <- function(params) {
        mu <- params[1:iT]
        resid <- outer(Zo_mat[, 1L], mu, "-")
        -0.5 * resid^2 - 0.5 * log(2 * pi)
      }
      starting.lm <- lm(Zo_mat[, 1L] ~ -1 + as.factor(max.col(w.is_dis)))
      beta_init <- coef(starting.lm)
    }

    if (use.bch) {
      D <- qr.solve(res_adj$p.wx_mat)
      w.it <- res_adj$w.is %*% D

      neg.ll <- function(params) {
        -sum(w.it * p.zx(params))
      }
    } else {
      neg.ll <- function(params) {
        pzx <- exp(pmax(p.zx(params), -500))
        # classification error probabilities for each person's assignment
        assignment_errors <- res_adj$w.is %*% res_adj$p.wx_mat
        -sum(log(rowSums(pi_adj * pzx * assignment_errors)))
      }
    }

    s3.distal <- lca_step3.distal(
      neg.ll = neg.ll,
      em.maxIter = em.maxIter,
      pwx = res_adj$p.wx_mat,
      w.is_cc = res_adj$w.is,
      Zo_cc = Zo_mat[, 1L],
      use.bch = use.bch,
      covariate.tol = covariate.tol,
      iT = iT,
      beta_init = beta_init,
      family = family,
      p.zx = p.zx,
      vPi = fit0$vPi,
      pi_mat = pi_adj
    )

    #Variance-covariance
    Sigma.3.distal <- lca_vcov_distal(
      mu_hat = s3.distal$res$par,
      three_step.score = s3.distal$three_step.score,
      pi_adj = pi_adj,
      w.is = res_adj$w.is,
      p.wx_mat = res_adj$p.wx_mat,
      p.zx = p.zx,
      family = family,
      H.3.inv = s3.distal$H.3.inv,
      Sigma.1 = if (use.simple.cov || use.bch) {
        NULL
      } else {
        lca_indiv_varmat(
          step1_Y$Y.exp,
          step1_Y$mDesign,
          fit0,
          step1_Y$ivItemcat,
          boundary.tol = boundary.tol,
          u_post = step1_Y$u_post
        )$Varmat
      },
      s2 = s2_for_dis,
      Sigma.3 = if (!is.null(Zp.names)) Sigma.3 else NULL,
      s3.par = if (!is.null(Zp.names)) s3$res$par else NULL,
      p.xz.cov = if (!is.null(Zp.names) && exists("p.xz_dis")) {
        p.xz_dis
      } else if (!is.null(Zp.names)) {
        p.xz
      } else {
        NULL
      },
      Z_mat_cov = if (!is.null(Zp.names) && !is.null(Z_mat_dis)) {
        Z_mat_dis
      } else if (!is.null(Zp.names)) {
        Z_mat
      } else {
        NULL
      },
      iT = iT,
      use.simple.cov = use.simple.cov,
      use.bch = use.bch
    )

    distal_par <- s3.distal$res$par
    names(distal_par) <- paste0("mu_C", seq_len(iT))

    # -- Distal log-likelihood, AIC, BIC ----------------------------------------
    # Step-3 llik: log P(Zo|X=t) weighted by class assignments (from optim).
    # Total joint llik: sum_i log[ sum_t P(X=t|Zp_i) P(Zo_i|X=t) P(Y_i|X=t) ]
    distal.llik <- -s3.distal$res$value

    # log P(Zo_i | X=t) at converged mu_hat: N_dis x iT
    log_pZo_t <- p.zx(distal_par)

    Y_dis <- Y.obs[keep_step3_Zo_in_Y, , drop = FALSE]
    mDes_dis <- if (!is.null(mDesign)) {
      mDesign[keep_step3_Zo_in_Y, , drop = FALSE]
    } else {
      NULL
    }

    total.llik.dis <- joint_log_lik_distal(
      Y = Y_dis,
      mPhi = expand_Phi(fit0$mPhi, ivItemcat),
      log_pZo_t = log_pZo_t,
      pi_mat = pi_adj, # N x T: P(X=t|Zp_i) or flat vPi
      mDesign = mDes_dis
    )

    n_meas_params <- (iT - 1L) +
      sum(ifelse(ivItemcat == 2L, 1L, ivItemcat - 1L)) * iT
    n_distal_params <- iT
    total.k.dis <- n_meas_params + n_distal_params
    N_dis <- length(keep_step3_Zo)

    s3.distal.list <- list(
      three_step = distal_par,
      three_step_vcov = Sigma.3.distal,
      three_step.llik = distal.llik,
      llik = total.llik.dis,
      AIC = -2 * total.llik.dis + 2 * total.k.dis,
      BIC = -2 * total.llik.dis + total.k.dis * log(N_dis)
    )
  }

  if (!is.null(Zo_mat) && is.null(Z_mat)) {
    out <- s3.distal.list
    out$family <- family
    out$n_classes <- iT
    out$estimator <- if (use.bch) "BCH" else "ML"
    out$posteriors <- s2$p.xy
    out$classifications <- max.col(s2$p.xy)
    class(out) <- c("tseLCA_distal", "tseLCA")
    return(out)
  }
  if (!is.null(Z_mat) && is.null(Zo_mat)) {
    class(s3.covariate) <- c("tseLCA_covariate", "tseLCA")
    return(s3.covariate)
  }

  out <- list(
    covariate = s3.covariate,
    distal = s3.distal.list,
    family = family,
    n_classes = iT,
    estimator = if (use.bch) "BCH" else "ML",
    posteriors = s2$p.xy,
    classifications = max.col(s2$p.xy)
  )
  class(out) <- c("tseLCA_both", "tseLCA")
  return(out)
}

# -- S3 methods for tseLCA objects ----------------------------------------------
#
# Four subclasses:
#   tseLCA_measurement  - measurement model only (no Zp, no Zo)
#   tseLCA_covariate    - covariate model only   (Zp present, no Zo)
#   tseLCA_distal       - distal outcome only    (Zo present, no Zp)
#   tseLCA_both         - both covariate and distal

# -- helpers -------------------------------------------------------------------

#' Format covariate coefficients as a printable data frame
#' @noRd
.covariate_table <- function(x) {
  # x is a tseLCA_covariate or x$covariate for tseLCA_both
  est <- as.vector(x$three_step)
  se <- sqrt(diag(x$three_step_vcov))
  zval <- est / se
  pval <- 2 * pnorm(-abs(zval))
  nms <- rownames(x$three_step_vcov)
  data.frame(
    Estimate = est,
    Std.Error = se,
    z.value = zval,
    p.value = pval,
    row.names = nms,
    check.names = FALSE
  )
}

#' Format distal outcome estimates as a printable data frame
#' @noRd
.distal_table <- function(x, family) {
  # x is a tseLCA_distal or x$distal for tseLCA_both
  est <- x$three_step
  se <- sqrt(diag(x$three_step_vcov))
  zval <- est / se
  pval <- 2 * pnorm(-abs(zval))
  scale_label <- switch(
    family,
    gaussian = "(mean)",
    poisson = "(log mean)",
    binomial = "(logit)",
    "(parameter)"
  )
  data.frame(
    Estimate = est,
    Std.Error = se,
    z.value = zval,
    p.value = pval,
    row.names = paste0(names(est), " ", scale_label),
    check.names = FALSE
  )
}

#' Format p-values with significance stars, fixed width
#' @noRd
.format_pval <- function(p) {
  # All entries formatted to the same width so decimal points align in print():
  #   "< 0.001 ***"   (special case, 11 chars)
  #   "0.0099  ** "   (6-char number + space + 4-char marker)
  #   "0.0526  *  "
  #   "0.0800  .  "
  #   "0.2275     "
  ifelse(
    p < 0.001,
    "< 0.001 ***",
    ifelse(
      p < 0.01,
      sprintf("%.4f  ** ", p),
      ifelse(
        p < 0.05,
        sprintf("%.4f  *  ", p),
        ifelse(p < 0.10, sprintf("%.4f  .  ", p), sprintf("%.4f     ", p))
      )
    )
  )
}

#' Print a coefficient data frame with rounded values and significance codes
#' @noRd
.print_table <- function(df, digits = 4) {
  df_fmt <- df
  df_fmt$Estimate <- round(df$Estimate, digits)
  df_fmt$Std.Error <- round(df$Std.Error, digits)
  df_fmt$z.value <- round(df$z.value, digits)
  df_fmt$p.value <- .format_pval(df$p.value)
  print(df_fmt, quote = FALSE, right = TRUE)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}


# -- print ---------------------------------------------------------------------

#' Print a tseLCA model object
#'
#' Compact one-line or table summary printed to the console.
#'
#' @param x   A `tseLCA` object returned by [tseLCA::three_step()].
#' @param digits Integer. Number of decimal places for coefficient tables.
#' @param ... Further arguments.
#' @return Invisibly returns `x`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' print(fit_m)
#' @export
print.tseLCA_measurement <- function(x, ...) {
  cat("tseLCA -- measurement model\n")
  cat(sprintf(
    "  Classes: %d   Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
    x$n_classes,
    x$llik,
    x$AIC,
    x$BIC
  ))
  if (!is.null(x$R2entr)) {
    cat(sprintf("  Entropy R\u00b2: %.4f\n", x$R2entr))
  }
  invisible(x)
}

#' @rdname print.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' print(fit)
#' }
#' @export
print.tseLCA_covariate <- function(x, digits = 4, ...) {
  est <- if (!is.null(x$estimator)) x$estimator else "ML"
  cat("tseLCA -- three-step covariate model\n")
  cat(sprintf(
    "  Classes: %d   Estimator: %s   Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
    x$n_classes,
    est,
    x$llik,
    x$AIC,
    x$BIC
  ))
  if (!is.null(x$entropy.R2)) {
    cat(sprintf("  Entropy R\u00b2 (covariate-adjusted): %.4f\n", x$entropy.R2))
  }
  cat("\nCovariate coefficients (three-step):\n")
  .print_table(.covariate_table(x), digits = digits)
  invisible(x)
}

#' @rdname print.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' print(fit)
#' }
#' @export
print.tseLCA_distal <- function(x, digits = 4, ...) {
  fam <- if (!is.null(x$family)) x$family else "gaussian"
  est <- if (!is.null(x$estimator)) x$estimator else "ML"
  cat("tseLCA -- three-step distal outcome model\n")
  cat(sprintf(
    "  Classes: %d   Estimator: %s   Family: %s\n",
    x$n_classes,
    est,
    fam
  ))
  if (!is.null(x$llik)) {
    cat(sprintf(
      "  Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
      x$llik,
      x$AIC,
      x$BIC
    ))
  }
  cat("\nDistal outcome means by class:\n")
  .print_table(.distal_table(x, fam), digits = digits)
  invisible(x)
}

#' @rdname print.tseLCA_measurement
#' @export
print.tseLCA_both <- function(x, digits = 4, ...) {
  fam <- if (!is.null(x$family)) x$family else "gaussian"
  est <- if (!is.null(x$estimator)) x$estimator else "ML"
  cat("tseLCA -- three-step model with covariate and distal outcome\n")
  cat(sprintf(
    "  Classes: %d   Estimator: %s   Family: %s\n",
    x$n_classes,
    est,
    fam
  ))
  cat(sprintf(
    "  Log-lik: %.4f   AIC: %.2f   BIC: %.2f\n",
    x$covariate$llik,
    x$covariate$AIC,
    x$covariate$BIC
  ))
  cat("\nCovariate coefficients (three-step):\n")
  .print_table(.covariate_table(x$covariate), digits = digits)
  cat("\nDistal outcome means by class:\n")
  .print_table(.distal_table(x$distal, fam), digits = digits)
  invisible(x)
}


# -- summary -------------------------------------------------------------------

#' Summarize a tseLCA model object
#'
#' Verbose summary including model fit, class prevalences, item-response
#' probabilities, and coefficient tables with standard errors and p-values.
#'
#' @param object A `tseLCA` object returned by [tseLCA::three_step()].
#' @param digits Integer. Number of decimal places for coefficient tables.
#' @param ... Further arguments (currently unused).
#' @return Invisibly returns `object`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' summary(fit_m)
#' @export
summary.tseLCA_measurement <- function(object, ...) {
  cat("-- tseLCA Measurement Model --------------------------------\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Log-likelihood : %.4f\n", object$llik))
  cat(sprintf("AIC            : %.4f\n", object$AIC))
  cat(sprintf("BIC            : %.4f\n", object$BIC))
  if (!is.null(object$R2entr)) {
    cat(sprintf("Entropy R\u00b2     : %.4f\n", object$R2entr))
  }
  cat("\nClass prevalences:\n")
  vPi <- object$measurement_model$fit0$vPi
  names(vPi) <- paste0("C", seq_along(vPi))
  print(round(vPi, 4))
  cat("\nItem-response probabilities (P(Y=1|class)):\n")
  print(round(object$measurement_model$fit0$mPhi, 4))
  invisible(object)
}

#' @rdname summary.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' }
#' @export
summary.tseLCA_covariate <- function(object, digits = 4, ...) {
  est <- if (!is.null(object$estimator)) object$estimator else "ML"
  cat("-- tseLCA Three-step Covariate Model -----------------------\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Estimator      : %s\n", est))
  cat(sprintf("Log-likelihood : %.4f\n", object$llik))
  cat(sprintf("AIC            : %.4f\n", object$AIC))
  cat(sprintf("BIC            : %.4f\n", object$BIC))
  if (!is.null(object$entropy.R2)) {
    cat(sprintf(
      "Entropy R\u00b2     : %.4f  (covariate-adjusted)\n",
      object$entropy.R2
    ))
  }

  if (!is.null(object$two_step)) {
    ts <- object$two_step
    if (nrow(ts) > 0L && (is.null(rownames(ts)) || rownames(ts)[1L] == "")) {
      rownames(ts)[1L] <- "Intercept"
    }
    cat("\nTwo-step (starting) estimates:\n")
    print(round(ts, digits))
  }

  cat("\nThree-step estimates:\n")
  .print_table(.covariate_table(object), digits = digits)

  invisible(object)
}

#' @rdname summary.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' summary(fit)
#' }
#' @export
summary.tseLCA_distal <- function(object, digits = 4, ...) {
  fam <- if (!is.null(object$family)) object$family else "gaussian"
  est <- if (!is.null(object$estimator)) object$estimator else "ML"
  cat("-- tseLCA Three-step Distal Outcome Model -------------------\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Estimator      : %s\n", est))
  cat(sprintf("Family         : %s\n", fam))
  if (!is.null(object$llik)) {
    cat(sprintf("Log-likelihood : %.4f\n", object$llik))
    cat(sprintf("AIC            : %.4f\n", object$AIC))
    cat(sprintf("BIC            : %.4f\n", object$BIC))
  }
  cat("\nDistal outcome estimates by class:\n")
  .print_table(.distal_table(object, fam), digits = digits)
  invisible(object)
}

#' @rdname summary.tseLCA_measurement
#' @export
summary.tseLCA_both <- function(object, digits = 4, ...) {
  fam <- if (!is.null(object$family)) object$family else "gaussian"
  est <- if (!is.null(object$estimator)) object$estimator else "ML"
  cat("-- tseLCA Three-step Model: Covariate + Distal Outcome -----\n")
  cat(sprintf("Latent classes : %d\n", object$n_classes))
  cat(sprintf("Estimator      : %s\n", est))
  cat(sprintf("Family         : %s\n", fam))
  cat(sprintf("Log-likelihood : %.4f\n", object$covariate$llik))
  cat(sprintf("AIC            : %.4f\n", object$covariate$AIC))
  cat(sprintf("BIC            : %.4f\n", object$covariate$BIC))

  if (!is.null(object$covariate$two_step)) {
    ts <- object$covariate$two_step
    if (nrow(ts) > 0L && (is.null(rownames(ts)) || rownames(ts)[1L] == "")) {
      rownames(ts)[1L] <- "Intercept"
    }
    cat("\nCovariate -- two-step (starting) estimates:\n")
    print(round(ts, digits))
  }

  cat("\nCovariate -- three-step estimates:\n")
  .print_table(.covariate_table(object$covariate), digits = digits)

  cat("\nDistal outcome -- three-step estimates:\n")
  .print_table(.distal_table(object$distal, fam), digits = digits)

  invisible(object)
}


# -- coef ----------------------------------------------------------------------

#' Extract coefficients from a tseLCA model object
#'
#' @param object A `tseLCA` object returned by [tseLCA::three_step()].
#' @param which Character. For covariate and both models: `"three_step"`
#'   (default) or `"two_step"`. For both models also accepts `"covariate"`,
#'   `"distal"`, or `"both"`.
#' @param step  Character. For `tseLCA_both`: `"three_step"` (default) or
#'   `"two_step"`.
#' @param ... Further arguments (currently unused).
#' @return The coefficient matrix (covariate models), named numeric vector
#'   (distal models), or a named list of both (measurement or both models).
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' coef(fit_m)   # returns list with $prevalences and $item_probs
#' @export
coef.tseLCA_measurement <- function(object, ...) {
  list(
    prevalences = object$measurement_model$fit0$vPi,
    item_probs = object$measurement_model$fit0$mPhi
  )
}

#' @rdname coef.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' coef(fit)                      # three-step estimates
#' coef(fit, which = "two_step")  # two-step starting values
#' }
#' @export
coef.tseLCA_covariate <- function(
  object,
  which = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  if (which == "two_step") {
    if (is.null(object$two_step)) {
      stop("No two-step estimates available in this object.", call. = FALSE)
    }
    return(object$two_step)
  }
  object$three_step
}

#' @rdname coef.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' coef(fit)   # named vector of class means
#' }
#' @export
coef.tseLCA_distal <- function(object, ...) {
  object$three_step
}

#' @rdname coef.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", Zo.name = "Zo",
#'                   use.simple.cov = TRUE)
#' coef(fit, which = "covariate")
#' coef(fit, which = "distal")
#' coef(fit, which = "both")
#' }
#' @export
coef.tseLCA_both <- function(
  object,
  which = c("both", "covariate", "distal"),
  step = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  step <- match.arg(step)
  cov_coef <- if (step == "two_step") {
    if (is.null(object$covariate$two_step)) {
      stop("No two-step covariate estimates available.", call. = FALSE)
    }
    object$covariate$two_step
  } else {
    object$covariate$three_step
  }
  if (which == "covariate") {
    return(cov_coef)
  }
  if (which == "distal") {
    return(object$distal$three_step)
  }
  list(covariate = cov_coef, distal = object$distal$three_step)
}


# -- vcov ----------------------------------------------------------------------

#' Extract the variance-covariance matrix from a tseLCA model object
#'
#' For measurement models, returns the BHHH variance-covariance matrix in
#' the unconstrained log-ratio parameterization (NOT the probability scale).
#' Row and column names identify each parameter as
#' `log(pi_t/pi_1)` (class prevalences) or
#' `log(P(Y=k|C_t)/P(Y=0|C_t))` (item-response probabilities).
#' An attribute `"parameterization"` is attached to remind the user of the
#' scale.
#'
#' @param object A `tseLCA` object returned by [tseLCA::three_step()].
#' @param boundary.tol Scalar. Parameters within this tolerance of 0 or 1
#'   are treated as fixed. Default \code{1e-2}.
#' @param which Character. `"three_step"` (default) or `"two_step"` for
#'   covariate models; `"covariate"`, `"distal"`, or `"both"` for both models.
#' @param step  Character. For `tseLCA_both`: `"three_step"` (default) or
#'   `"two_step"`.
#' @param ... Further arguments (currently unused).
#' @return A named square matrix in the unconstrained log-ratio
#'   parameterization. Row/column names identify each parameter as
#'   `log(pi_t/pi_1)` or `log(P(Y=k|C_t)/P(Y=0|C_t))`. An attribute
#'   `"parameterization"` is attached as a reminder. Returns `NULL`
#'   invisibly if `fit0$mU` is not available. For structural models,
#'   returns the Step-3 vcov matrix; the two-step vcov is only available
#'   when `get.twostep.vcov = TRUE`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' V <- vcov(fit_m)
#' # Names show log-ratio parameterization:
#' rownames(V)
#' attr(V, "parameterization")
#' @export
vcov.tseLCA_measurement <- function(object, boundary.tol = 1e-2, ...) {
  fit0 <- object$measurement_model$fit0
  ivItemcat <- object$measurement_model$ivItemcat
  Y.names <- object$measurement_model$Y.names

  if (is.null(fit0)) {
    stop("No fit0 found in measurement_model.", call. = FALSE)
  }

  # ---- Compute Varmat from mU if available, else return NULL -----------------
  # Use stored ref_idx to reorder u_post columns to match rebased fit0.
  ref_idx <- if (!is.null(object$measurement_model$ref_idx)) {
    object$measurement_model$ref_idx
  } else {
    1L
  }

  step1_Y <- if (!is.null(fit0$mU)) {
    raw <- extract_Y_from_mU(fit0, ivItemcat)
    if (ref_idx != 1L) {
      T_ <- length(fit0$vPi)
      ord <- c(ref_idx, seq_len(T_)[-ref_idx])
      raw$u_post <- raw$u_post[, ord, drop = FALSE]
    }
    raw
  } else {
    message(
      "vcov.tseLCA_measurement: fit0$mU not found. ",
      "Re-estimate with multilevLCA to obtain the measurement vcov."
    )
    return(invisible(NULL))
  }

  V <- lca_indiv_varmat(
    step1_Y$Y.exp,
    step1_Y$mDesign,
    fit0,
    step1_Y$ivItemcat,
    boundary.tol = boundary.tol,
    u_post = step1_Y$u_post
  )$Varmat

  # ---- Build parameter names -------------------------------------------------
  # Parameter ordering matches lca_indiv_varmat:
  #   Pi block   (T-1 params): log(pi_t / pi_1),  t = 2..T
  #   Phi block  (n_free_phi * T params, class-major):
  #     for t=1..T: for each item h, log(P(Y=k|C_t) / P(Y=0|C_t)), k=1..K_h-1
  #
  # Note: parameters are log-ratio transforms of the probability parameters,
  # NOT the probabilities themselves.

  iT <- length(fit0$vPi)
  H <- length(ivItemcat)
  classes <- paste0("C", seq_len(iT))

  # Pi names: log(pi_t/pi_1) for t=2..T
  pi_names <- paste0("log(pi_", classes[-1L], "/pi_", classes[1L], ")")

  # Phi names per item: log(P(Y=k|C_t)/P(Y=0|C_t)) for free k, all classes
  item_labels <- if (!is.null(Y.names)) Y.names else paste0("Y", seq_len(H))

  phi_names <- character(0L)
  for (t in seq_len(iT)) {
    for (h in seq_len(H)) {
      K_h <- ivItemcat[h]
      for (k in seq_len(K_h - 1L)) {
        phi_names <- c(
          phi_names,
          sprintf(
            "log(P(%s=%d|%s)/P(%s=0|%s))",
            item_labels[h],
            k,
            classes[t],
            item_labels[h],
            classes[t]
          )
        )
      }
    }
  }

  all_names <- c(pi_names, phi_names)
  rownames(V) <- all_names
  colnames(V) <- all_names

  attr(V, "parameterization") <-
    "log-ratio (unconstrained); NOT probabilities"

  V
}

#' @rdname vcov.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' vcov(fit)   # Q*(T-1) x Q*(T-1) vcov matrix with named rows/cols
#' }
#' @export
vcov.tseLCA_covariate <- function(
  object,
  which = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  if (which == "two_step") {
    if (is.null(object$two_step_vcov)) {
      stop(
        "No two-step vcov available. Set get.twostep.vcov = TRUE in three_step().",
        call. = FALSE
      )
    }
    return(object$two_step_vcov)
  }
  object$three_step_vcov
}

#' @rdname vcov.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "distal", seed = 2)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zo.name = "Zo", use.simple.cov = TRUE)
#' vcov(fit)   # T x T vcov matrix with mu_C1..mu_CT row/col names
#' }
#' @export
vcov.tseLCA_distal <- function(object, ...) {
  object$three_step_vcov
}

#' @rdname vcov.tseLCA_measurement
#' @examples
#' \donttest{
#' d   <- generate_data(200, "high", "covariate", seed = 1)
#' d$Zo <- rnorm(200, mean = c(-1, 0, 1)[d$X], sd = 0.5)
#' fit <- three_step(d, paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", Zo.name = "Zo",
#'                   use.simple.cov = TRUE)
#' vcov(fit, which = "covariate")
#' vcov(fit, which = "distal")
#' }
#' @export
vcov.tseLCA_both <- function(
  object,
  which = c("both", "covariate", "distal"),
  step = c("three_step", "two_step"),
  ...
) {
  which <- match.arg(which)
  step <- match.arg(step)
  cov_vcov <- if (step == "two_step") {
    if (is.null(object$covariate$two_step_vcov)) {
      stop(
        "No two-step vcov available. Set get.twostep.vcov = TRUE in three_step().",
        call. = FALSE
      )
    }
    object$covariate$two_step_vcov
  } else {
    object$covariate$three_step_vcov
  }
  if (which == "covariate") {
    return(cov_vcov)
  }
  if (which == "distal") {
    return(object$distal$three_step_vcov)
  }
  list(covariate = cov_vcov, distal = object$distal$three_step_vcov)
}

# -- plot methods --------------------------------------------------------------
#
# All four subclasses delegate to multilevLCA's plot.multiLCA, which draws
# item-response probability profiles across classes from the Step-1 fit.
# Extra arguments (horiz, clab, ...) are passed straight through.

#' Extract the raw multiLCA fit0 object from any tseLCA subclass
#' @noRd
.get_fit0 <- function(x) {
  # Extract the raw multiLCA fit0 object from any tseLCA subclass.
  if (inherits(x, "tseLCA_both")) {
    x$covariate$measurement_model$fit0
  } else if (inherits(x, c("tseLCA_covariate", "tseLCA_distal"))) {
    x$measurement_model$fit0
  } else {
    # tseLCA_measurement
    x$measurement_model$fit0
  }
}

#' Plot item-response probability profiles for a tseLCA model
#'
#' Delegates to `plot.multiLCA` from \pkg{multilevLCA}, which draws the
#' class-specific item-response probability profiles from the Step-1
#' measurement model.
#'
#' @param x    A `tseLCA` object returned by [tseLCA::three_step()].
#' @param horiz Logical. If `TRUE`, item labels are drawn horizontally.
#' @param clab  Optional character vector of length T giving class labels.
#' @param ...  Further arguments passed to `plot.multiLCA`.
#'
#' @return Called for its side effect (a base-graphics plot). Invisibly
#'   returns `NULL`.
#' @examples
#' d    <- generate_data(100, "high", "covariate", seed = 1)
#' fit_m <- three_step(d, paste0("Y", 1:6), n_classes = 3)
#' plot(fit_m)
#'
#' \donttest{
#' # Custom class labels
#' plot(fit_m, clab = c("Low risk", "Mixed", "High risk"))
#' }
#' @export
plot.tseLCA_measurement <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

#' @rdname plot.tseLCA_measurement
#' @export
plot.tseLCA_covariate <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

#' @rdname plot.tseLCA_measurement
#' @export
plot.tseLCA_distal <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}

#' @rdname plot.tseLCA_measurement
#' @export
plot.tseLCA_both <- function(x, horiz = FALSE, clab = NULL, ...) {
  plot(.get_fit0(x), horiz = horiz, clab = clab, ...)
  invisible(NULL)
}
