# tseLCA/R/step2.R
#
# Step 2 (classification): posterior assignment weights, the
# classification-error matrix P(W = s | X = t), and BCH weights.

#' Compute classification-error matrix with optional covariate-adjusted prior
#'
#' Returns posteriors, modal/soft assignments (w.is), and the T x T
#' classification-error probability matrix p.wx_mat = P(W=s|X=t).
#' When pi_adj (n x T) is supplied, uses person-specific class priors from the
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

      # gamma_st = log(sum_i w_is p_it) - log(sum_i w_it p_it), with w the
      # assignment weights: the posteriors themselves (proportional; they
      # depend on the parameters too) or the modal assignment (held fixed).
      # With dp_it/du_c = p_it (1(t = c) - p_ic) g_ic, where g_ic is 1 for a
      # class-size parameter and (Y - phi) for an item parameter of class c,
      # d gamma_st / du_c = sum_i Q_stc,i g_ic below.
      w_ik <- if (use.modal.assignment) {
        w <- matrix(0, N, T_classes)
        w[cbind(seq_len(N), max.col(p_ik))] <- 1
        w
      } else {
        p_ik
      }
      A <- t(w_ik) %*% p_ik
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

      n_free_phi <- sum(ivItemcat - 1L)
      for (t in seq_len(T_classes)) {
        P_tt <- (w_ik[, t] * p_ik[, t]) / A[t, t]

        for (s in seq_len(T_classes)) {
          if (s == t) {
            next
          }
          row_J <- st_map[s, t]
          P_st <- (w_ik[, s] * p_ik[, t]) / A[s, t]

          for (c_prime in seq_len(T_classes)) {
            I_s <- if (s == c_prime) 1.0 else 0.0
            I_t <- if (t == c_prime) 1.0 else 0.0

            #shared derivative component for class c_prime
            Q_stc <- if (use.modal.assignment) {
              (P_st - P_tt) * (I_t - p_ik[, c_prime])
            } else {
              P_st * (I_s + I_t - 2 * p_ik[, c_prime]) -
                2 * P_tt * (I_t - p_ik[, c_prime])
            }

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

              # columns in the order of the Step-1 variance: class by class,
              # item by item within a class
              col_J_start <- L_rho + (c_prime - 1L) * n_free_phi + phi_row_start - 1L

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

#' BCH classification-error-corrected weight matrix
#'
#' Computes the n x T BCH weight matrix used throughout the BCH estimators:
#' \code{w.it = w.is \%*\% t(pwx)^-1}, where \code{pwx[s, t] = P(W = s | X =
#' t)} is the column-stochastic classification-error matrix from
#' \code{compute_pwx_adj()}/\code{lca_step2()} (\code{colSums(pwx) == 1}).
#' @noRd
bch_weight_matrix <- function(w.is, pwx) {
  w.is %*% t(qr.solve(pwx))
}

#' Restrict Step-2 output to the rows used by one Step-3 model
#'
#' Step 2 is estimated on every row with indicator data; each Step-3 model
#' uses only the rows where its structural variables are observed. `rows`
#' indexes those rows within the Step-2 sample. The Step-2 Jacobian (J.2)
#' is recomputed on the subset when uncertainty propagation needs it.
#' @noRd
.subset_step2 <- function(s2, rows, dat, iT) {
  J.2 <- if (!is.null(s2$compute_J_unc)) {
    Y_sub <- dat$Y.obs[rows, , drop = FALSE]
    mDes_sub <- if (!is.null(dat$mDesign)) {
      dat$mDesign[rows, , drop = FALSE]
    } else {
      matrix(1L, nrow(Y_sub), ncol(Y_sub))
    }
    s2$compute_J_unc(
      s2$p.xy[rows, , drop = FALSE],
      Y_sub,
      mDes_sub,
      s2$theta1,
      dat$ivItemcat,
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
    J.2 = J.2,
    w.is = s2$w.is[rows, , drop = FALSE],
    post = if (!is.null(s2$post)) s2$post[rows, , drop = FALSE] else NULL
  )
}

#' Step 2 for three_step(): classification
#'
#' Posterior class probabilities, assignment weights, and the
#' classification-error matrix on the Step-2 sample, plus their restrictions
#' to the covariate-model (`$cov`) and distal-model (`$dis`) rows.
#' @noRd
.step2 <- function(dat, fit0, n_classes, opts) {
  s2 <- lca_step2(
    dat$Y.obs,
    fit0,
    n_classes,
    opts$use.modal.assignment,
    opts$boundary.tol,
    opts$use.simple.cov || opts$use.bch,
    ivItemcat = dat$ivItemcat,
    mDesign = dat$mDesign
  )
  list(
    all = s2,
    cov = if (!is.null(dat$Z_mat)) {
      .subset_step2(s2, dat$keep_step3_Z_in_Y, dat, n_classes)
    },
    dis = if (!is.null(dat$Zo_mat)) {
      .subset_step2(s2, dat$keep_step3_Zo_in_Y, dat, n_classes)
    }
  )
}
