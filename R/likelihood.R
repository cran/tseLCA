# tseLCA/R/likelihood.R
#
# Indicator expansion and latent class log-likelihood building blocks shared
# by all three steps.

#' One-hot expand an integer response matrix
#'
#' Converts an n x K matrix of 0-based integer category values into an
#' n x sum(ivItemcat) binary indicator matrix, one column per category per item.
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
#' Returns an n x T matrix where entry `[i, t]` is the conditional log-likelihood
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
#' @param Y       n x sum(R_k) expanded one-hot response matrix.
#' @param Zo      Length-n distal outcome vector.
#' @param mPhi    sum(R_k) x T expanded item-response probability matrix.
#' @param p.zx    n x T matrix of log-densities log P(Zo_i|X=t).
#' @param pi_mat  n x T matrix of class priors P(X=t|Zp_i). If NULL, uses
#'   flat prevalences from the row means of p.zx (not used; vPi supplied).
#' @param vPi     Length-T flat prevalences, used when pi_mat is NULL.
#' @param mDesign n x sum(R_k) design matrix (NULL for complete data).
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
#' n x T soft posterior matrix P(X=t|Y_i).
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
