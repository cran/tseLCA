# tseLCA/R/vcov.R
#
# Variance estimation: the Step-1 measurement-model variance and the Step-3
# sandwich variances with Step-1 (and Step-2 covariate) uncertainty
# propagation.

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
#'   \item Dichotomous item k (\code{ivItemcat[k] == 2}): 1 row =
#'     \eqn{P(Y=1|C)}; the base level \eqn{P(Y=0|C)} is excluded.
#'   \item Polytomous item k (\code{ivItemcat[k] > 2}): \code{R_k} rows =
#'     \eqn{P(Y=0|C), \ldots, P(Y=R_k-1|C)}; the base level is included.
#' }
#' \code{expand_Y} produces one-hot columns in the same order so that
#' \code{expand_Phi(fit0$mPhi, ivItemcat)} aligns column-wise with
#' \code{expand_Y(mY, ivItemcat)}.  Free (estimable) parameters per item are
#' the single \eqn{P(Y=1|C)} row for dichotomous items, and rows 2 through
#' \eqn{R_k} for polytomous items (row 1, \eqn{P(Y=0|C)}, is the reference).
#' Boundary parameters (within \code{boundary.tol} of 0 or 1) are treated as
#' fixed: their score columns are zeroed and they do not contribute to the
#' information matrix.
#'
#' @param Y.exp       n x sum(R_k) expanded one-hot indicator matrix.
#' @param mDesign.exp Expanded design matrix (same dimensions as \code{Y.exp}),
#'   or \code{NULL} for complete data.
#' @param fit0        Step-1 fit object with \code{$vPi} and \code{$mPhi}.
#' @param ivItemcat   Integer vector of category counts per item.
#' @param boundary.tol Scalar tolerance for boundary detection. Default
#'   \code{1e-2}.
#' @param use.freq    Logical. Collapse duplicate score rows before computing
#'   the cross-product, weighting by frequency. Default \code{TRUE}.
#' @param u_post      Optional n x T matrix of posterior class probabilities.
#'   When supplied (e.g. extracted from \code{fit0$mU} with
#'   \code{extract_Y_from_mU}), \code{compute_posteriors} is skipped.
#'   Default \code{NULL}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`Infomat`}{Square BHHH information matrix of dimension p x p,
#'       where p = (iT-1) + sum(ivItemcat - 1) * iT is the total number of free
#'       parameters. Boundary parameters have zero rows and columns.}
#'     \item{`Varmat`}{Inverse of \code{Infomat} divided by n, giving the
#'       asymptotic variance-covariance matrix on the same scale as
#'       \pkg{multilevLCA}'s \code{$Varmat}. Boundary parameters have zero
#'       rows and columns.}
#'     \item{`SEs`}{Numeric vector of length p. Square root of the diagonal of
#'       \code{Varmat}; zero for boundary parameters.}
#'     \item{`mScore`}{n x p matrix of individual score contributions in the
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

  # ---- Reference categories on the boundary ----------------------------------
  # The free parameters of a polytomous item are log-ratios against its first
  # category. When P(Y = first | class t) is on the boundary, the scores of that
  # item's free parameters in class t sum to (almost) zero: a common shift of
  # its log-ratios, i.e. log P(Y = first | t), is not informed by the data.
  # Treat that direction as fixed, like other boundary parameters: invert the
  # information on its orthogonal complement.
  null_dirs <- list()
  for (h in which(ivItemcat > 2L)) {
    rows_h <- which(free_idx >= starts[h] & free_idx < starts[h] + ivItemcat[h])
    for (t in seq_len(iT)) {
      if (!phi_bdry[starts[h], t]) next
      cols <- (iT - 1L) + (t - 1L) * n_free_phi + rows_h
      cols <- match(cols, active)
      cols <- cols[!is.na(cols)]
      if (length(cols) < 1L) next
      v <- numeric(length(active))
      v[cols] <- 1
      null_dirs[[length(null_dirs) + 1L]] <- v
    }
  }
  invert_info <- function(I) {
    if (length(null_dirs) == 0L) {
      return(qr.solve(I))
    }
    Nd <- do.call(cbind, null_dirs)
    B <- qr.Q(qr(Nd), complete = TRUE)[, -seq_len(qr(Nd)$rank), drop = FALSE]
    V <- B %*% qr.solve(crossprod(B, I %*% B)) %*% t(B)
    (V + t(V)) / 2
  }

  Varmat_active <- tryCatch(
    invert_info(Infomat_active) / N,
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
  use.bch,
  unit_scores = NULL
) {
  # mu_hat is the full Step-3 parameter vector: the T class parameters, plus
  # sigma2 for the gaussian ML model. `unit_scores(mu_hat)`, when supplied,
  # returns a list with one N x T matrix per parameter p holding
  # d log f(z_i | X = t) / d theta_p; otherwise each class parameter only
  # moves its own class and the unit scores are recovered from the score.
  J.3 <- three_step.score(mu_hat)
  meat <- crossprod(J.3)
  P <- ncol(J.3)
  par_names <- c(paste0("mu_C", seq_len(iT)), if (P > iT) "sigma2")

  if (use.bch || use.simple.cov) {
    result <- H.3.inv %*% meat %*% H.3.inv
    dimnames(result) <- list(par_names, par_names)
    return(result)
  }

  if (is.null(unit_scores)) {
    stop("lca_vcov_distal(): `unit_scores` is required for ML uncertainty propagation.")
  }
  # Per-record posteriors and unit scores at the estimates; the Step-1 (C1)
  # and Step-2 (C_mat) cross-derivatives are derived in R/distal-ml.R.
  rec <- distal_records(p.zx(mu_hat), pi_adj, p.wx_mat)
  cross <- distal_cross_derivs(
    rec,
    w.is,
    p.wx_mat,
    unit_scores(mu_hat),
    Z_mat = if (!is.null(s3.par) && !is.null(p.xz.cov) && !is.null(Z_mat_cov)) {
      Z_mat_cov
    } else {
      NULL
    }
  )

  step1.uncertainty <- cross$C1 %*%
    s2$J.2 %*%
    Sigma.1 %*%
    t(s2$J.2) %*%
    t(cross$C1)

  step2.uncertainty <- if (!is.null(cross$C_mat)) {
    cross$C_mat %*% Sigma.3 %*% t(cross$C_mat)
  } else {
    matrix(0, P, P)
  }

  result <- H.3.inv %*%
    (meat + step1.uncertainty + step2.uncertainty) %*%
    H.3.inv
  dimnames(result) <- list(par_names, par_names)
  result
}

#' Distal outcome variance-covariance for the multinomial family
#'
#' Multinomial analog of \code{lca_vcov_distal()}: assembles the
#' \code{(iT*C) x (iT*C)} sandwich variance matrix for the flattened T x C
#' class-conditional probability matrix \code{pi_hat}
#' (\code{matrix(theta_hat, nrow = iT, ncol = C)} recovers it), propagating
#' Step-1 measurement uncertainty (\code{C1_mat}/\code{step1.uncertainty})
#' and, when both a covariate and distal model are fitted, Step-3 covariate
#' uncertainty (\code{C_mat}/\code{step2.uncertainty}) when
#' \code{use.bch = FALSE} and \code{use.simple.cov = FALSE}. The
#' generalization from scalar \code{mu_t} (one parameter per class, as in
#' \code{lca_vcov_distal()}) to \code{pi_hat[t, ]} (C parameters per class)
#' only touches the "unit score" \code{g_it}: in place of dividing the
#' length-\code{iT} score by \code{r_it} once, the length-\code{iT*C} score
#' is divided by \code{r_it} replicated across the C categories, since
#' neither chain-rule term (\code{dr}, through \code{d ae/d theta2}; nor
#' \code{inner_t}, through \code{d r_it/d gamma}) depends on the category
#' dimension at all -- only on which class \code{t} a given column belongs
#' to. Both terms were cross-validated against independent numerical
#' differentiation of the case-wise estimating equation (see the "full
#' propagation" tests in test-integration.R).
#' @noRd
lca_vcov_distal_multinomial <- function(
  theta_hat,
  three_step.score,
  pi_adj,
  w.is,
  p.wx_mat,
  Y_cat,
  C,
  H.3.inv,
  Sigma.1,
  s2,
  iT,
  use.simple.cov,
  use.bch,
  Sigma.3 = NULL,
  s3.par = NULL,
  p.xz.cov = NULL,
  Z_mat_cov = NULL
) {
  J.3 <- three_step.score(theta_hat) # N x (iT*C)
  meat <- crossprod(J.3)

  if (use.bch || use.simple.cov) {
    return(H.3.inv %*% meat %*% H.3.inv)
  }

  pi_hat <- matrix(theta_hat, nrow = iT, ncol = C)
  n_par <- iT * C
  rec <- distal_records(
    log(pmax(t(pi_hat[, Y_cat, drop = FALSE]), 1e-300)),
    pi_adj,
    p.wx_mat
  )
  cross <- distal_cross_derivs(
    rec,
    w.is,
    p.wx_mat,
    distal_unit_derivs(theta_hat, Y_cat, iT, "multinomial", C = C)$G,
    Z_mat = if (!is.null(s3.par) && !is.null(p.xz.cov) && !is.null(Z_mat_cov)) {
      Z_mat_cov
    } else {
      NULL
    }
  )

  step1.uncertainty <- cross$C1 %*% s2$J.2 %*% Sigma.1 %*% t(s2$J.2) %*% t(cross$C1)
  step2.uncertainty <- if (!is.null(cross$C_mat)) {
    cross$C_mat %*% Sigma.3 %*% t(cross$C_mat)
  } else {
    matrix(0, n_par, n_par)
  }

  H.3.inv %*% (meat + step1.uncertainty + step2.uncertainty) %*% H.3.inv
}
