# tseLCA/R/distal-ml.R
#
# Three-step ML estimation for distal outcomes (Bakk, Tekle & Vermunt 2013).
#
# Log-likelihood over the expanded data (one record per person i and assigned
# class s, weighted by w_is = P(W = s | Y_i)):
#
#   l = sum_i sum_s w_is log M_is,   M_is = sum_t a_it f(z_i | t) P(W = s | X = t),
#
# where a_it is the class prior: P(X = t), or P(X = t | Zp_i) when covariates
# are also modeled. With modal assignment (a single w_is = 1 per person) this
# is Eq. 12 of Bakk et al.; with proportional assignment the weights stay
# outside the log, as for the covariate model (Vermunt 2010).
#
# Per-record posteriors R_ist = a_it f(z_i | t) P(W = s | t) / M_is give the
# E-step weights lambda_it = sum_s w_is R_ist.

#' Per-record posteriors and log M_is of the three-step ML distal model
#'
#' @param log_f n x T matrix of log f(z_i | X = t).
#' @param a n x T matrix of class priors a_it.
#' @param pwx T x T classification-error matrix, pwx\[s, t\] = P(W = s | X = t).
#' @return list(R = list of T matrices (n x T), R\[\[s\]\]\[i, t\] = R_ist;
#'   logM = n x T matrix of log M_is).
#' @noRd
distal_records <- function(log_f, a, pwx) {
  row_max <- apply(log_f, 1L, max)
  af <- a * exp(log_f - row_max)
  S <- nrow(pwx)
  R <- vector("list", S)
  logM <- matrix(0, nrow(log_f), S)
  for (s in seq_len(S)) {
    num <- sweep(af, 2L, pwx[s, ], "*")
    M <- pmax(rowSums(num), 1e-300)
    R[[s]] <- num / M
    logM[, s] <- log(M) + row_max
  }
  list(R = R, logM = logM)
}

#' E-step weights lambda_it = sum_s w_is R_ist
#' @noRd
distal_lambda <- function(R, w.is) {
  Reduce(`+`, lapply(seq_along(R), function(s) R[[s]] * w.is[, s]))
}

#' Log-likelihood sum_i sum_s w_is log M_is
#' @noRd
distal_loglik <- function(logM, w.is) {
  keep <- w.is > 0
  sum(w.is[keep] * logM[keep])
}

#' Unit scores and second derivatives of log f(z_i | X = t)
#'
#' Parameters theta:
#' * gaussian: (mu_1, ..., mu_T, sigma2);
#' * poisson / binomial: class linear predictors (log means / logits);
#' * multinomial: the T x C class-conditional category probabilities,
#'   column-major (index (c - 1) * T + t), with z holding categories 1..C.
#'   Its unit score is 1(z_i = c) - pi_tc (the gradient in the softmax
#'   parameterization, as used throughout); the Hessian is obtained from the
#'   Jacobian of the estimating equation (distal_multinomial_jacobian).
#'
#' @return list(G = list of P matrices (n x T), G\[\[p\]\]\[i, t\] =
#'   d log f(z_i | t) / d theta_p; D2 = function(p, q) returning the n x T
#'   matrix of second derivatives, or NULL for multinomial).
#' @noRd
distal_unit_derivs <- function(theta, z, iT, family, C = NULL) {
  N <- length(z)
  in_class <- function(t, v) {
    m <- matrix(0, N, iT)
    m[, t] <- v
    m
  }
  if (family == "gaussian") {
    mu <- theta[1:iT]
    s2 <- theta[iT + 1L]
    resid <- outer(z, mu, "-")
    P <- iT + 1L
    G <- c(
      lapply(seq_len(iT), function(t) in_class(t, resid[, t] / s2)),
      list((resid^2 - s2) / (2 * s2^2))
    )
    D2 <- function(p, q) {
      if (p <= iT && q <= iT) {
        return(if (p == q) in_class(p, -1 / s2) else matrix(0, N, iT))
      }
      if (p == P && q == P) {
        return(1 / (2 * s2^2) - resid^2 / s2^3)
      }
      t <- min(p, q)
      in_class(t, -resid[, t] / s2^2)
    }
    return(list(G = G, D2 = D2))
  }
  if (family == "multinomial") {
    pi_hat <- matrix(theta, nrow = iT, ncol = C)
    G <- vector("list", iT * C)
    for (c in seq_len(C)) {
      for (t in seq_len(iT)) {
        G[[(c - 1L) * iT + t]] <- in_class(t, as.numeric(z == c) - pi_hat[t, c])
      }
    }
    return(list(G = G, D2 = NULL))
  }
  # canonical-link families: d log f / d eta_t = z - E[z | t],
  # d2 log f / d eta_t2 = -Var(z | t)
  eta <- theta[1:iT]
  if (family == "poisson") {
    m <- exp(eta)
    v <- m
  } else {
    m <- 1 / (1 + exp(-eta))
    v <- m * (1 - m)
  }
  g <- outer(z, m, "-")
  G <- lapply(seq_len(iT), function(t) in_class(t, g[, t]))
  D2 <- function(p, q) if (p == q) in_class(p, -v[p]) else matrix(0, N, iT)
  list(G = G, D2 = D2)
}

#' Case-wise score (n x P): sum_t lambda_it G^p_it
#' @noRd
distal_score <- function(lambda, G) {
  vapply(G, function(g) rowSums(lambda * g), numeric(nrow(lambda)))
}

#' K^p_is = sum_t R_ist G^p_it for every assigned class s (list of n x P)
#' @noRd
distal_record_scores <- function(R, G) {
  lapply(R, function(Rs) distal_score(Rs, G))
}

#' Hessian of the negative ML distal log-likelihood
#'
#'   d2 l = sum_i sum_t lambda_it (d2 log f_it + G_it G_it')
#'          - sum_i sum_s w_is K_is K_is'
#' @noRd
distal_neg_hessian <- function(rec, w.is, derivs) {
  G <- derivs$G
  P <- length(G)
  lambda <- distal_lambda(rec$R, w.is)
  H <- matrix(0, P, P)
  for (p in seq_len(P)) {
    for (q in p:P) {
      H[p, q] <- H[q, p] <- sum(lambda * (derivs$D2(p, q) + G[[p]] * G[[q]]))
    }
  }
  K <- distal_record_scores(rec$R, G)
  for (s in seq_along(K)) {
    H <- H - crossprod(K[[s]], K[[s]] * w.is[, s])
  }
  -H
}

#' Jacobian of the multinomial ML distal estimating equation
#'
#' Psi_(t,c) = sum_i lambda_it (1(z_i = c) - pi_tc), differentiated with
#' respect to pi_(t',c'):
#'   (1(c = c') - pi_tc) / pi_t'c' *
#'     \[1(t = t') sum_{i: z_i = c'} lambda_it
#'      - sum_{i: z_i = c'} sum_s w_is R_ist R_ist']
#'   - 1(t = t', c = c') sum_i lambda_it.
#' @noRd
distal_multinomial_jacobian <- function(pi_hat, rec, w.is, Y_cat) {
  iT <- nrow(pi_hat)
  C <- ncol(pi_hat)
  n_par <- iT * C
  lambda <- distal_lambda(rec$R, w.is)
  Jac <- matrix(0, n_par, n_par)
  lambda_colsums <- colSums(lambda)
  for (cprime in seq_len(C)) {
    idx_i <- which(Y_cat == cprime)
    if (length(idx_i) == 0L) {
      next
    }
    R1 <- colSums(lambda[idx_i, , drop = FALSE])
    R2 <- Reduce(
      `+`,
      lapply(seq_along(rec$R), function(s) {
        Rs <- rec$R[[s]][idx_i, , drop = FALSE]
        crossprod(Rs, Rs * w.is[idx_i, s])
      })
    )
    for (tprime in seq_len(iT)) {
      col_idx <- (cprime - 1L) * iT + tprime
      denom <- pmax(pi_hat[tprime, cprime], 1e-300)
      for (c in seq_len(C)) {
        for (t in seq_len(iT)) {
          row_idx <- (c - 1L) * iT + t
          bracket <- (if (t == tprime) R1[t] else 0) - R2[t, tprime]
          term1 <- (as.numeric(cprime == c) - pi_hat[t, c]) * bracket / denom
          term2 <- if (t == tprime && c == cprime) lambda_colsums[t] else 0
          Jac[row_idx, col_idx] <- term1 - term2
        }
      }
    }
  }
  Jac
}

#' Cross-derivatives of the ML distal estimating equation with respect to
#' the classification-error parameters (Step-1 term) and the covariate
#' coefficients of the class prior (Step-2 term)
#'
#' theta2 are the off-diagonal log-ratios of pwx (column softmax), ordered
#' t0 then s0 != t0. With c_s = 1(s = s0) - pwx\[s0, t0\]:
#'   C1\[p, (s0, t0)\] = sum_i sum_s w_is c_s R_ist0 (G^p_it0 - K^p_is).
#' For gamma_(l, q) of a_it = P(X = t | Zp_i) (multinomial logit, class l + 1
#' against the reference) the prior's own derivative cancels:
#'   C\[p, (l, q)\] = sum_i z_iq sum_s w_is R_is,l+1 (G^p_i,l+1 - K^p_is).
#' @return list(C1 = P x T(T-1) matrix, C_mat = P x (Q+1)(T-1) matrix or NULL).
#' @noRd
distal_cross_derivs <- function(rec, w.is, pwx, G, Z_mat = NULL) {
  iT <- ncol(pwx)
  S <- nrow(pwx)
  P <- length(G)
  K <- distal_record_scores(rec$R, G)
  resid_s <- function(s, t) {
    # N x P matrix of G^p_it - K^p_is
    vapply(
      seq_len(P),
      function(p) G[[p]][, t] - K[[s]][, p],
      numeric(nrow(w.is))
    )
  }

  C1 <- matrix(0, P, iT * (iT - 1L))
  col <- 0L
  for (t0 in seq_len(iT)) {
    for (s0 in seq_len(iT)) {
      if (s0 == t0) {
        next
      }
      col <- col + 1L
      for (s in seq_len(S)) {
        c_s <- as.numeric(s == s0) - pwx[s0, t0]
        wr <- w.is[, s] * rec$R[[s]][, t0] * c_s
        C1[, col] <- C1[, col] + colSums(wr * resid_s(s, t0))
      }
    }
  }

  C_mat <- NULL
  if (!is.null(Z_mat)) {
    Q <- ncol(Z_mat)
    C_mat <- matrix(0, P, (iT - 1L) * Q)
    for (l in seq_len(iT - 1L)) {
      idx <- ((l - 1L) * Q + 1L):(l * Q)
      for (s in seq_len(S)) {
        wr <- w.is[, s] * rec$R[[s]][, l + 1L]
        C_mat[, idx] <- C_mat[, idx] + crossprod(wr * resid_s(s, l + 1L), Z_mat)
      }
    }
  }
  list(C1 = C1, C_mat = C_mat)
}
