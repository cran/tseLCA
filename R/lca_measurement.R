# tseLCA/R/lca_measurement.R
#
# Step-1 measurement model via multilevLCA::multiLCA.
#
# Exports:
#   lca_step1()          - measurement model fit + optional two-step fitZ
#   fitZ_from_fit0()     - pure-R EM for gamma with mPhi fixed (default fitZ path)
#   fitZ_from_multiLCA() - two-step estimation with multiLCA(fixedpars=1, Z=...) (used when get.twostep.vcov = TRUE in three_step())

# -- lca_step1 -----------------------------------------------------------------

#' Fit the LCA measurement model (Step 1)
#'
#' Estimates the latent class measurement model with \pkg{multilevLCA} and
#' optionally, fixes `mPhi` and estimates covariate effects (two-step
#' initialization) with `fitZ_from_fit0()`.
#'
#' @param data A data.frame containing at minimum the indicator columns.
#' @param Y.names Character vector of item column names.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names, or `NULL`.
#' @param maxIter.measurement Maximum EM iterations before giving up on convergence. Default `5000L`.
#' @param measurement.tol Convergence tolerance. Default `1e-8`.
#' @param covariate.tol Convergence tolerance for the `fitZ` M-step. Default `1e-6`.
#' @param iter.measurement Number of random restarts when entropy R\eqn{^2} is low. Default `10`.
#' @param R2.threshold Entropy R\eqn{^2} below which restarts are triggered. Default `0.7`.
#' @param use.two.step Logical. If `TRUE`, also estimate `fitZ` with `fitZ_from_fit0()` if `Zp.names` is applied. Default `TRUE`.
#' @param estimate.one.step Logical. If `FALSE`, skip the unconditional EM and only compute `fitZ`. Default `TRUE`.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}. Default `FALSE`.
#' @param maxIter.fitZ Maximum EM iterations for `fitZ_from_fit0()`. Default `200`.
#' @param include.intercept Logical. Prepend intercept to covariate design matrix. Default `TRUE`.
#' @param rebase Character or integer specifying the reference latent class.
#'   Use `"C1"`, `"C2"`, etc. or an integer index. Default `"C1"`. The
#'   measurement model is permuted so this class becomes column 1, making it
#'   the reference for all downstream multinomial logit parameterizations.
#' @param verbose Logical. Print progress messages. Default `FALSE`.
#'
#' @return A list with `$fit0` ([multilevLCA::multiLCA()] measurement model) and `$fitZ`
#'   (two-step covariate model from [fitZ_from_fit0()], or `NULL`).
#' @examples
#' \donttest{
#' d <- generate_data(200, "high", "covariate", seed = 1)
#'
#' # Measurement model only
#' s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' s1$fit0$vPi    # estimated class prevalences
#' s1$fit0$mPhi   # item-response probabilities
#'
#' # With two-step covariate initialization
#' s1z <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                  Zp.names = "Zp", use.two.step = TRUE, verbose = TRUE)
#' s1z$fitZ$mGamma   # two-step gamma estimates
#' }
#' @export
lca_step1 <- function(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  maxIter.measurement = 5000L,
  measurement.tol = 1e-8,
  covariate.tol = 1e-6,
  iter.measurement = 10L,
  R2.threshold = 0.70,
  use.two.step = TRUE,
  estimate.one.step = TRUE,
  incomplete = FALSE,
  maxIter.fitZ = 200L,
  include.intercept = TRUE,
  rebase = "C1",
  verbose = FALSE
) {
  run_measurement_fit <- function(extra_args = list()) {
    args <- c(
      list(
        data,
        Y.names,
        n_classes,
        extout = TRUE,
        incomplete = incomplete,
        maxIter = maxIter.measurement,
        tol = measurement.tol,
        verbose = FALSE
      ),
      extra_args
    )
    fit <- do.call(multilevLCA::multiLCA, args)
    if (nrow(fit$LLKSeries) == maxIter.measurement) {
      args$maxIter <- 2L * maxIter.measurement
      fit <- do.call(multilevLCA::multiLCA, args)
      if (verbose) {
        warning(sprintf(
          "Measurement model hit %d iterations; retried with %d. Low separation is likely the cause.",
          maxIter.measurement,
          2L * maxIter.measurement
        ))
      }
      if (nrow(fit$LLKSeries) == 2L * maxIter.measurement) {
        warning(
          "Measurement model still failed to converge even after running more iterations. Consider increasing maxIter.measurement and or measurement.tol"
        )
      }
    }
    fit
  }

  best_fit <- function(initial, run_fn) {
    ll0 <- initial$LLKSeries[nrow(initial$LLKSeries), 1L]
    if (is.null(initial$R2entr) || initial$R2entr >= R2.threshold) {
      return(initial)
    }
    if (verbose) {
      warning(sprintf(
        "Measurement model has low entropy R\u00b2 (%.3f < %.3f). Running %d additional random restarts.",
        initial$R2entr,
        R2.threshold,
        iter.measurement
      ))
    }
    if (iter.measurement > 0L) {
      cands <- lapply(seq_len(iter.measurement), function(r) run_fn())
      cand_lls <- vapply(
        cands,
        function(f) f$LLKSeries[nrow(f$LLKSeries), 1L],
        numeric(1L)
      )
      best_r <- which.max(cand_lls)
      if (cand_lls[best_r] > ll0) {
        if (verbose) {
          message(sprintf(
            "Restart %d improved log-likelihood to %.4f.",
            best_r,
            cand_lls[best_r]
          ))
        }
        cands[[best_r]]
      } else {
        if (verbose) {
          message("No restart improved on the initial measurement model.")
        }
        initial
      }
    } else {
      initial
    }
  }

  fit0 <- if (estimate.one.step) {
    best_fit(initial = run_measurement_fit(), run_fn = run_measurement_fit)
  } else {
    NULL
  }

  #Permute classes so the desired reference is column 1
  if (!is.null(fit0)) {
    ref_idx <- parse_rebase(rebase, n_classes)
    fit0 <- permute_fit0_classes(fit0, ref_idx)
  }

  fitZ <- if (use.two.step && !is.null(Zp.names) && !is.null(fit0)) {
    fitZ_from_fit0(
      fit0 = fit0,
      data = data,
      Y.names = Y.names,
      Zp.names = Zp.names,
      tol = covariate.tol,
      maxIter = maxIter.fitZ,
      incomplete = incomplete,
      include.intercept = include.intercept,
      rebase = rebase,
      verbose = verbose
    )
  } else {
    NULL
  }

  list(fit0 = fit0, fitZ = fitZ)
}


# -- fitZ_from_fit0 ------------------------------------------------------------

#' Estimate covariate effects with measurement parameters fixed (two-step EM)
#'
#' Fixes `mPhi` at `fit0$mPhi` and estimates multinomial logit coefficients
#' `mGamma` (Q x (T-1)) via an EM algorithm with a BFGS M-step.
#'
#' @param fit0 Output of `lca_step1()$fit0`.
#' @param data A data.frame.
#' @param Y.names Character vector of item column names.
#' @param Zp.names Character vector of covariate column names.
#' @param tol Convergence tolerance. Default `1e-6`.
#' @param maxIter Maximum EM iterations. Default `200`.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}. Default `FALSE`.
#' @param include.intercept Logical. Prepend intercept to covariate design matrix. Default `TRUE`.
#' @param rebase Character or integer. Reference class for the multinomial logit
#'   parameterization (e.g. `"C1"`, `"C2"`, or an integer). Default `"C1"`.
#'   Must match the `rebase` used in `lca_step1()` so class column ordering
#'   is consistent.
#' @param starting_val Optional Q x (T-1) starting value matrix for `mGamma`.
#' @param verbose Logical. Print convergence messages. Default `FALSE`.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`mGamma`}{Q x (T-1) numeric matrix of multinomial logit
#'       coefficients, where Q is the number of columns in the covariate design
#'       matrix (including intercept if `include.intercept = TRUE`). Rows are
#'       named by covariate, columns by non-reference class (e.g. `"C2"`,
#'       `"C3"`).}
#'     \item{`mPhi`}{Expanded item parameter matrix (items x classes), fixed at
#'       `fit0$mPhi` throughout estimation.}
#'     \item{`vOmega`}{Length-T vector of marginal class proportions implied by
#'       the final `mGamma`, computed as column means of the fitted class
#'       probability matrix.}
#'     \item{`LLKSeries`}{Single-column matrix of observed-data log-likelihoods,
#'       one row per EM iteration. Useful for diagnosing convergence.}
#'     \item{`converged`}{Logical. `TRUE` if the EM loop exited before
#'       `maxIter` iterations or if the final log-likelihood change was below
#'       `tol`.}
#'     \item{`n_obs`}{Integer. Number of observations used in estimation after
#'       listwise deletion on covariates.}
#'   }
#' @examples
#' \donttest{
#' d  <- generate_data(200, "high", "covariate", seed = 1)
#' s1 <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#'
#' # Estimate two-step gamma with mPhi fixed at Step-1 values
#' fZ <- fitZ_from_fit0(
#'   fit0     = s1$fit0,
#'   data     = d,
#'   Y.names  = paste0("Y", 1:6),
#'   Zp.names = "Zp",
#'   verbose  = TRUE
#' )
#' fZ$mGamma   # Q x (T-1) coefficient matrix
#' fZ$converged
#' }
#' @export
fitZ_from_fit0 <- function(
  fit0,
  data,
  Y.names,
  Zp.names,
  tol = 1e-6,
  maxIter = 200L,
  incomplete = FALSE,
  include.intercept = TRUE,
  rebase = "C1",
  starting_val = NULL,
  verbose = FALSE
) {
  cd <- clean_data(
    data = data,
    Y.names = Y.names,
    Zp.names = Zp.names,
    incomplete = incomplete,
    include.intercept = include.intercept,
    verbose = verbose
  )
  mY <- cd$Y.obs # expanded N_Y x K
  mDesign <- cd$mDesign
  ivItemcat <- cd$ivItemcat
  # For fitZ we need the Z rows that overlap with the Y-kept rows
  mZ <- cd$Z_mat # N_Z x Q, already complete-case

  mY <- mY[cd$keep_step3_Z_in_Y, , drop = FALSE]
  if (!is.null(mDesign)) {
    mDesign <- mDesign[cd$keep_step3_Z_in_Y, , drop = FALSE]
  }

  mPhi <- expand_Phi(fit0$mPhi, ivItemcat)
  iT <- ncol(mPhi)
  iN <- nrow(mY)
  iP <- ncol(mZ)

  #Prevent parameter estimates on the boundary of the support (prevents NA's in posteriors)
  phi_clamped <- pmax(pmin(mPhi, 1 - 1e-10), 1e-10)
  log_p_it <- if (is.null(mDesign)) {
    mY %*% log(phi_clamped)
  } else {
    (mDesign * mY) %*% log(phi_clamped)
  }

  softmax_rows <- function(mat) {
    mat <- mat - apply(mat, 1L, max)
    ex <- exp(mat)
    ex / rowSums(ex)
  }

  gamma <- matrix(0, nrow = iP, ncol = iT - 1L)
  if (!is.null(starting_val)) {
    if (!isTRUE(all.equal(dim(gamma), dim(starting_val)))) {
      warning(sprintf("starting_val dimensions must be %d x %d.", iP, iT - 1L))
    } else {
      gamma <- starting_val
    }
  }

  ll_prev <- -Inf
  LLKSeries <- numeric(0L)

  for (iter in seq_len(maxIter)) {
    eta_full <- cbind(0, mZ %*% gamma)
    pi_mat <- softmax_rows(eta_full)
    log_joint <- log_p_it + log(pi_mat)
    log_marg <- apply(log_joint, 1L, function(row) {
      mx <- max(row)
      mx + log(sum(exp(row - mx)))
    })
    ll_curr <- sum(log_marg)
    LLKSeries <- c(LLKSeries, ll_curr)
    w_mat <- exp(log_joint - log_marg)

    gamma_new <- tryCatch(
      {
        obj <- function(g_vec) {
          g_mat <- matrix(g_vec, nrow = iP, ncol = iT - 1L)
          pi_ <- softmax_rows(cbind(0, mZ %*% g_mat))
          -sum(w_mat * log(pi_))
        }
        gr <- function(g_vec) {
          g_mat <- matrix(g_vec, nrow = iP, ncol = iT - 1L)
          pi_ <- softmax_rows(cbind(0, mZ %*% g_mat))
          resid <- w_mat[, -1L, drop = FALSE] - pi_[, -1L, drop = FALSE]
          -as.vector(t(mZ) %*% resid)
        }
        res <- optim(par = as.vector(gamma), fn = obj, gr = gr, method = "BFGS")
        matrix(res$par, nrow = iP, ncol = iT - 1L)
      },
      error = function(e) {
        warning(
          "fitZ_from_fit0: optim failed at iter ",
          iter,
          ": ",
          conditionMessage(e),
          ". Keeping previous gamma."
        )
        gamma
      }
    )

    if (iter > 1L && abs(ll_curr - ll_prev) < tol) {
      gamma <- gamma_new
      if (verbose) {
        message(sprintf("fitZ EM converged in %d iterations.", iter))
      }
      break
    }
    gamma <- gamma_new
    ll_prev <- ll_curr
  }

  converged <- (length(LLKSeries) < maxIter) ||
    (abs(LLKSeries[length(LLKSeries)] - LLKSeries[length(LLKSeries) - 1L]) <
      tol)
  if (!converged) {
    warning(
      "fitZ_from_fit0: gamma EM did not converge in ",
      maxIter,
      " iterations."
    )
  }

  pi_final <- softmax_rows(cbind(0, mZ %*% gamma))
  vOmega <- colMeans(pi_final)

  rownames(gamma) <- colnames(mZ)
  # Column names reflect the non-reference classes
  # (all classes except the reference, in ascending order)
  ref_idx <- parse_rebase(rebase, iT)
  non_ref_classes <- seq_len(iT)[-ref_idx]
  colnames(gamma) <- paste0("C", non_ref_classes)

  list(
    mGamma = gamma,
    mPhi = mPhi,
    vOmega = vOmega,
    LLKSeries = matrix(LLKSeries, ncol = 1L),
    converged = converged,
    n_obs = iN
  )
}


# -- fitZ_from_multiLCA --------------------------------------------------------

#' Estimate two-step covariate model via multilevLCA (optional reference path)
#'
#' Calls `multilevLCA::multiLCA` with `fixedpars = 1` and `Z = Zp.names` to
#' fit the two-step covariate model.  This is the original multilevLCA approach
#' and is used when `get.twostep.vcov = TRUE` in [tseLCA::three_step()] to obtain
#' multilevLCA's corrected standard errors for the two-step gamma estimates.
#'
#' @param data A data.frame.
#' @param Y.names Character vector of item column names.
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names.
#' @param maxIter.measurement Maximum EM iterations.
#' @param measurement.tol Convergence tolerance.
#' @param covariate.tol NR tolerance for the covariate model.
#' @param iter.measurement Number of random restarts.
#' @param R2.threshold Entropy R\eqn{^2} restart threshold.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}.
#'   Default `FALSE`.
#' @param rebase Character or integer. Reference class for column naming of
#'   `$mGamma`. Must match the `rebase` used in [tseLCA::three_step()] so
#'   coefficient labels are consistent. Default `"C1"`.
#' @param verbose Logical.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{`mGamma`}{Q x (T-1) numeric matrix of multinomial logit
#'       coefficients. Rows are named by covariate (including `"Intercept"`),
#'       columns by non-reference class (e.g. `"C2"`, `"C3"`).}
#'     \item{`mPhi`}{Item parameter matrix (items x classes) from the
#'       fixed-parameter multilevLCA fit.}
#'     \item{`vOmega`}{Length-T vector of marginal class proportions, computed
#'       as the average of the fitted class probability matrix (`vPi_avg` in
#'       multilevLCA output).}
#'     \item{`LLKSeries`}{Matrix of observed-data log-likelihoods across EM
#'       iterations, passed through directly from the multilevLCA fit.}
#'     \item{`raw_fit`}{The full [multilevLCA::multiLCA()] output object,
#'       including `$Varmat_cor` (corrected variance matrix) and
#'       `$SEs_cor_gamma` (corrected standard errors for `mGamma`) if
#'       available.}
#'   }
#' @examples
#' \donttest{
#' d <- generate_data(200, "high", "covariate", seed = 1)
#'
#' # Two-step estimation via multiLCA (fixedpars = 1)
#' fZ_ml <- fitZ_from_multiLCA(
#'   data                = d,
#'   Y.names             = paste0("Y", 1:6),
#'   n_classes           = 3,
#'   Zp.names            = "Zp",
#'   maxIter.measurement = 5000L,
#'   measurement.tol     = 1e-8,
#'   covariate.tol       = 1e-6,
#'   iter.measurement    = 10L,
#'   R2.threshold        = 0.70
#' )
#' fZ_ml$mGamma           # two-step estimates
#' fZ_ml$raw_fit$Varmat_cor   # multilevLCA corrected vcov
#' }
#' @export
fitZ_from_multiLCA <- function(
  data,
  Y.names,
  n_classes,
  Zp.names,
  maxIter.measurement,
  measurement.tol,
  covariate.tol,
  iter.measurement,
  R2.threshold,
  incomplete = FALSE,
  rebase = "C1",
  verbose = FALSE
) {
  run_fit <- function() {
    args <- list(
      data,
      Y.names,
      n_classes,
      Z = Zp.names,
      extout = TRUE,
      incomplete = incomplete,
      maxIter = maxIter.measurement,
      tol = measurement.tol,
      NRtol = covariate.tol,
      fixedpars = 1L,
      verbose = FALSE
    )
    fit <- do.call(multilevLCA::multiLCA, args)
    if (nrow(fit$LLKSeries) == maxIter.measurement) {
      args$maxIter <- 2L * maxIter.measurement
      fit <- do.call(multilevLCA::multiLCA, args)
      if (verbose) {
        warning(sprintf(
          "fitZ multiLCA hit %d iterations; retried with %d.",
          maxIter.measurement,
          2L * maxIter.measurement
        ))
      }
    }
    fit
  }

  initial <- run_fit()
  ll0 <- initial$LLKSeries[nrow(initial$LLKSeries), 1L]

  if (!is.null(initial$R2entr) && initial$R2entr < R2.threshold) {
    if (verbose) {
      warning(sprintf(
        "fitZ multiLCA has low entropy R\u00b2 (%.3f < %.3f). Running %d additional random restarts.",
        initial$R2entr,
        R2.threshold,
        iter.measurement
      ))
    }
    if (iter.measurement > 0L) {
      cands <- lapply(seq_len(iter.measurement), function(r) run_fit())
      cand_lls <- vapply(
        cands,
        function(f) f$LLKSeries[nrow(f$LLKSeries), 1L],
        numeric(1L)
      )
      best_r <- which.max(cand_lls)
      if (cand_lls[best_r] > ll0) {
        if (verbose) {
          message(sprintf(
            "fitZ restart %d improved log-likelihood to %.4f.",
            best_r,
            cand_lls[best_r]
          ))
        }
        initial <- cands[[best_r]]
      } else {
        if (verbose) {
          message(
            "No fitZ restart improved on the initial multiLCA covariate fit."
          )
        }
      }
    }
  }

  raw <- initial
  mGamma <- raw$mGamma
  rownames(mGamma) <- c("Intercept", Zp.names)
  ref_idx <- parse_rebase(rebase, n_classes)
  non_ref_classes <- seq_len(n_classes)[-ref_idx]
  colnames(mGamma) <- paste0("C", non_ref_classes)

  list(
    mGamma = mGamma,
    mPhi = raw$mPhi,
    vOmega = as.vector(raw$vPi_avg),
    LLKSeries = raw$LLKSeries,
    raw_fit = raw
  )
}
