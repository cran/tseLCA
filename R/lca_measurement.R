# tseLCA/R/lca_measurement.R
#
# Step-1 measurement model with multilevLCA::multiLCA.
#
# Exports:
#   lca_step1()          - measurement model fit + optional two-step fitZ
#   lca_step1_startval() - measurement model fit from an externally supplied
#                          starting classification or item-response
#                          probability matrix (bypasses multilevLCA's
#                          k-means-on-PCA initialization; see the "External
#                          starting values" section below)
#   fitZ_from_fit0()     - pure-R EM for gamma with mPhi fixed (default fitZ path)
#   fitZ_from_multiLCA() - two-step estimation with multiLCA(fixedpars=1, Z=...) (used when get.twostep.vcov = TRUE in three_step())
#
# External starting values
# -------------------------
# multilevLCA::multiLCA()'s default initialization (k-means on principal
# components) is deterministic given the data and, on some datasets, lands
# on a local rather than global optimum of the Step-1 likelihood.
#
# `startval` (on lca_step1(), lca_step1_startval(), fitZ_from_multiLCA(), and
# three_step()) works around this by accepting either:
#   - an integer vector of length nrow(data), a per-row class assignment
#     (1..n_classes), used as-is; or
#   - a numeric matrix of conditional item-response probabilities
#     P(Y_h = k | X = t) -- one row per (item, category) pair (matching
#     expand_Y()'s column order) and one column per class -- from which a
#     per-row classification is derived (classify_from_phi()) by naive-Bayes
#     argmax under a flat class prior. This is the natural format for an
#     externally estimated Step-1 solution (e.g. poLCA's `probs`, or a
#     published item-response table), which gives item-response profiles
#     rather than a classification for this specific sample.
# Either form is injected into a data column and passed to
# multiLCA(..., startval = <that column>, kmea = FALSE), which skips k-means
# entirely and initializes the EM algorithm from the supplied classification.
#
# `n_init` (on lca_step1(), fitZ_from_multiLCA(), and three_step()) is the
# unconditional multi-random-start analog of StepMix's `n_init` and poLCA's
# `nrep`: it fits the measurement model `n_init` times from independent
# random classifications (not multilevLCA's deterministic k-means-on-PCA
# path) and keeps the highest-log-likelihood fit. This is a separate
# argument from `iter.measurement`, which instead reruns multilevLCA's own
# (kmea = TRUE) k-means initialization, and only when entropy R^2 is low.

# -- External starting values ---------------------------------------------------

#' Derive a Step-1 classification from an item-response probability matrix
#'
#' Computes, for each row of `data`, the class `t` maximizing the naive-Bayes
#' conditional log-likelihood `sum_h log P(Y_ih | X=t)` under a user-supplied
#' item-response probability matrix `phi` (a flat class prior is assumed,
#' since no prevalences are supplied). `phi` must have one row per
#' (item, category) pair -- in the same order as `expand_Y()`'s columns,
#' i.e. items in `Y.names` order, categories `0..K_h-1` within each item --
#' and one column per class. Rows within an item block are renormalized
#' (with a warning-free tolerance of 5%) to sum to 1 per class, to absorb
#' minor rounding in hand-transcribed probability tables.
#' @noRd
classify_from_phi <- function(data, Y.names, n_classes, phi, incomplete = FALSE) {
  cd <- clean_data(data = data, Y.names = Y.names, incomplete = incomplete)
  ivItemcat <- cd$ivItemcat

  if (!is.matrix(phi)) {
    phi <- as.matrix(phi)
  }
  if (nrow(phi) != sum(ivItemcat)) {
    stop(
      sprintf(
        "`startval` matrix must have sum(category counts) = %d rows (one row per item-category pair, in `Y.names` order), got %d.",
        sum(ivItemcat),
        nrow(phi)
      ),
      call. = FALSE
    )
  }
  if (ncol(phi) != n_classes) {
    stop(
      sprintf(
        "`startval` matrix must have n_classes = %d columns, got %d.",
        n_classes,
        ncol(phi)
      ),
      call. = FALSE
    )
  }
  if (anyNA(phi) || any(phi < -1e-8 | phi > 1 + 1e-8)) {
    stop(
      "`startval` matrix must contain probabilities in [0, 1] with no NAs.",
      call. = FALSE
    )
  }

  col_start <- 1L
  for (h in seq_along(ivItemcat)) {
    K_h <- ivItemcat[h]
    idx <- col_start:(col_start + K_h - 1L)
    block <- phi[idx, , drop = FALSE]
    block_sums <- colSums(block)
    if (any(abs(block_sums - 1) > 0.05)) {
      stop(
        sprintf(
          "`startval` matrix rows %d:%d (item '%s') must be conditional probabilities that sum to ~1 within each class column; got column sums ranging [%.3f, %.3f].",
          idx[1L],
          idx[length(idx)],
          Y.names[h],
          min(block_sums),
          max(block_sums)
        ),
        call. = FALSE
      )
    }
    phi[idx, ] <- sweep(block, 2L, block_sums, "/")
    col_start <- col_start + K_h
  }

  phi_clamped <- pmax(pmin(phi, 1 - 1e-10), 1e-10)
  ll <- log_lik_matrix(cd$Y.obs, phi_clamped, cd$mDesign)
  cls_kept <- max.col(ll, ties.method = "first")

  startval <- rep(1L, nrow(data))
  startval[cd$keep_Y] <- cls_kept
  startval
}

#' Validate a startval input and attach it to `data` as a reserved column
#'
#' `multilevLCA::multiLCA()`'s `startval` argument takes the *name* of a
#' `data` column, not a vector, so a user-supplied classification has to be
#' written into `data` before it can be passed through. `startval` may be
#' supplied as an integer classification vector (used as-is) or a numeric
#' matrix of conditional item-response probabilities (converted to a
#' classification with `classify_from_phi()`); either way this helper
#' validates the result (length, range, no missing values) and returns
#' `data` with the classification attached under a reserved column name.
#' @noRd
attach_startval_column <- function(
  data,
  startval,
  n_classes,
  Y.names = NULL,
  incomplete = FALSE
) {
  if (is.matrix(startval) || is.data.frame(startval)) {
    if (is.null(Y.names)) {
      stop(
        "`Y.names` must be supplied when `startval` is an item-response probability matrix.",
        call. = FALSE
      )
    }
    startval <- classify_from_phi(data, Y.names, n_classes, startval, incomplete)
  }

  if (length(startval) != nrow(data)) {
    stop(
      sprintf(
        "`startval` must have length nrow(data) = %d, got %d.",
        nrow(data),
        length(startval)
      ),
      call. = FALSE
    )
  }
  startval_int <- suppressWarnings(as.integer(startval))
  if (anyNA(startval_int)) {
    stop(
      "`startval` must not contain NA (or non-integer-coercible) values.",
      call. = FALSE
    )
  }
  if (any(startval_int < 1L | startval_int > n_classes)) {
    stop(
      sprintf(
        "`startval` must contain integers between 1 and n_classes = %d.",
        n_classes
      ),
      call. = FALSE
    )
  }

  col_name <- ".tseLCA_startval"
  if (col_name %in% names(data)) {
    stop(
      sprintf(
        "`data` already contains a reserved column '%s'. Rename or remove it before supplying `startval`.",
        col_name
      ),
      call. = FALSE
    )
  }
  data[[col_name]] <- startval_int

  list(data = data, col_name = col_name)
}

#' Fit multiLCA from an externally supplied starting classification
#'
#' Injects `startval` into `data` and calls `multilevLCA::multiLCA()` with
#' `startval = <injected column>, kmea = FALSE`, bypassing multilevLCA's
#' default k-means-on-principal-components initialization entirely. Returns
#' the raw (not yet rebase-permuted) `multiLCA` fit object. Because the
#' starting classification is user-supplied, no automatic random restarts
#' are attempted (unlike `run_measurement_fit()`'s low-entropy restart path).
#' @noRd
run_measurement_fit_startval <- function(
  data,
  Y.names,
  n_classes,
  startval,
  maxIter.measurement,
  measurement.tol,
  incomplete,
  verbose
) {
  attached <- attach_startval_column(
    data,
    startval,
    n_classes,
    Y.names = Y.names,
    incomplete = incomplete
  )

  args <- list(
    attached$data,
    Y.names,
    n_classes,
    startval = attached$col_name,
    kmea = FALSE,
    extout = TRUE,
    incomplete = incomplete,
    maxIter = maxIter.measurement,
    tol = measurement.tol,
    verbose = FALSE
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
        "Measurement model still failed to converge even after running more iterations. Consider increasing maxIter.measurement and/or measurement.tol"
      )
    }
  }

  if (verbose && !is.null(fit$R2entr) && fit$R2entr < 0.70) {
    warning(sprintf(
      "Measurement model fit from `startval` has low entropy R\u00b2 (%.3f). Because `startval` is user-supplied, tseLCA does not run automatic random restarts on top of it; inspect the external solution if this is unexpected.",
      fit$R2entr
    ))
  }

  fit
}

#' Fit multiLCA from `n_init` independent random-classification restarts
#'
#' The unconditional multi-random-start analog of StepMix's `n_init` and
#' poLCA's `nrep`: fits the measurement model `n_init` times, each from an
#' independent uniform-random classification (`sample.int(n_classes, N,
#' replace = TRUE)`) injected with `run_measurement_fit_startval()` (so each
#' restart also uses `kmea = FALSE`, not multilevLCA's k-means-on-PCA path),
#' and returns the fit with the highest final log-likelihood. Unlike
#' `run_measurement_fit()`'s `iter.measurement`/`R2.threshold` restarts,
#' this always runs all `n_init` fits regardless of entropy.
#' @noRd
run_measurement_fit_random_restarts <- function(
  data,
  Y.names,
  n_classes,
  n_init,
  maxIter.measurement,
  measurement.tol,
  incomplete,
  verbose
) {
  if (!is.numeric(n_init) || length(n_init) != 1L || n_init < 1L) {
    stop("`n_init` must be a single positive integer.", call. = FALSE)
  }
  n_init <- as.integer(n_init)
  N <- nrow(data)

  fits <- vector("list", n_init)
  lls <- numeric(n_init)
  for (r in seq_len(n_init)) {
    rand_start <- sample.int(n_classes, N, replace = TRUE)
    fits[[r]] <- run_measurement_fit_startval(
      data = data,
      Y.names = Y.names,
      n_classes = n_classes,
      startval = rand_start,
      maxIter.measurement = maxIter.measurement,
      measurement.tol = measurement.tol,
      incomplete = incomplete,
      verbose = FALSE
    )
    lls[r] <- fits[[r]]$LLKSeries[nrow(fits[[r]]$LLKSeries), 1L]
  }

  best <- which.max(lls)
  if (verbose) {
    message(sprintf(
      "Best of %d random-start Step-1 fits: run %d with log-likelihood %.4f (range [%.4f, %.4f]).",
      n_init,
      best,
      lls[best],
      min(lls),
      max(lls)
    ))
  }
  fits[[best]]
}

#' Fit the LCA measurement model from an externally supplied classification
#'
#' A thin wrapper around \pkg{multilevLCA}'s deterministic initialization
#' path. \code{multilevLCA::multiLCA()}'s default Step-1 initialization
#' (k-means on principal components) is deterministic given the data and, on
#' some datasets, converges to a local rather than global optimum of the
#' Step-1 log-likelihood. If you have already found a better solution with
#' an external solver run with many random starts (e.g. \pkg{StepMix},
#' \pkg{poLCA}, or similar), this function lets you inject that
#' classification directly: it writes \code{startval} into a temporary
#' column of \code{data} and calls \code{multiLCA(..., startval = <that
#' column>, kmea = FALSE)}, which skips k-means entirely and initializes the
#' EM algorithm from the supplied classification instead.
#'
#' Most users should not need to call this function directly. Pass
#' \code{startval} to \code{\link{three_step}()} (for structural estimation)
#' or \code{\link{lca_step1}()} (for a measurement-only fit) instead --
#' both implement the same mechanism and return the fitted measurement
#' model as \code{$measurement_model$fit0} / \code{$fit0} respectively. This
#' function is documented mainly to describe what `startval` accepts and
#' how it is used internally.
#'
#' @param data A data.frame containing at minimum the indicator columns.
#' @param Y.names Character vector of item column names.
#' @param n_classes Integer. Number of latent classes.
#' @param startval Either of the following, giving a starting classification
#'   for the measurement model:
#'   \describe{
#'     \item{An integer vector}{Length \code{nrow(data)}, a starting class
#'       assignment (\code{1..n_classes}) for every row of \code{data},
#'       typically obtained from an external latent class solver run with
#'       many random starts (e.g. the modal class from many-random-start
#'       posterior probabilities, as in \pkg{StepMix} or \pkg{poLCA}).}
#'     \item{A numeric matrix}{A conditional item-response probability
#'       matrix \eqn{P(Y_h = k \mid X = t)} with one row per (item, category)
#'       pair -- items in \code{Y.names} order, categories \code{0..K_h-1}
#'       within each item, matching the column order of
#'       \code{expand_Y(data[, Y.names], ivItemcat)} -- and one column per
#'       class. A per-row classification is derived internally by
#'       naive-Bayes argmax under a flat class prior (see
#'       \code{classify_from_phi()}). This is the natural format for an
#'       externally estimated Step-1 solution that isn't tied to this
#'       specific sample, e.g. \pkg{poLCA}'s \code{probs} output or a
#'       published item-response table.}
#'   }
#'   No automatic random restarts are performed on top of this starting
#'   value (contrast \code{\link{lca_step1}()}'s
#'   \code{iter.measurement}/\code{R2.threshold} restart logic, which applies
#'   only to multilevLCA's own k-means initialization, and its \code{n_init}
#'   argument, which does run unconditional random restarts but from
#'   independent random classifications rather than a single fixed one).
#' @param maxIter.measurement Maximum EM iterations before giving up on
#'   convergence. Default `5000L`.
#' @param measurement.tol Convergence tolerance. Default `1e-8`.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}. Default `FALSE`.
#' @param rebase Character or integer specifying the reference latent class.
#'   Use `"C1"`, `"C2"`, etc. or an integer index. Default `"C1"`. The
#'   measurement model is permuted so this class becomes column 1, making it
#'   the reference for all downstream multinomial logit parameterizations.
#' @param verbose Logical. Print progress messages. Default `FALSE`.
#'
#' @return A list with `$fit0` ([multilevLCA::multiLCA()] measurement model,
#'   rebase-permuted per \code{rebase}) and \code{$fitZ = NULL}. This is the
#'   same shape as \code{\link{lca_step1}()}'s return value, and matches
#'   \code{\link{three_step}()}'s \code{$measurement_model} when
#'   \code{startval} is passed there directly.
#' @examples
#' \donttest{
#' d <- generate_data(200, "high", "covariate", seed = 1)
#'
#' # Recommended: pass `startval` to three_step() (or lca_step1() for a
#' # measurement-only fit) rather than calling this function directly --
#' # both use this same mechanism internally.
#'
#' # A starting classification from an external solver (here, the DGP's own
#' # true classes, standing in for e.g. a StepMix solution with many
#' # random starts):
#' fit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                   startval = d$X)
#' fit$measurement_model$fit0$vPi
#'
#' # Equivalently, supply a conditional item-response probability matrix
#' # (one row per item-category pair, in Y.names order -- since all 6 items
#' # here are binary, each contributes 2 rows: P(Y=0|C), P(Y=1|C)). In
#' # practice this would come from an external solver (e.g. poLCA's `probs`
#' # or a published item-response table); here a quick first-pass fit
#' # stands in for that external source.
#' fit_ref <- three_step(d, paste0("Y", 1:6), n_classes = 3)$measurement_model$fit0
#' phi <- matrix(0, nrow = 12, ncol = 3)
#' for (h in 1:6) {
#'   phi[2 * h - 1, ] <- 1 - fit_ref$mPhi[h, ] # P(Y_h = 0 | C)
#'   phi[2 * h,     ] <- fit_ref$mPhi[h, ]     # P(Y_h = 1 | C)
#' }
#' fit_phi <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       startval = phi)
#' }
#' @keywords internal
#' @export
lca_step1_startval <- function(
  data,
  Y.names,
  n_classes,
  startval,
  maxIter.measurement = 5000L,
  measurement.tol = 1e-8,
  incomplete = FALSE,
  rebase = "C1",
  verbose = FALSE
) {
  fit0 <- run_measurement_fit_startval(
    data = data,
    Y.names = Y.names,
    n_classes = n_classes,
    startval = startval,
    maxIter.measurement = maxIter.measurement,
    measurement.tol = measurement.tol,
    incomplete = incomplete,
    verbose = verbose
  )

  ref_idx <- parse_rebase(rebase, n_classes)
  fit0 <- permute_fit0_classes(fit0, ref_idx)

  list(fit0 = fit0, fitZ = NULL)
}


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
#' @param startval Optional starting classification for the Step-1
#'   measurement model: either an integer vector of length `nrow(data)`
#'   (`1..n_classes` per row) or a numeric matrix of conditional
#'   item-response probabilities from which a classification is derived. See
#'   [lca_step1_startval()] for the full description of both forms. When
#'   supplied, `lca_step1()` fits the measurement model with
#'   [lca_step1_startval()] instead of multilevLCA's default
#'   k-means-on-principal-components initialization, and `estimate.one.step`,
#'   `iter.measurement`, and `R2.threshold` (which govern the default
#'   restart-on-low-entropy behavior) are ignored. Mutually exclusive with
#'   `n_init`. Default `NULL`.
#' @param n_init Optional positive integer. If supplied, fits the
#'   measurement model `n_init` times from independent uniform-random
#'   classifications (each through `startval`-style injection with
#'   `kmea = FALSE`, not multilevLCA's k-means-on-PCA path) and keeps the
#'   fit with the highest log-likelihood -- the unconditional multi-start
#'   analog of `n_init` in \pkg{StepMix} or `nrep` in \pkg{poLCA}. Unlike
#'   `iter.measurement` (which reruns multilevLCA's own k-means
#'   initialization, and only when entropy R\eqn{^2} is low), all `n_init`
#'   fits are always run. `estimate.one.step`, `iter.measurement`, and
#'   `R2.threshold` are ignored when `n_init` is supplied. Mutually exclusive
#'   with `startval`. Default `NULL`.
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
#'
#' # Many random-classification restarts, keeping the best (analogous to
#' # n_init in StepMix or nrep in poLCA)
#' s1r <- lca_step1(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                  n_init = 20L, verbose = TRUE)
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
  startval = NULL,
  n_init = NULL,
  verbose = FALSE
) {
  if (!is.null(startval) && !is.null(n_init)) {
    stop(
      "`startval` and `n_init` are mutually exclusive ways of controlling ",
      "Step-1 initialization: supply a fixed starting classification with ",
      "`startval`, or a number of independent random restarts with ",
      "`n_init`, not both.",
      call. = FALSE
    )
  }

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

  fit0 <- if (!is.null(startval)) {
    run_measurement_fit_startval(
      data = data,
      Y.names = Y.names,
      n_classes = n_classes,
      startval = startval,
      maxIter.measurement = maxIter.measurement,
      measurement.tol = measurement.tol,
      incomplete = incomplete,
      verbose = verbose
    )
  } else if (!is.null(n_init)) {
    run_measurement_fit_random_restarts(
      data = data,
      Y.names = Y.names,
      n_classes = n_classes,
      n_init = n_init,
      maxIter.measurement = maxIter.measurement,
      measurement.tol = measurement.tol,
      incomplete = incomplete,
      verbose = verbose
    )
  } else if (estimate.one.step) {
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
#' `mGamma` (Q x (T-1)) with an EM algorithm using a BFGS M-step.
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

#' Estimate two-step covariate model with multilevLCA (optional reference path)
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
#' @param startval Optional starting classification for the measurement
#'   portion of this `multiLCA(fixedpars = 1)` fit -- an integer vector or a
#'   conditional item-response probability matrix, as described in
#'   [lca_step1_startval()] -- e.g. the same value passed to [lca_step1()]
#'   for the primary Step-1 fit. When supplied, `kmea = FALSE` is used and
#'   `iter.measurement`/`R2.threshold` restarts are skipped, for the same
#'   reasons as in [lca_step1_startval()]. Mutually exclusive with `n_init`.
#'   Default `NULL`.
#' @param n_init Optional positive integer. If supplied, fits this
#'   `multiLCA(fixedpars = 1)` model `n_init` times from independent
#'   uniform-random classifications (`kmea = FALSE`) and keeps the fit with
#'   the highest log-likelihood, as in [lca_step1()]'s `n_init` argument.
#'   `iter.measurement`/`R2.threshold` restarts are skipped. Mutually
#'   exclusive with `startval`. Default `NULL`.
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
#' # Two-step estimation with multiLCA (fixedpars = 1)
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
  startval = NULL,
  n_init = NULL,
  verbose = FALSE
) {
  if (!is.null(startval) && !is.null(n_init)) {
    stop(
      "`startval` and `n_init` are mutually exclusive ways of controlling ",
      "this fit's initialization: supply a fixed starting classification ",
      "with `startval`, or a number of independent random restarts with ",
      "`n_init`, not both.",
      call. = FALSE
    )
  }

  startval_col <- NULL
  if (!is.null(startval)) {
    attached <- attach_startval_column(
      data,
      startval,
      n_classes,
      Y.names = Y.names,
      incomplete = incomplete
    )
    data <- attached$data
    startval_col <- attached$col_name
  }

  run_fit <- function(fit_data = data, fit_startval_col = startval_col) {
    args <- list(
      fit_data,
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
    if (!is.null(fit_startval_col)) {
      args$startval <- fit_startval_col
      args$kmea <- FALSE
    }
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

  if (!is.null(n_init)) {
    if (!is.numeric(n_init) || length(n_init) != 1L || n_init < 1L) {
      stop("`n_init` must be a single positive integer.", call. = FALSE)
    }
    n_init_int <- as.integer(n_init)
    N <- nrow(data)
    fits <- vector("list", n_init_int)
    lls <- numeric(n_init_int)
    for (r in seq_len(n_init_int)) {
      rand_start <- sample.int(n_classes, N, replace = TRUE)
      attached_r <- attach_startval_column(
        data,
        rand_start,
        n_classes,
        Y.names = Y.names,
        incomplete = incomplete
      )
      fits[[r]] <- run_fit(
        fit_data = attached_r$data,
        fit_startval_col = attached_r$col_name
      )
      lls[r] <- fits[[r]]$LLKSeries[nrow(fits[[r]]$LLKSeries), 1L]
    }
    best <- which.max(lls)
    if (verbose) {
      message(sprintf(
        "Best of %d random-start fitZ multiLCA fits: run %d with log-likelihood %.4f (range [%.4f, %.4f]).",
        n_init_int,
        best,
        lls[best],
        min(lls),
        max(lls)
      ))
    }
    initial <- fits[[best]]
  } else {
    initial <- run_fit()
    ll0 <- initial$LLKSeries[nrow(initial$LLKSeries), 1L]

    if (!is.null(startval_col)) {
      if (
        verbose && !is.null(initial$R2entr) && initial$R2entr < R2.threshold
      ) {
        warning(sprintf(
          "fitZ multiLCA fit from `startval` has low entropy R\u00b2 (%.3f). No automatic restarts performed because `startval` was user-supplied.",
          initial$R2entr
        ))
      }
    } else if (!is.null(initial$R2entr) && initial$R2entr < R2.threshold) {
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
