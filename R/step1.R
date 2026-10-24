# tseLCA/R/step1.R
#
# Step 1 (measurement model): internal helpers. The multilevLCA-based
# estimation routines themselves live in R/lca_measurement.R.

#' Posterior class probabilities from a Step-1 fit, in data-row order
#'
#' P(X = t | Y_i) for the rows of `Y.exp` (one-hot expanded, as returned by
#' clean_data()), using the item-response probabilities and class sizes in
#' `fit0` as they stand (i.e. after any rebase permutation).
#' @noRd
step1_posteriors <- function(Y.exp, mDesign, fit0, ivItemcat) {
  log_joint <- sweep(
    log_lik_matrix(Y.exp, expand_Phi(fit0$mPhi, ivItemcat), mDesign),
    2L,
    log(fit0$vPi),
    "+"
  )
  row_max <- apply(log_joint, 1L, max)
  post <- exp(log_joint - (row_max + log(rowSums(exp(log_joint - row_max)))))
  colnames(post) <- paste0("C", seq_len(ncol(post)))
  post
}

#' Step-1 sample used for the measurement-model variance (Sigma.1)
#'
#' Prefers the data stored with the measurement model at fit time
#' (`s1$Y.exp`, `s1$mDesign.exp`, in data-row order). Falls back to decoding
#' multilevLCA's `fit0$mU` for measurement models that do not carry their
#' data (e.g. raw lca_step1() output); `mU` is sorted by response pattern,
#' which is harmless here because its rows are used consistently. Returns
#' NULL if neither is available.
#' @noRd
step1_sample <- function(s1, ivItemcat, ref_idx = 1L) {
  if (!is.null(s1$Y.exp)) {
    return(list(
      Y.exp = s1$Y.exp,
      mDesign = s1$mDesign.exp,
      ivItemcat = ivItemcat,
      u_post = NULL
    ))
  }
  fit0 <- s1$fit0
  if (is.null(fit0$mU)) {
    return(NULL)
  }
  raw <- extract_Y_from_mU(fit0, ivItemcat)
  if (ref_idx != 1L) {
    iT <- length(fit0$vPi)
    raw$u_post <- raw$u_post[, c(ref_idx, seq_len(iT)[-ref_idx]), drop = FALSE]
  }
  raw
}

#' Step 1 for three_step(): the measurement model
#'
#' Uses a pre-fitted measurement model (`step1`: a tseLCA object or raw
#' lca_step1() output), rebased so that class `ref_idx` is the reference, or
#' fits one. When covariates are modeled and `opts$use.two.step` is TRUE,
#' also attaches two-step starting values (`$fitZ`) with the measurement
#' parameters held fixed.
#' @noRd
.fit_step1 <- function(data, Y.names, n_classes, Zp.names, step1, ref_idx, opts,
                       Y.levels = NULL) {
  if (!is.null(step1)) {
    # Normalize: accept raw lca_step1() list or any tseLCA object
    s1 <- if (inherits(step1, "tseLCA")) step1$measurement_model else step1
    # Apply rebase permutation so the desired reference class is column 1.
    s1$fit0 <- permute_fit0_classes(s1$fit0, ref_idx)
    if (!is.null(s1$fitZ)) {
      s1$fitZ <- normalize_fitZ_names(s1$fitZ, n_classes = n_classes)
      s1$fitZ <- permute_fitZ_classes(s1$fitZ, ref_idx)
    }
  } else {
    s1 <- lca_step1(
      data,
      Y.names,
      n_classes,
      Zp.names,
      opts$maxIter.measurement,
      opts$measurement.tol,
      opts$covariate.tol,
      opts$iter.measurement,
      opts$R2.threshold,
      opts$get.twostep.vcov,
      incomplete = opts$incomplete,
      include.intercept = opts$include.intercept,
      rebase = opts$rebase,
      startval = opts$startval,
      n_init = opts$n_init,
      verbose = opts$verbose
    )
  }

  # Two-step coefficients with Step 1 held fixed, when not already computed
  # (e.g. when step1 was passed in from a measurement-only fit).
  if (opts$use.two.step && is.null(s1$fitZ) && !is.null(Zp.names)) {
    s1$fitZ <- fitZ_from_fit0(
      fit0 = s1$fit0,
      data = data,
      Y.names = Y.names,
      Zp.names = Zp.names,
      tol = opts$covariate.tol,
      maxIter = opts$em.maxIter,
      incomplete = opts$incomplete,
      include.intercept = opts$include.intercept,
      rebase = opts$rebase,
      verbose = opts$verbose,
      Y.levels = Y.levels
    )
  }
  s1
}

#' Attach the prepared Step-1 sample and bookkeeping to a measurement model
#'
#' The Step-1 sample is kept in data-row order: it is needed for posteriors
#' and for the Step-1 variance when the model is reused (possibly on another
#' sample) through `step1`. multilevLCA's fit0$mU is not used as the data
#' source because it is sorted by response pattern and its polytomous coding
#' differs between the listwise and FIML paths. A pre-fitted model keeps the
#' sample it was estimated on.
#' @noRd
.attach_step1_data <- function(s1, dat, ref_idx, fitted_here) {
  s1$Y.names <- dat$Y.names
  s1$ivItemcat <- dat$ivItemcat
  s1$Y.levels <- dat$Y.levels
  s1$ref_idx <- ref_idx
  # multilevLCA's N x p score matrix is not used (tseLCA computes its own
  # scores) and would dominate the size of every fitted object.
  if (!is.null(s1$fit0)) s1$fit0$mScore <- NULL
  if (fitted_here) {
    s1$Y.exp <- dat$Y.obs
    s1$mDesign.exp <- dat$mDesign
  }
  s1
}

#' Validate a user-supplied Step-1 variance matrix (log-ratio parameters)
#' @noRd
.check_step1_vcov <- function(V, npar) {
  if (!is.matrix(V) || !is.numeric(V) || !identical(dim(V), c(npar, npar))) {
    stop(sprintf(
      "`vcov` must be a numeric %d x %d matrix (the variance of the Step-1 log-ratio parameters, as in vcov() of a tse_lca() model).",
      npar, npar
    ), call. = FALSE)
  }
  if (anyNA(V) || any(!is.finite(V))) {
    stop("`vcov` must not contain missing or infinite values.", call. = FALSE)
  }
  if (max(abs(V - t(V))) > 1e-6 * max(1, max(abs(V)))) {
    stop("`vcov` must be symmetric.", call. = FALSE)
  }
  if (any(diag(V) < 0)) {
    stop("`vcov` must have non-negative variances on its diagonal.", call. = FALSE)
  }
  V <- (V + t(V)) / 2
  dimnames(V) <- NULL
  V
}

#' Does a measurement model lack the variance needed for corrected SEs?
#'
#' TRUE for a model from as_tse_lca() given without `vcov`: its parameters
#' were not estimated on the stored sample, so Sigma.1 cannot be computed.
#' @noRd
.lacks_step1_vcov <- function(s1) {
  isTRUE(s1$step1_external) && is.null(s1$step1_vcov)
}

.no_step1_vcov_warning <- function() {
  paste0(
    "The measurement model was created by as_tse_lca() without `vcov`, so the ",
    "standard errors cannot be corrected for the Step-1 uncertainty; robust ",
    "standard errors are reported. Supply the variance matrix of the Step-1 ",
    "parameters with as_tse_lca(vcov = )."
  )
}

#' Step-1 variance with the classes reordered so that `ref_idx` comes first
#'
#' `V` is the variance of (log(pi_t / pi_1), t = 2..T; then, class by class,
#' the m free item log-ratios), as in vcov() of a measurement model. After
#' permute_fit0_classes() the parameters are the same functions of the
#' permuted classes, so the variance transforms as A V A'.
#' @noRd
.rebase_step1_vcov <- function(V, ref_idx, iT) {
  if (ref_idx == 1L) {
    return(V)
  }
  m <- (nrow(V) - (iT - 1L)) %/% iT
  ord <- c(ref_idx, seq_len(iT)[-ref_idx])
  A <- matrix(0, nrow(V), ncol(V))
  for (j in 2:iT) {
    if (ord[j] != 1L) A[j - 1L, ord[j] - 1L] <- 1
    if (ref_idx != 1L) A[j - 1L, ref_idx - 1L] <- A[j - 1L, ref_idx - 1L] - 1
  }
  for (j in seq_len(iT)) {
    to <- (iT - 1L) + (j - 1L) * m + seq_len(m)
    from <- (iT - 1L) + (ord[j] - 1L) * m + seq_len(m)
    A[cbind(to, from)] <- 1
  }
  A %*% V %*% t(A)
}

#' Step-1 variance matrix (Sigma.1) of a measurement model
#'
#' Computed on the Step-1 sample stored with the model (step1_sample()),
#' falling back to the current sample. A model with externally given
#' parameters (as_tse_lca()) uses the variance matrix supplied with it.
#' @noRd
.step1_varmat <- function(s1, dat, ref_idx, boundary.tol) {
  if (isTRUE(s1$step1_external)) {
    return(.rebase_step1_vcov(s1$step1_vcov, ref_idx, length(s1$fit0$vPi)))
  }
  sample1 <- step1_sample(s1, dat$ivItemcat, ref_idx)
  if (is.null(sample1)) {
    sample1 <- list(
      Y.exp = dat$Y.obs,
      mDesign = dat$mDesign,
      ivItemcat = dat$ivItemcat,
      u_post = NULL
    )
  }
  lca_indiv_varmat(
    sample1$Y.exp,
    sample1$mDesign,
    s1$fit0,
    sample1$ivItemcat,
    boundary.tol = boundary.tol,
    u_post = sample1$u_post
  )$Varmat
}

#' Indicator categories stored with a measurement model
#'
#' `step1` is a tseLCA object, raw lca_step1() output, or NULL. Returns NULL
#' when no categories are stored (NULL input, or a model fitted before they
#' were recorded), in which case they are derived from the data.
#' @noRd
.measurement_levels <- function(step1) {
  if (is.null(step1)) {
    return(NULL)
  }
  s1 <- if (inherits(step1, "tseLCA")) step1$measurement_model else step1
  s1$Y.levels
}
