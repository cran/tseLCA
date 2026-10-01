# tseLCA/R/tse_structural.R
#
# Step 3 building blocks: covariate models, distal outcome models, and the
# two-step estimator.

#' Relate latent classes to covariates (Step 3)
#'
#' Estimates a multinomial logistic regression of latent class membership on
#' covariates, \eqn{P(X = t \mid Z) \propto \exp(Z\gamma_t)}, correcting for
#' the classification error of the Step-2 class assignments. The measurement
#' model is held fixed, so the covariates cannot change the classes.
#'
#' @details
#' Estimators (`method`):
#' * `"ML"` (default): the bias-adjusted maximum likelihood estimator of
#'   Vermunt (2010), which treats the assigned class as an indicator of the
#'   true class with known classification-error probabilities.
#' * `"BCH"`: the Bolck-Croon-Hagenaars estimator (Bolck, Croon, and
#'   Hagenaars 2004; Vermunt 2010), which reweights the assignments with the
#'   inverse of the classification-error matrix.
#' * `"none"`: the uncorrected three-step estimator (a weighted multinomial
#'   logit of the assigned classes), which is biased toward zero; provided for
#'   comparison.
#'
#' Standard errors (`se`): `"corrected"` (default) adds the uncertainty of
#' the Step-1 measurement model (Bakk, Oberski, and Vermunt 2014) to the
#' robust (sandwich) Step-3 variance; `"robust"` omits it. For `"BCH"` the
#' robust variance is used, which accounts for the Step-1 uncertainty through
#' the weights (Vermunt 2010); for `"none"` the robust variance is used.
#'
#' @param object A classification from [tse_classify()].
#' @param formula One-sided formula for the covariates, e.g. `~ age + sex`.
#'   Factors, interactions, and transformations are allowed.
#' @param method Step-3 estimator: `"ML"`, `"BCH"`, or `"none"` (see Details).
#' @param se Standard errors: `"corrected"` or `"robust"` (see Details).
#' @param ref Reference class of the multinomial logit: a class number or
#'   label such as `"C2"`.
#' @param start Optional starting values: a (Q+1) x (T-1) coefficient matrix
#'   (Q covariates plus the intercept, T classes). By default, the two-step estimates
#'   (Bakk and Kuha 2018) are used.
#' @param control Estimation settings; default: those of the measurement
#'   model. See [tse_control()].
#' @param data Optional data frame with the covariates: the classified data
#'   (the same rows, in the same order) with any additional columns. Omitted:
#'   the data stored in `object`.
#'
#' @return A `tseLCA_covariate` object; see [coef.tseLCA_structural()],
#'   [summary.tseLCA_structural()], [predict.tseLCA_covariate()], and
#'   [anova.tseLCA_covariate()]. Pass it to [tse_distal()] to also model a
#'   distal outcome.
#'
#' @references
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#'   classes and external variables. \emph{Psychometrika}, 83(4), 871--892.
#'   \doi{10.1007/s11336-017-9592-7}
#'
#' Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
#'   assignments to external variables: Standard errors for correct
#'   inference. \emph{Political Analysis}, 22(4), 520--540.
#'   \doi{10.1093/pan/mpu003}
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#'   models with categorical variables: One-step versus three-step
#'   estimators. \emph{Political Analysis}, 12(1), 3--27.
#'   \doi{10.1093/pan/mph001}
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#'   three-step approaches. \emph{Political Analysis}, 18(4), 450--469.
#'   \doi{10.1093/pan/mpq025}
#'
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' cl <- tse_classify(m)
#' fit <- tse_covariate(cl, ~ Zp)
#' summary(fit)
#' confint(fit)
#' predict(fit, newdata = data.frame(Zp = 1:5))
#'
#' # BCH, and the uncorrected estimator for comparison
#' coef(tse_covariate(cl, ~ Zp, method = "BCH"))
#' coef(tse_covariate(cl, ~ Zp, method = "none"))
#' @export
tse_covariate <- function(
  object,
  formula,
  method = c("ML", "BCH", "none"),
  se = c("corrected", "robust"),
  ref = 1,
  start = NULL,
  control = NULL,
  data = NULL
) {
  cl <- match.call()
  if (!inherits(object, "tseLCA_classify")) {
    stop("`object` must be a classification from tse_classify().", call. = FALSE)
  }
  method <- match.arg(method)
  se <- match.arg(se)
  formula <- .covariate_formula(formula, NULL, TRUE)
  object <- .with_data(object, data)
  setup <- .structural_setup(object, formula, NULL, "gaussian", ref, method, se, control)
  iT <- object$n_classes
  opts <- setup$opts
  s1 <- setup$s1
  Q <- ncol(setup$dat$Z_mat)

  if (!is.null(start)) {
    if (!is.matrix(start) || !identical(dim(start), c(Q, iT - 1L))) {
      stop(sprintf(
        "`start` must be a %d x %d matrix (design columns x non-reference classes).",
        Q, iT - 1L
      ), call. = FALSE)
    }
    opts$gamma_start <- start
    opts$use.two.step <- FALSE
  } else if (method != "none") {
    s1$fitZ <- .two_step_start(setup, formula, ref)
  } else {
    opts$use.two.step <- FALSE
  }

  res <- .fit_covariate(setup$dat, s1, setup$s2, setup$Sigma.1, iT, opts)
  .finish_structural(res$fit, cl, setup, method, se, object, formula = formula)
}

#' Relate latent classes to a distal outcome (Step 3)
#'
#' Estimates the class-specific distribution of a distal outcome, correcting
#' for the classification error of the Step-2 class assignments (Bakk,
#' Tekle, and Vermunt 2013). Given a covariate model from [tse_covariate()],
#' the class prior depends on the covariates and the covariate-model
#' uncertainty is propagated to the distal estimates.
#'
#' @details
#' The class-specific parameters are means (`gaussian`, with a common
#' within-class variance, reported as `$sigma2`), log means (`poisson`),
#' logits (`binomial`), or category probabilities (`"multinomial"`). The
#' estimators and standard errors are as for [tse_covariate()]. Use
#' [omnibus_test()] to test whether the outcome differs across classes.
#'
#' @param object A classification from [tse_classify()], or a covariate
#'   model from [tse_covariate()] (combined model).
#' @param formula `Zo ~ 1`, with the distal outcome on the left-hand side.
#' @param family Distribution of the outcome within classes: `"gaussian"`
#'   (default), `"poisson"`, `"binomial"`, `"multinomial"` (nominal
#'   outcome), or the corresponding family object (`gaussian()`,
#'   `poisson()`, `binomial()`; canonical links only).
#' @param method,se As for [tse_covariate()]. For a combined model they
#'   default to those of the covariate model.
#' @param control Estimation settings; default: those of `object`.
#' @param data Optional data frame with the distal outcome, as in
#'   [tse_covariate()].
#'
#' @return A `tseLCA_distal` object, or a `tseLCA_both` object when `object`
#'   is a covariate model.
#'
#' @references
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
#'   association between latent class membership and external variables
#'   using bias-adjusted three-step approaches. \emph{Sociological
#'   Methodology}, 43(1), 272--311. \doi{10.1177/0081175012470644}
#'
#' @examples
#' d <- generate_data(500, "high", "distal", seed = 2)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' fd <- tse_distal(tse_classify(m, assignment = "proportional"), Zo ~ 1)
#' summary(fd)
#' omnibus_test(fd)
#' @export
tse_distal <- function(
  object,
  formula,
  family = "gaussian",
  method = NULL,
  se = NULL,
  control = NULL,
  data = NULL
) {
  cl <- match.call()
  family <- .distal_family(family)
  Zo.name <- .distal_outcome_from_formula(formula)

  if (inherits(object, "tseLCA_covariate")) {
    classify <- object[["classification", exact = TRUE]]
    if (is.null(classify)) {
      stop("The covariate model must come from tse_covariate().", call. = FALSE)
    }
    classify <- .with_data(classify, data)
    method <- .match_or_inherit(method, object$method, c("ML", "BCH", "none"), "method")
    se <- .match_or_inherit(se, object$se_requested, c("corrected", "robust"), "se")
    if (method != object$method) {
      stop("A combined model uses one estimator: `method` must match the covariate model.",
           call. = FALSE)
    }
    setup <- .structural_setup(classify, object$formula, Zo.name, family, object$ref,
                               method, se, control)
    Z_mat <- setup$dat$Z_mat
    cov <- list(
      par = as.vector(object$three_step),
      Sigma.3 = object$three_step_vcov,
      p.xz = function(params) {
        eta <- cbind(0, Z_mat %*% params)
        e <- exp(eta - apply(eta, 1L, max))
        e / rowSums(e)
      }
    )
    cov_fit <- object
  } else if (inherits(object, "tseLCA_classify")) {
    classify <- .with_data(object, data)
    method <- .match_or_inherit(method, "ML", c("ML", "BCH", "none"), "method")
    se <- .match_or_inherit(se, "corrected", c("corrected", "robust"), "se")
    setup <- .structural_setup(classify, NULL, Zo.name, family, 1, method, se, control)
    cov <- NULL
    cov_fit <- NULL
  } else {
    stop(
      "`object` must be a classification from tse_classify() or a covariate ",
      "model from tse_covariate().",
      call. = FALSE
    )
  }

  dis <- .fit_distal(setup$dat, setup$s1, setup$s2, setup$Sigma.1, cov,
                     classify$n_classes, family, setup$opts)
  dis <- .distal_original_order(dis, setup$ref_idx, classify$n_classes)
  fit <- .new_structural_fit(setup$s1, setup$s2, cov_fit, dis, classify$n_classes,
                             family, setup$opts$use.bch,
                             estimator = .estimator_label(setup$opts))
  fit <- .finish_structural(fit, cl, setup, method, se, classify)
  fit$outcome <- Zo.name
  fit
}

#' Two-step estimates of covariate effects
#'
#' Estimates the multinomial logit of class membership on covariates with
#' the measurement-model parameters held fixed at their Step-1 values (Bakk
#' and Kuha 2018). Unlike the three-step estimators, the indicators enter the
#' Step-2 likelihood directly, so no classification step is needed.
#'
#' @param object A measurement model from [tse_lca()] (it must keep its
#'   data).
#' @param formula One-sided covariate formula.
#' @param ref Reference class of the multinomial logit.
#' @param se Logical. If `TRUE`, the estimates and their standard errors
#'   (corrected for the Step-1 uncertainty) are obtained with the two-step
#'   estimator of \pkg{multilevLCA}, initialized at this model's classes; its
#'   measurement model is checked against `object`. If `FALSE` (default),
#'   only the estimates are computed, and the variance is `NA`.
#' @param control Estimation settings; default: those of `object`.
#'
#' @return A `tseLCA_twostep` object (also a `tseLCA_covariate`).
#'
#' @references
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#'   classes and external variables. \emph{Psychometrika}, 83(4), 871--892.
#'   \doi{10.1007/s11336-017-9592-7}
#'
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' coef(tse_twostep(m, ~ Zp))
#' @export
tse_twostep <- function(object, formula, ref = 1, se = FALSE, control = NULL) {
  cl <- match.call()
  if (!inherits(object, "tseLCA_measurement") || object$n_classes < 2L) {
    stop("`object` must be a measurement model from tse_lca() with two or more classes.",
         call. = FALSE)
  }
  if (is.null(object$data)) {
    stop("This measurement model does not store its data; refit it with tse_lca().",
         call. = FALSE)
  }
  formula <- .covariate_formula(formula, NULL, TRUE)
  classify <- tse_classify(object, control = control)
  setup <- .structural_setup(classify, formula, NULL, "gaussian", ref, "ML", "robust", control)
  iT <- object$n_classes
  s1 <- setup$s1
  dat <- setup$dat
  fz <- .two_step_start(setup, formula, ref)
  est <- fz$mGamma
  V <- matrix(NA_real_, length(est), length(est))

  if (se) {
    ml <- .two_step_multilevLCA(object, setup, formula, ref)
    est <- ml$mGamma
    V <- ml$vcov
  }
  nms <- as.vector(outer(rownames(est), colnames(est), paste, sep = ":"))
  dimnames(V) <- list(nms, nms)

  rows <- dat$keep_step3_Z_in_Y
  Y_cc <- dat$Y.obs[rows, , drop = FALSE]
  mDes_cc <- if (!is.null(dat$mDesign)) dat$mDesign[rows, , drop = FALSE] else NULL
  llik <- joint_log_lik(Y_cc, dat$Z_mat, expand_Phi(s1$fit0$mPhi, dat$ivItemcat), est, mDes_cc)
  k <- iT * sum(dat$ivItemcat - 1L) + ncol(dat$Z_mat) * (iT - 1L)
  fit <- structure(
    list(
      measurement_model = s1,
      three_step = est,
      three_step_vcov = V,
      llik = llik,
      AIC = -2 * llik + 2 * k,
      BIC = -2 * llik + k * log(nrow(Y_cc)),
      npar = k,
      nobs = nrow(Y_cc),
      n_classes = iT,
      estimator = "two-step",
      posteriors = classify$posteriors,
      classifications = classify$classifications
    ),
    class = c("tseLCA_twostep", "tseLCA_covariate", "tseLCA_structural", "tseLCA")
  )
  fit <- .finish_structural(fit, cl, setup, "two-step", if (se) "corrected" else "none",
                            classify, formula = formula)
  fit
}

# -- internals ---------------------------------------------------------------------

#' Use a data frame with the classified rows (and possibly more columns)
#'
#' `data` must hold the rows of the classification's data, in the same order,
#' with the same indicator values; other columns may be added.
#' @noRd
.with_data <- function(classify, data) {
  if (is.null(data)) {
    return(classify)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  old <- classify$data
  items <- classify$measurement_model$Y.names
  same <- nrow(data) == nrow(old) && all(items %in% names(data)) &&
    isTRUE(all.equal(
      lapply(data[items], as.character),
      lapply(old[items], as.character),
      check.attributes = FALSE
    ))
  if (!same) {
    stop(
      "`data` must contain the classified data (same rows, in the same order, ",
      "with the same indicator values), possibly with additional columns. To ",
      "classify other data, use tse_classify(newdata = ).",
      call. = FALSE
    )
  }
  classify$data <- data
  classify
}

#' Common set-up of a Step-3 model on a classification
#'
#' Prepares the data with the structural variables, rebases the measurement
#' model to the reference class (recomputing Step 2 if needed), computes the
#' Step-1 variance when corrected standard errors are requested, and, for the
#' uncorrected estimator, replaces the classification-error matrix by the
#' identity.
#' @noRd
.structural_setup <- function(classify, Zp.formula, Zo.name, family, ref, method, se, control) {
  s1 <- classify$measurement_model
  iT <- classify$n_classes
  if (is.null(control)) control <- classify$control
  if (is.null(control)) control <- tse_control()
  ref_idx <- parse_rebase(ref, iT)
  opts <- .opts_from_control(
    control,
    incomplete = classify$missing == "fiml",
    include.intercept = TRUE,
    use.modal.assignment = classify$assignment == "modal",
    use.bch = method == "BCH",
    use.simple.cov = se == "robust" || method %in% c("BCH", "none"),
    uncorrected = method == "none",
    use.two.step = TRUE,
    get.twostep.vcov = FALSE,
    rebase = ref,
    startval = NULL
  )

  rec <- .recode_indicators(classify$data, s1$Y.names, s1$Y.levels)
  dat <- .prepare_data(rec$data, s1$Y.names, NULL, Zo.name, family, opts, s1$Y.levels,
                       Zp.formula = Zp.formula)
  if (!identical(as.integer(dat$keep_Y), as.integer(classify$rows))) {
    stop("Internal error: Step-2 rows do not match the classification.", call. = FALSE)
  }

  if (ref_idx != 1L) {
    s1$fit0 <- permute_fit0_classes(s1$fit0, ref_idx)
    s1$ref_idx <- ref_idx
    s2 <- .step2(dat, s1$fit0, iT, opts)
  } else {
    s2 <- list(
      all = classify$step2,
      cov = if (!is.null(dat$Z_mat)) .subset_step2(classify$step2, dat$keep_step3_Z_in_Y, dat, iT),
      dis = if (!is.null(dat$Zo_mat)) .subset_step2(classify$step2, dat$keep_step3_Zo_in_Y, dat, iT)
    )
  }
  if (method == "none") {
    if (!is.null(s2$cov)) s2$cov$p.wx_mat <- diag(iT)
    if (!is.null(s2$dis)) s2$dis$p.wx_mat <- diag(iT)
  }
  Sigma.1 <- if (opts$use.simple.cov || opts$use.bch) {
    NULL
  } else {
    .step1_varmat(s1, dat, ref_idx, opts$boundary.tol)
  }
  se_fallback <- !is.null(Sigma.1) && anyNA(Sigma.1)
  if (se_fallback) {
    warning(
      "The Step-1 information matrix is singular, so the standard errors cannot be ",
      "corrected for the Step-1 uncertainty; robust standard errors are reported.",
      call. = FALSE
    )
    Sigma.1 <- NULL
    opts$use.simple.cov <- TRUE
  }
  list(s1 = s1, s2 = s2, dat = dat, opts = opts, Sigma.1 = Sigma.1, ref_idx = ref_idx,
       ref = ref, se_fallback = se_fallback)
}

#' Two-step estimates (measurement model fixed) used as Step-3 starting values
#' @noRd
.two_step_start <- function(setup, formula, ref) {
  opts <- setup$opts
  fitZ_from_fit0(
    fit0 = setup$s1$fit0,
    data = setup$dat$data,
    Y.names = setup$s1$Y.names,
    Zp.names = NULL,
    tol = opts$covariate.tol,
    maxIter = opts$em.maxIter,
    incomplete = opts$incomplete,
    rebase = ref,
    verbose = opts$verbose,
    Y.levels = setup$s1$Y.levels,
    Zp.formula = formula
  )
}

#' multilevLCA two-step estimates and corrected variance, initialized at the
#' classes of `object` and checked against its measurement model
#' @noRd
.two_step_multilevLCA <- function(object, setup, formula, ref) {
  dat <- setup$dat
  opts <- setup$opts
  # covariate design columns (without intercept) as plain data columns
  Z <- dat$Z_mat
  has_int <- "(Intercept)" %in% colnames(Z)
  Zx <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  if (!has_int || ncol(Zx) == 0L) {
    stop("`se = TRUE` needs a covariate formula with an intercept and at least one covariate.",
         call. = FALSE)
  }
  zn <- paste0(".tseLCA_z", seq_len(ncol(Zx)))
  rows <- dat$keep_Y[dat$keep_step3_Z_in_Y]
  d2 <- dat$data[rows, setup$s1$Y.names, drop = FALSE]
  d2[zn] <- as.data.frame(Zx)
  start <- max.col(step1_posteriors(
    dat$Y.obs[dat$keep_step3_Z_in_Y, , drop = FALSE],
    if (!is.null(dat$mDesign)) dat$mDesign[dat$keep_step3_Z_in_Y, , drop = FALSE] else NULL,
    object$measurement_model$fit0,
    dat$ivItemcat
  ))
  fz <- fitZ_from_multiLCA(
    data = d2,
    Y.names = setup$s1$Y.names,
    n_classes = object$n_classes,
    Zp.names = zn,
    maxIter.measurement = opts$maxIter.measurement,
    measurement.tol = opts$measurement.tol,
    covariate.tol = opts$covariate.tol,
    iter.measurement = opts$iter.measurement,
    R2.threshold = opts$R2.threshold,
    incomplete = opts$incomplete,
    startval = start,
    verbose = opts$verbose
  )
  # multilevLCA starts from `object`'s classes in their original order and
  # uses class 1 as the reference; compare with that order, then move the
  # estimates to `ref`.
  phi_ml <- fz$raw_fit$mPhi
  if (!is.null(phi_ml) && !isTRUE(all.equal(unname(as.matrix(phi_ml)),
                                            unname(as.matrix(object$measurement_model$fit0$mPhi)),
                                            tolerance = 1e-3))) {
    warning(
      "multilevLCA's two-step measurement model differs from `object`; ",
      "the standard errors refer to its measurement model.",
      call. = FALSE
    )
  }
  est <- fz$mGamma
  V <- fz$Varmat_cor
  if (is.null(V)) V <- fz$raw_fit$Varmat_cor
  if (is.null(V)) V <- matrix(NA_real_, length(est), length(est))
  moved <- .rebase_logit(est, V, setup$ref_idx)
  rownames(moved$mGamma) <- colnames(Z)
  moved
}

#' Change the reference class of multinomial logit coefficients
#'
#' `est` is (Q+1) x (T-1) with class 1 as the reference and `V` the variance of
#' `vec(est)`. Returns the coefficients against class `ref` (columns: the
#' other classes in their original order) and the transformed variance.
#' @noRd
.rebase_logit <- function(est, V, ref) {
  iT <- ncol(est) + 1L
  Q <- nrow(est)
  M <- matrix(0, iT - 1L, iT - 1L) # new non-reference classes x old classes 2..iT
  for (i in seq_len(iT - 1L)) {
    t <- seq_len(iT)[-ref][i]
    if (t != 1L) M[i, t - 1L] <- 1
    if (ref != 1L) M[i, ref - 1L] <- -1
  }
  A <- kronecker(M, diag(Q))
  new <- matrix(A %*% as.vector(est), Q, iT - 1L,
                dimnames = list(rownames(est), paste0("C", seq_len(iT)[-ref])))
  list(mGamma = new, vcov = A %*% V %*% t(A))
}

#' Distal parameters in the measurement model's class order
#'
#' A combined model is estimated with the classes rebased to the covariate
#' model's reference class (reference first, then the others in order). Put
#' the class-specific distal parameters, and their variance, back in the
#' original class order, so that `C<t>` means class t as everywhere else.
#' @noRd
.distal_original_order <- function(dis, ref_idx, iT) {
  if (ref_idx == 1L) {
    return(dis)
  }
  back <- match(seq_len(iT), c(ref_idx, seq_len(iT)[-ref_idx]))
  if (is.matrix(dis$three_step)) { # multinomial: classes x categories
    idx <- as.vector(outer(back, (seq_len(ncol(dis$three_step)) - 1L) * iT, "+"))
    dis$three_step[] <- dis$three_step[back, , drop = FALSE]
  } else {
    idx <- back
    dis$three_step[] <- dis$three_step[back]
  }
  dis$three_step_vcov[] <- dis$three_step_vcov[idx, idx, drop = FALSE]
  dis
}

#' Record the specification on a fitted Step-3 object
#' @noRd
.finish_structural <- function(fit, cl, setup, method, se, classify, formula = NULL) {
  fit$call <- cl
  if (!is.null(formula)) {
    fit$formula <- formula
    fit$terms <- setup$dat$Z_terms
    fit$xlevels <- setup$dat$Z_xlevels
  }
  fit$method <- method
  fit$se_requested <- se
  fit$se <- if (method %in% c("BCH", "none") || isTRUE(setup$se_fallback)) "robust" else se
  fit$assignment <- classify$assignment
  fit$ref <- paste0("C", setup$ref_idx)
  fit$classification <- classify
  # The model was estimated with the classes rebased to `ref`; report the
  # measurement model and the Step-2 posteriors in the original class order,
  # so that class t is the same class in every accessor.
  fit$measurement_model <- classify$measurement_model
  fit$posteriors <- classify$posteriors
  fit$classifications <- classify$classifications
  fit
}

#' Validate a distal outcome family (string or family object)
#' @noRd
.distal_family <- function(family) {
  families <- c("gaussian", "poisson", "binomial", "multinomial")
  if (inherits(family, "family")) {
    canonical <- c(gaussian = "identity", poisson = "log", binomial = "logit")
    if (!family$family %in% names(canonical) || family$link != canonical[[family$family]]) {
      stop(
        "Supported family objects: gaussian(\"identity\"), poisson(\"log\"), binomial(\"logit\").",
        call. = FALSE
      )
    }
    return(family$family)
  }
  if (!is.character(family) || length(family) != 1L || !family %in% families) {
    stop("`family` must be one of ", paste(dQuote(families, FALSE), collapse = ", "),
         " or a family object.", call. = FALSE)
  }
  family
}

#' Distal outcome name from `Zo ~ 1`
#' @noRd
.distal_outcome_from_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.name(formula[[2L]])) {
    stop("`formula` must be `outcome ~ 1`, with the distal outcome's column name on the left.",
         call. = FALSE)
  }
  if (!identical(formula[[3L]], 1) && !identical(formula[[3L]], 1L)) {
    stop(
      "The distal model's right-hand side must be `1`. To model covariates as ",
      "well, fit them with tse_covariate() and pass that model to tse_distal().",
      call. = FALSE
    )
  }
  as.character(formula[[2L]])
}

#' match.arg(), or the inherited value when `x` is NULL
#' @noRd
.match_or_inherit <- function(x, inherited, choices, name) {
  if (is.null(x)) {
    return(inherited)
  }
  if (!is.character(x) || length(x) != 1L || !x %in% choices) {
    stop(sprintf("`%s` must be one of %s.", name, paste(dQuote(choices, FALSE), collapse = ", ")),
         call. = FALSE)
  }
  x
}

# -- covariate model methods -----------------------------------------------------

#' Class-membership probabilities from a covariate model
#'
#' The fitted class prior \eqn{P(X = t \mid Z)} for the rows of `newdata`, or
#' of the estimation data.
#'
#' @param object A `tseLCA_covariate` object from [tse_covariate()].
#' @param newdata Optional data frame with the covariates.
#' @param type `"prob"` (default) for the n x T probability matrix, or
#'   `"class"` for the most likely class.
#' @param ... Unused.
#' @return A matrix (rows with missing covariates are `NA`) or integer vector.
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' fit <- tse_covariate(tse_classify(m), ~ Zp)
#' predict(fit, newdata = data.frame(Zp = 1:5))
#' @export
predict.tseLCA_covariate <- function(object, newdata = NULL, type = c("prob", "class"), ...) {
  type <- match.arg(type)
  if (is.null(object$terms)) {
    stop("predict() needs a model fitted with tse_covariate().", call. = FALSE)
  }
  if (is.null(newdata)) newdata <- object[["classification", exact = TRUE]]$data
  tt <- stats::delete.response(object$terms)
  mf <- stats::model.frame(tt, newdata, xlev = object$xlevels, na.action = stats::na.pass)
  Z <- stats::model.matrix(tt, mf)
  coefs <- object$three_step
  eta <- cbind(0, Z %*% coefs)
  p <- exp(eta - apply(eta, 1L, max))
  p <- p / rowSums(p)
  colnames(p) <- c(object$ref, colnames(coefs))
  p <- p[, order(as.integer(sub("^C", "", colnames(p)))), drop = FALSE]
  rownames(p) <- rownames(newdata)
  if (type == "class") max.col(p, ties.method = "first") else p
}

#' Change the reference class of a covariate model
#'
#' Refits the model with another reference class of the multinomial logit.
#'
#' @param x A `tseLCA_covariate` object from [tse_covariate()].
#' @param ref The new reference class (number or label such as `"C2"`).
#' @param ... Unused.
#' @return The refitted `tseLCA_covariate` object.
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' fit <- tse_covariate(tse_classify(m), ~ Zp)
#' coef(stats::relevel(fit, ref = "C3"))
#' @export
relevel.tseLCA_covariate <- function(x, ref, ...) {
  if (is.null(x[["classification", exact = TRUE]])) {
    stop("relevel() needs a model fitted with tse_covariate().", call. = FALSE)
  }
  call <- x$call
  call$object <- x[["classification", exact = TRUE]]
  call$ref <- ref
  eval(call, parent.frame())
}

#' Wald tests of covariate terms
#'
#' Tests, term by term, that all multinomial-logit coefficients of a
#' covariate term are zero (for all non-reference classes), using the model's
#' variance matrix.
#'
#' @param object A `tseLCA_covariate` object from [tse_covariate()].
#' @param ... Unused.
#' @return An `anova` table with the Wald statistic, degrees of freedom, and
#'   p-value of each term.
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' anova(tse_covariate(tse_classify(m), ~ Zp))
#' @export
anova.tseLCA_covariate <- function(object, ...) {
  if (is.null(object$terms)) {
    stop("anova() needs a model fitted with tse_covariate().", call. = FALSE)
  }
  tt <- object$terms
  labels <- attr(tt, "term.labels")
  mf <- stats::model.frame(stats::delete.response(tt), object[["classification", exact = TRUE]]$data,
                           xlev = object$xlevels)
  assign <- attr(stats::model.matrix(tt, mf), "assign")
  b <- coef(object)
  V <- vcov(object)
  n_logit <- object$n_classes - 1L
  rows <- lapply(seq_along(labels), function(j) {
    idx <- which(rep(assign, n_logit) == j)
    bj <- b[idx]
    W <- as.numeric(t(bj) %*% solve(V[idx, idx, drop = FALSE], bj))
    c(Df = length(idx), Chisq = W, `Pr(>Chisq)` = stats::pchisq(W, length(idx), lower.tail = FALSE))
  })
  tab <- as.data.frame(do.call(rbind, rows))
  rownames(tab) <- labels
  structure(
    tab,
    heading = "Wald tests of covariate terms (all class contrasts)\n",
    class = c("anova", "data.frame")
  )
}
