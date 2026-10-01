# tseLCA/R/tseLCA.R
#
# One-call interface chaining the step-wise building blocks, and accessors
# for the components of a fitted model.

#' Three-step latent class analysis in one call
#'
#' Fits the measurement model, classifies the observations, and relates the
#' classes to covariates and/or a distal outcome, all from one formula. This
#' is a convenience wrapper around the step-wise functions [tse_lca()],
#' [tse_classify()], [tse_covariate()], and [tse_distal()], which remain
#' available for inspecting each step (see [measurement()]).
#'
#' @details
#' The formula has up to three parts: `indicators ~ covariates | outcome`.
#' * Left-hand side: the indicators, `cbind(Y1, Y2, ...)`.
#' * First right-hand side part: the covariates of class membership, with the
#'   usual formula syntax (`1` for none).
#' * Optional second part: the name of a distal outcome.
#'
#' For example, `cbind(Y1, Y2, Y3) ~ age + sex | income` relates the classes
#' to the covariates `age` and `sex` and to the distal outcome `income`;
#' `cbind(Y1, Y2, Y3) ~ 1 | income` has only the distal outcome; and
#' `cbind(Y1, Y2, Y3) ~ 1` fits the measurement model alone.
#'
#' The number of classes is chosen beforehand from the measurement model, for
#' example with `tse_lca(..., nclass = 1:6)`.
#'
#' @param formula `cbind(indicators) ~ covariates | distal outcome`; see
#'   Details.
#' @param data A data frame.
#' @param nclass Number of latent classes.
#' @param family Distribution of the distal outcome; see [tse_distal()].
#' @param method Step-3 estimator: `"ML"`, `"BCH"`, or `"none"`; see
#'   [tse_covariate()].
#' @param se Standard errors: `"corrected"` or `"robust"`.
#' @param assignment Step-2 class assignment: `"modal"` or `"proportional"`.
#' @param ref Reference class of the covariate model.
#' @param missing Missing indicator values: `"listwise"` or `"fiml"`; see
#'   [tse_lca()].
#' @param start Optional Step-1 starting values; see [tse_lca()].
#' @param control Estimation settings, see [tse_control()].
#'
#' @return A `tseLCA_measurement`, `tseLCA_covariate`, `tseLCA_distal`, or
#'   `tseLCA_both` object, depending on the formula. Its components are
#'   available with [measurement()], [classification()], [covariate()], and
#'   [distal()].
#'
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' d$Zo <- draw_Zo(d$X, bk2018_params$distal_params)
#'
#' fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
#' summary(fit)
#' measurement(fit)
#' classification(fit)
#'
#' # the same model, step by step
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' fc <- tse_covariate(tse_classify(m), ~ Zp)
#' fb <- tse_distal(fc, Zo ~ 1)
#' @export
tseLCA <- function(
  formula,
  data,
  nclass,
  family = "gaussian",
  method = c("ML", "BCH", "none"),
  se = c("corrected", "robust"),
  assignment = c("modal", "proportional"),
  ref = 1,
  missing = c("listwise", "fiml"),
  start = NULL,
  control = tse_control()
) {
  cl <- match.call()
  method <- match.arg(method)
  se <- match.arg(se)
  assignment <- match.arg(assignment)
  missing <- match.arg(missing)
  parts <- .split_tseLCA_formula(formula)
  if (length(nclass) != 1L) {
    stop(
      "`nclass` must be a single number of classes. Compare numbers of ",
      "classes with tse_lca(..., nclass = 1:6), for example.",
      call. = FALSE
    )
  }

  m <- tse_lca(parts$measurement, data = data, nclass = nclass, start = start,
               missing = missing, control = control)
  if (is.null(parts$covariates) && is.null(parts$outcome)) {
    m$call <- cl
    return(m)
  }
  clf <- tse_classify(m, assignment = assignment)
  fit <- if (!is.null(parts$covariates)) {
    tse_covariate(clf, parts$covariates, method = method, se = se, ref = ref)
  } else {
    clf
  }
  if (!is.null(parts$outcome)) {
    fit <- tse_distal(fit, stats::reformulate("1", response = parts$outcome),
                      family = family, method = method, se = se)
  }
  fit$call <- cl
  fit$formula <- formula
  fit
}

#' Split `cbind(Y...) ~ covariates | outcome` into its parts
#' @noRd
.split_tseLCA_formula <- function(formula) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula such as `cbind(Y1, Y2) ~ x | z`.", call. = FALSE)
  }
  F <- Formula::Formula(formula)
  len <- length(F)
  if (len[1L] != 1L || len[2L] > 2L) {
    stop(
      "`formula` must be `cbind(indicators) ~ covariates | distal outcome`, ",
      "with at most two right-hand side parts.",
      call. = FALSE
    )
  }
  measurement <- stats::formula(F, lhs = 1L, rhs = 0L)
  measurement[[3L]] <- 1
  cov <- stats::formula(F, lhs = 0L, rhs = 1L)
  covariates <- if (length(attr(stats::terms(cov), "term.labels")) > 0L) cov else NULL

  outcome <- NULL
  if (len[2L] == 2L) {
    out <- stats::formula(F, lhs = 0L, rhs = 2L)[[2L]]
    if (!is.name(out)) {
      stop("The distal outcome (after `|`) must be a single column name.", call. = FALSE)
    }
    outcome <- as.character(out)
  }
  list(measurement = measurement, covariates = covariates, outcome = outcome)
}

# -- accessors ---------------------------------------------------------------------

#' Components of a fitted tseLCA model
#'
#' Extract the step-wise components of a model fitted with [tseLCA()] or the
#' step-wise functions: the Step-1 measurement model, the Step-2
#' classification, and the Step-3 covariate and distal outcome models.
#'
#' @param x A fitted tseLCA object.
#' @param ... Unused.
#' @return `measurement()`: a `tseLCA_measurement` object.
#'   `classification()`: a `tseLCA_classify` object. `covariate()`: a
#'   `tseLCA_covariate` object. `distal()`: a `tseLCA_distal` object.
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' d$Zo <- draw_Zo(d$X, bk2018_params$distal_params)
#' fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
#' class_sizes(measurement(fit))
#' classification(fit)$D
#' coef(covariate(fit))
#' omnibus_test(distal(fit))
#' @export
measurement <- function(x, ...) UseMethod("measurement")

#' @rdname measurement
#' @export
measurement.tseLCA <- function(x, ...) {
  if (inherits(x, "tseLCA_measurement")) {
    return(x)
  }
  clf <- .classification_of(x)
  if (is.null(clf$measurement)) {
    stop("This model does not keep its measurement model; fit it with tse_lca().",
         call. = FALSE)
  }
  clf$measurement
}

#' @rdname measurement
#' @export
classification <- function(x, ...) UseMethod("classification")

#' @rdname measurement
#' @export
classification.tseLCA <- function(x, ...) .classification_of(x)

#' @rdname measurement
#' @export
covariate <- function(x, ...) UseMethod("covariate")

#' @rdname measurement
#' @export
covariate.tseLCA <- function(x, ...) {
  if (inherits(x, "tseLCA_both")) return(x$covariate)
  if (inherits(x, "tseLCA_covariate")) return(x)
  stop("This model has no covariate component.", call. = FALSE)
}

#' @rdname measurement
#' @export
distal <- function(x, ...) UseMethod("distal")

#' @rdname measurement
#' @export
distal.tseLCA <- function(x, ...) {
  if (inherits(x, "tseLCA_distal")) return(x)
  if (!inherits(x, "tseLCA_both")) {
    stop("This model has no distal outcome component.", call. = FALSE)
  }
  out <- c(
    x$distal,
    list(
      measurement_model = x$measurement_model,
      family = x$family,
      n_classes = x$n_classes,
      estimator = x$estimator,
      posteriors = x$posteriors,
      classifications = x$classifications,
      outcome = x$outcome
    )
  )
  class(out) <- c("tseLCA_distal", "tseLCA_structural", "tseLCA")
  out
}

#' The tseLCA_classify object behind a fitted model
#' @noRd
.classification_of <- function(x) {
  if (inherits(x, "tseLCA_classify")) return(x)
  clf <- x[["classification", exact = TRUE]]
  if (is.null(clf) && inherits(x, "tseLCA_both")) clf <- x$covariate[["classification", exact = TRUE]]
  if (is.null(clf)) {
    stop(
      "This model does not keep its classification; fit it with the step-wise ",
      "functions or tseLCA().",
      call. = FALSE
    )
  }
  clf
}
