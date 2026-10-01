# tseLCA/R/control.R
#
# Estimation tuning settings, collected in one object (as glm.control() does
# for glm()).

#' Estimation settings for tseLCA models
#'
#' Collects the numerical settings of the three estimation steps. Pass the
#' result as the `control` argument of the model-fitting functions.
#'
#' @param step1.maxit Maximum number of EM iterations for the Step-1
#'   measurement model. A fit that reaches it is retried with twice as many.
#' @param step1.tol Convergence tolerance of the Step-1 EM algorithm (change
#'   in log-likelihood).
#' @param step1.restarts Number of additional random starts tried when the
#'   default Step-1 fit has an entropy R\eqn{^2} below `step1.restart.R2`;
#'   the fit with the highest log-likelihood is kept.
#' @param step1.restart.R2 Entropy R\eqn{^2} threshold that triggers
#'   `step1.restarts`.
#' @param n_init Optional number of Step-1 fits from independent random
#'   classifications, bypassing the default k-means initialization; the fit
#'   with the highest log-likelihood is kept. Unlike `step1.restarts`, these
#'   always run. `NULL` (default) uses the default initialization.
#' @param step3.maxit Maximum number of iterations for the Step-3
#'   (structural model) EM or Newton-Raphson algorithm.
#' @param step3.tol Convergence tolerance of the Step-3 algorithm.
#' @param boundary.tol Step-1 probabilities within this distance of 0 or 1
#'   are treated as fixed when computing the Step-1 variance.
#' @param hessian Step-3 information matrix used for standard errors.
#'   `"observed"` (default) uses the analytic Hessian, giving a sandwich
#'   variance that is robust to misspecification of the Step-3 model. `"opg"`
#'   uses the outer product of the case-wise scores, which relies on the
#'   information-matrix equality and is valid only when the Step-3 model is
#'   correctly specified. If the Hessian cannot be inverted, `"opg"` is used
#'   with a warning. Applies to ML covariate models.
#' @param verbose Logical. Print progress and convergence messages.
#'
#' @return A list of class `"tse_control"`.
#' @examples
#' tse_control()
#' tse_control(step1.maxit = 10000, n_init = 20)
#' @export
tse_control <- function(
  step1.maxit = 5000L,
  step1.tol = 1e-8,
  step1.restarts = 10L,
  step1.restart.R2 = 0.70,
  n_init = NULL,
  step3.maxit = 200L,
  step3.tol = 1e-6,
  boundary.tol = 1e-2,
  hessian = c("observed", "opg"),
  verbose = FALSE
) {
  count <- function(x, nm, allow_zero = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x == round(x) &&
      (x > 0 || (allow_zero && x == 0))
    if (!ok) {
      stop(sprintf(
        "`%s` must be a %s whole number.", nm,
        if (allow_zero) "non-negative" else "positive"
      ), call. = FALSE)
    }
    as.integer(x)
  }
  positive <- function(x, nm) {
    if (!(is.numeric(x) && length(x) == 1L && !is.na(x) && x > 0)) {
      stop(sprintf("`%s` must be a positive number.", nm), call. = FALSE)
    }
    x
  }
  if (!(is.numeric(step1.restart.R2) && length(step1.restart.R2) == 1L &&
        !is.na(step1.restart.R2) && step1.restart.R2 >= 0 && step1.restart.R2 <= 1)) {
    stop("`step1.restart.R2` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!(is.numeric(boundary.tol) && length(boundary.tol) == 1L &&
        !is.na(boundary.tol) && boundary.tol >= 0 && boundary.tol < 0.5)) {
    stop("`boundary.tol` must be a number in [0, 0.5).", call. = FALSE)
  }
  if (!(is.logical(verbose) && length(verbose) == 1L && !is.na(verbose))) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  structure(
    list(
      step1.maxit = count(step1.maxit, "step1.maxit"),
      step1.tol = positive(step1.tol, "step1.tol"),
      step1.restarts = count(step1.restarts, "step1.restarts", allow_zero = TRUE),
      step1.restart.R2 = step1.restart.R2,
      n_init = if (is.null(n_init)) NULL else count(n_init, "n_init"),
      step3.maxit = count(step3.maxit, "step3.maxit"),
      step3.tol = positive(step3.tol, "step3.tol"),
      boundary.tol = boundary.tol,
      hessian = match.arg(hessian),
      verbose = verbose
    ),
    class = "tse_control"
  )
}

#' @export
print.tse_control <- function(x, ...) {
  cat("tseLCA estimation settings\n")
  vals <- vapply(x, function(v) if (is.null(v)) "NULL" else format(v), character(1))
  print(noquote(cbind(value = vals)), right = FALSE)
  invisible(x)
}

#' Internal estimation options from a tse_control object
#'
#' Maps the user-facing control names onto the option names used by the
#' estimation internals, and adds the model-specification settings.
#' @noRd
.opts_from_control <- function(control, ...) {
  if (!inherits(control, "tse_control")) {
    stop("`control` must be created with tse_control().", call. = FALSE)
  }
  c(
    list(
      maxIter.measurement = control$step1.maxit,
      measurement.tol = control$step1.tol,
      iter.measurement = control$step1.restarts,
      R2.threshold = control$step1.restart.R2,
      n_init = control$n_init,
      em.maxIter = control$step3.maxit,
      covariate.tol = control$step3.tol,
      boundary.tol = control$boundary.tol,
      correct.spec = control$hessian == "opg",
      verbose = control$verbose
    ),
    list(...)
  )
}
