# tseLCA/R/tse_classify.R
#
# Step 2 building block: class assignment and classification error.

#' Assign observations to latent classes (Step 2)
#'
#' Computes posterior class-membership probabilities from a fitted
#' measurement model, assigns observations to classes, and estimates the
#' classification-error probabilities \eqn{D_{ts} = P(W = s \mid X = t)}
#' between the true class \eqn{X} and the assigned class \eqn{W}. These error
#' probabilities are what the bias-adjusted Step-3 estimators (BCH and ML)
#' correct for.
#'
#' The measurement model is held fixed. With `newdata`, observations from
#' another sample are classified with it, e.g. to relate the classes to
#' covariates observed only in a subsample; the uncertainty of the
#' measurement model is then still that of the sample it was estimated on.
#'
#' @param object A measurement model from [tse_lca()] (two or more classes).
#' @param newdata Optional data frame to classify. Omitted: the data the
#'   measurement model was estimated on.
#' @param assignment `"modal"` (default): each observation is assigned to its
#'   most likely class. `"proportional"`: each observation is assigned to
#'   every class with its posterior probability as weight.
#' @param control Estimation settings; default: those of `object`. See
#'   [tse_control()].
#'
#' @return A `tseLCA_classify` object with components `posteriors` (n x T),
#'   `classifications` (modal classes), `weights` (the assignment weights
#'   \eqn{P(W = s \mid Y_i)}), `D` (the T x T classification-error matrix),
#'   `entropy.R2` (computed from these posteriors), and `data`. Pass it to
#'   the Step-3 functions.
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#'   three-step approaches. \emph{Political Analysis}, 18(4), 450--469.
#'   \doi{10.1093/pan/mpq025}
#'
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' cl <- tse_classify(m)
#' cl
#' tse_classify(m, assignment = "proportional")
#'
#' # classify another sample with the same measurement model
#' tse_classify(m, newdata = d[1:200, ])
#' @export
tse_classify <- function(
  object,
  newdata = NULL,
  assignment = c("modal", "proportional"),
  control = NULL
) {
  cl <- match.call()
  if (!inherits(object, "tseLCA_measurement")) {
    stop("`object` must be a measurement model from tse_lca().", call. = FALSE)
  }
  if (object$n_classes < 2L) {
    stop("Classification needs a model with at least two classes.", call. = FALSE)
  }
  assignment <- match.arg(assignment)
  data <- if (!is.null(newdata)) newdata else object$data
  if (is.null(data)) {
    stop(
      "This measurement model does not store its data; supply `newdata`.",
      call. = FALSE
    )
  }
  if (!is.data.frame(data)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  if (is.null(control)) {
    control <- if (!is.null(object$control)) object$control else tse_control()
  }
  s1 <- object$measurement_model
  missing <- .measurement_missing(object)

  opts <- .opts_from_control(
    control,
    incomplete = missing == "fiml",
    include.intercept = TRUE,
    use.modal.assignment = assignment == "modal",
    use.simple.cov = FALSE,
    use.bch = FALSE
  )
  rec <- .recode_indicators(data, s1$Y.names, s1$Y.levels)
  dat <- .prepare_data(rec$data, s1$Y.names, NULL, NULL, "gaussian", opts, s1$Y.levels)
  s2 <- .step2(dat, s1$fit0, object$n_classes, opts)$all

  post <- s2$p.xy
  colnames(post) <- paste0("C", seq_len(object$n_classes))
  structure(
    list(
      measurement_model = s1,
      measurement = object,
      step2 = s2,
      data = data,
      rows = dat$keep_Y,
      assignment = assignment,
      missing = missing,
      control = control,
      n_classes = object$n_classes,
      posteriors = post,
      classifications = max.col(post),
      weights = s2$w.is,
      D = .classification_error_matrix(s2$p.wx_mat),
      entropy.R2 = .entropy_R2(post),
      nobs = nrow(post),
      call = cl
    ),
    class = c("tseLCA_classify", "tseLCA")
  )
}

#' Missing-data handling of a measurement model ("listwise" or "fiml")
#'
#' Recorded by tse_lca(); for models from three_step(), inferred from whether
#' a FIML design mask was stored.
#' @noRd
.measurement_missing <- function(object) {
  if (!is.null(object$missing)) {
    return(object$missing)
  }
  if (!is.null(object$measurement_model$mDesign.exp)) "fiml" else "listwise"
}

#' T x T matrix `D[t, s]` = P(W = s | X = t) from the internal `pwx[s, t]`
#' @noRd
.classification_error_matrix <- function(pwx) {
  iT <- ncol(pwx)
  D <- t(pwx)
  dimnames(D) <- list(paste0("X=C", seq_len(iT)), paste0("W=C", seq_len(iT)))
  D
}

#' Entropy R^2 of a matrix of posterior class probabilities
#' @noRd
.entropy_R2 <- function(post) {
  iT <- ncol(post)
  p <- post[post > 0]
  1 - sum(-p * log(p)) / (nrow(post) * log(iT))
}

#' @rdname tse_classify
#' @param x A `tseLCA_classify` object.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#' @export
print.tseLCA_classify <- function(x, digits = max(3L, getOption("digits") - 4L), ...) {
  cat("Latent class assignment (Step 2)\n")
  cat(sprintf(
    "  Classes: %d   Assignment: %s   N: %d   Entropy R\u00b2: %.4f\n",
    x$n_classes, x$assignment, x$nobs, x$entropy.R2
  ))
  cat("\nClassification error probabilities P(W = s | X = t)\n")
  cat("(rows: true class X; columns: assigned class W)\n")
  print(round(x$D, digits))
  sizes <- rbind(
    estimated = class_sizes(x),
    assigned = colMeans(x$weights)
  )
  cat("\nClass proportions: estimated (measurement model) and assigned\n")
  print(round(sizes, digits))
  invisible(x)
}
