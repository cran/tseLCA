# tseLCA/R/tse_lca.R
#
# Step 1 building block: the latent class measurement model, and class
# enumeration over several numbers of classes.

#' Fit a latent class measurement model (Step 1)
#'
#' Estimates the measurement model: the class sizes and the class-conditional
#' item-response probabilities of the indicators on the left-hand side of
#' `formula`. This is the first step of three-step estimation; covariates and
#' distal outcomes are related to the classes afterwards, holding this model
#' fixed.
#'
#' With a vector `nclass`, a model is fitted for each number of classes and
#' the result is a class-enumeration table of fit statistics (see Details).
#'
#' @details
#' The number of classes is chosen from the measurement model alone, before
#' any structural variables are considered, typically by the BIC, the
#' interpretability of the classes, and their separation; see Nylund, Asparouhov, and
#' \enc{Muthén}{Muthen} (2007) and Masyn (2013). The enumeration table reports, for each
#' number of classes, the log-likelihood, number of free parameters, AIC,
#' BIC, sample-size adjusted BIC (SABIC; Sclove 1987), entropy R\eqn{^2}, and
#' the smallest estimated class proportion. The one-class model is the
#' independence model, fitted in closed form.
#'
#' Indicators may be factors, logicals, character, or numeric codes; their
#' categories are stored with the model and reused when the model is applied
#' to new data (see [predict.tseLCA_measurement()]) or in later steps.
#'
#' @param formula A formula `cbind(Y1, Y2, ...) ~ 1` naming the indicators.
#'   The measurement model has no covariates.
#' @param data A data frame.
#' @param nclass Number of latent classes, or a vector of numbers of classes
#'   to compare (e.g. `1:6`).
#' @param start Optional fixed starting point for the EM algorithm (single
#'   `nclass` only): an integer vector with one class per row of `data`, or a
#'   matrix of item-response probabilities P(Y = k | X = t) with one row per
#'   item category (in the order of the indicators) and one column per class,
#'   such as [item_probs()] of a fitted model (whose binary items have one
#'   row, P(Y = 1 | X = t)). Bypasses the default k-means initialization.
#' @param missing Handling of missing indicator values: `"listwise"`
#'   (default) drops rows with any missing indicator; `"fiml"` keeps rows with
#'   at least one observed indicator (full-information maximum likelihood).
#' @param control Estimation settings, see [tse_control()].
#'
#' @return For a single `nclass`, a `tseLCA_measurement` object (see
#'   [class_sizes()], [item_probs()], [posterior()], [predict()]; it keeps
#'   `data` for [tse_classify()]); for several,
#'   a `tseLCA_select` object: the enumeration table with the fitted models,
#'   see [best_model()].
#'
#' @references
#' Masyn, K. E. (2013). Latent class analysis and finite mixture modeling. In
#'   T. D. Little (Ed.), \emph{The Oxford Handbook of Quantitative Methods},
#'   Vol. 2, 551--611. Oxford University Press.
#'
#' Nylund, K. L., Asparouhov, T., & \enc{Muthén}{Muthen}, B. O. (2007). Deciding on the
#'   number of classes in latent class analysis and growth mixture modeling:
#'   A Monte Carlo simulation study. \emph{Structural Equation Modeling},
#'   14(4), 535--569. \doi{10.1080/10705510701575396}
#'
#' Sclove, S. L. (1987). Application of model-selection criteria to some
#'   problems in multivariate analysis. \emph{Psychometrika}, 52(3),
#'   333--343. \doi{10.1007/BF02294360}
#'
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#'
#' # Class enumeration
#' sel <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 1:4)
#' sel
#' plot(sel)
#'
#' # The selected model
#' m <- best_model(sel, criterion = "BIC")
#' m
#' class_sizes(m)
#' item_probs(m)
#' head(predict(m, newdata = d[1:5, ]))
#' @export
tse_lca <- function(
  formula,
  data,
  nclass,
  start = NULL,
  missing = c("listwise", "fiml"),
  control = tse_control()
) {
  cl <- match.call()
  missing <- match.arg(missing)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  Y.names <- .indicators_from_formula(formula)
  if (!(is.numeric(nclass) && length(nclass) >= 1L && !anyNA(nclass) &&
        all(nclass >= 1) && all(nclass == round(nclass)))) {
    stop("`nclass` must be one or more positive whole numbers.", call. = FALSE)
  }
  nclass <- sort(unique(as.integer(nclass)))
  if (!is.null(start) && length(nclass) > 1L) {
    stop("`start` can only be supplied with a single `nclass`.", call. = FALSE)
  }
  if (!is.null(start) && !is.null(control$n_init)) {
    stop("Supply either `start` or `control$n_init`, not both.", call. = FALSE)
  }
  # multilevLCA cannot start its EM algorithm from a given classification when
  # indicator values are missing (it stops with "sort_index(): detected NaN"),
  # so random starts and `start` are not available with FIML on incomplete data.
  if (missing == "fiml" && anyNA(data[Y.names])) {
    if (!is.null(start)) {
      stop("`start` cannot be used with missing = \"fiml\" when indicator values ",
           "are missing (multilevLCA cannot initialize from it).", call. = FALSE)
    }
    if (!is.null(control$n_init)) {
      warning("Random starts (`n_init`) are not available with missing = \"fiml\" ",
              "when indicator values are missing; using multilevLCA's default ",
              "initialization.", call. = FALSE)
      control$n_init <- NULL
    }
  }

  opts <- .opts_from_control(
    control,
    incomplete = missing == "fiml",
    include.intercept = TRUE,
    use.two.step = FALSE,
    get.twostep.vcov = FALSE,
    rebase = "C1",
    startval = start
  )
  rec <- .recode_indicators(data, Y.names)
  dat <- .prepare_data(rec$data, Y.names, NULL, NULL, "gaussian", opts, rec$levels)

  fit_one <- function(iT) {
    fit <- if (iT == 1L) {
      .fit_independence(dat)
    } else {
      s1 <- .fit_step1(rec$data, Y.names, iT, NULL, NULL, 1L, opts, Y.levels = rec$levels)
      s1 <- .attach_step1_data(s1, dat, 1L, fitted_here = TRUE)
      .new_measurement_fit(s1, dat, iT)
    }
    cl_k <- cl
    cl_k$nclass <- iT
    fit$call <- cl_k
    fit$formula <- formula
    fit$missing <- missing
    fit$control <- control
    fit$data <- data
    fit
  }

  if (length(nclass) == 1L) {
    return(fit_one(nclass))
  }
  fits <- lapply(nclass, fit_one)
  names(fits) <- as.character(nclass)
  structure(
    list(
      table = .enumeration_table(fits),
      fits = fits,
      call = cl
    ),
    class = "tseLCA_select"
  )
}

#' Indicator names from a measurement formula `cbind(Y1, ...) ~ 1`
#' @noRd
.indicators_from_formula <- function(formula) {
  usage <- "e.g. `cbind(Y1, Y2, Y3) ~ 1`"
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be two-sided, ", usage, ".", call. = FALSE)
  }
  if (!identical(formula[[3L]], 1) && !identical(formula[[3L]], 1L)) {
    stop(
      "The measurement model has no covariates: its right-hand side must be `1` (",
      usage, "). Covariates and distal outcomes are related to the classes in ",
      "Step 3, with the measurement model held fixed.",
      call. = FALSE
    )
  }
  lhs <- formula[[2L]]
  if (!(is.call(lhs) && identical(lhs[[1L]], as.name("cbind")))) {
    stop("The indicators must be given as `cbind(...)` on the left-hand side, ",
         usage, ".", call. = FALSE)
  }
  args <- as.list(lhs)[-1L]
  if (!all(vapply(args, is.name, logical(1)))) {
    stop("The indicators in `cbind(...)` must be column names.", call. = FALSE)
  }
  items <- vapply(args, as.character, character(1))
  if (length(items) < 2L) {
    stop("A latent class model needs at least two indicators.", call. = FALSE)
  }
  if (anyDuplicated(items)) {
    stop("Duplicated indicators: ", paste(unique(items[duplicated(items)]), collapse = ", "),
         call. = FALSE)
  }
  items
}

#' One-class (independence) measurement model, in closed form
#'
#' Item-response probabilities are the observed category proportions of each
#' item (over the rows where it is observed). Returns a tseLCA_measurement
#' object whose `fit0` mimics the multilevLCA layout.
#' @noRd
.fit_independence <- function(dat) {
  Y <- dat$Y.obs
  M <- if (!is.null(dat$mDesign)) dat$mDesign else matrix(1, nrow(Y), ncol(Y))
  ivItemcat <- dat$ivItemcat
  counts <- colSums(Y * M)
  item_of <- rep(seq_along(ivItemcat), ivItemcat)
  n_h <- tapply(counts, item_of, sum)[item_of]
  p <- counts / n_h
  llik <- sum(ifelse(counts > 0, counts * log(p), 0))

  Y.names <- dat$Y.names
  cat_labels <- unlist(lapply(seq_along(ivItemcat), function(h) {
    if (ivItemcat[h] == 2L) {
      Y.names[h]
    } else {
      paste0(Y.names[h], ".", seq_len(ivItemcat[h]) - 1L)
    }
  }))
  keep <- unlist(lapply(seq_along(ivItemcat), function(h) {
    idx <- which(item_of == h)
    if (ivItemcat[h] == 2L) idx[2L] else idx
  }))
  mPhi <- matrix(p[keep], ncol = 1L, dimnames = list(sprintf("P(%s|C)", cat_labels), "C1"))

  npar <- sum(ivItemcat - 1L)
  n <- nrow(Y)
  fit0 <- list(
    vPi = matrix(1, 1L, 1L, dimnames = list("P(C1)", "")),
    mPhi = mPhi,
    LLKSeries = matrix(llik),
    AIC = -2 * llik + 2 * npar,
    BIC = -2 * llik + npar * log(n),
    R2entr = NA_real_
  )
  s1 <- list(fit0 = fit0, fitZ = NULL)
  s1 <- .attach_step1_data(s1, dat, 1L, fitted_here = TRUE)
  .new_measurement_fit(s1, dat, 1L)
}

#' Class-enumeration table from a list of measurement models
#' @noRd
.enumeration_table <- function(fits) {
  rows <- lapply(fits, function(f) {
    ll <- logLik(f)
    n <- attr(ll, "nobs")
    k <- attr(ll, "df")
    data.frame(
      nclass = f$n_classes,
      logLik = as.numeric(ll),
      npar = k,
      AIC = -2 * as.numeric(ll) + 2 * k,
      BIC = -2 * as.numeric(ll) + k * log(n),
      SABIC = -2 * as.numeric(ll) + k * log((n + 2) / 24),
      entropy.R2 = if (f$n_classes > 1L) f$R2entr else NA_real_,
      min.class = min(class_sizes(f)),
      nobs = n
    )
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  tab
}

# -- tseLCA_select methods -----------------------------------------------------

#' Class enumeration results
#'
#' Methods for the `tseLCA_select` object returned by [tse_lca()] with
#' several numbers of classes. `best_model()` returns the fitted model that
#' minimizes an information criterion; `x[[k]]` returns the `k`-class model.
#'
#' @param x,object A `tseLCA_select` object.
#' @param criterion Information criterion to minimize: `"BIC"` (default),
#'   `"AIC"`, or `"SABIC"`.
#' @param i Number of classes of the model to extract.
#' @param digits Number of significant digits to print.
#' @param which Criteria to plot.
#' @param ... Further arguments passed to [graphics::matplot()] (`plot`) or
#'   unused.
#' @return `best_model()` and `[[`: a `tseLCA_measurement` object.
#'   `as.data.frame()`: the enumeration table. `print()`, `plot()`: `x`,
#'   invisibly.
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' sel <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 1:4)
#' as.data.frame(sel)
#' best_model(sel)
#' sel[[2]]
#' @export
best_model <- function(object, ...) UseMethod("best_model")

#' @rdname best_model
#' @export
best_model.tseLCA_select <- function(object, criterion = c("BIC", "AIC", "SABIC"), ...) {
  criterion <- match.arg(criterion)
  tab <- .subset2(object, "table")
  .subset2(object, "fits")[[which.min(tab[[criterion]])]]
}

#' @rdname best_model
#' @export
`[[.tseLCA_select` <- function(x, i, ...) {
  fits <- .subset2(x, "fits")
  k <- as.character(i)
  if (!k %in% names(fits)) {
    stop(sprintf(
      "No %s-class model; fitted: %s.", k, paste(names(fits), collapse = ", ")
    ), call. = FALSE)
  }
  fits[[k]]
}

#' @rdname best_model
#' @export
as.data.frame.tseLCA_select <- function(x, ...) .subset2(x, "table")

#' @rdname best_model
#' @export
print.tseLCA_select <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  tab <- .subset2(x, "table")
  cat("Latent class enumeration (measurement model)\n\n")
  shown <- tab[, c("nclass", "logLik", "npar", "AIC", "BIC", "SABIC", "entropy.R2", "min.class")]
  best <- vapply(c("AIC", "BIC", "SABIC"), function(cr) which.min(tab[[cr]]), integer(1))
  for (cr in names(best)) {
    shown[[cr]] <- paste0(format(round(tab[[cr]], 2), nsmall = 2), ifelse(seq_len(nrow(tab)) == best[[cr]], "*", " "))
  }
  print(shown, digits = digits, row.names = FALSE, right = TRUE)
  cat("\n* smallest value of each criterion. N =", tab$nobs[1L], "\n")
  invisible(x)
}

#' @rdname best_model
#' @export
plot.tseLCA_select <- function(x, which = c("AIC", "BIC", "SABIC"), ...) {
  tab <- .subset2(x, "table")
  which <- match.arg(which, c("AIC", "BIC", "SABIC"), several.ok = TRUE)
  ic <- as.matrix(tab[, which, drop = FALSE])
  graphics::matplot(
    tab$nclass, ic, type = "b", pch = seq_along(which), lty = seq_along(which),
    col = 1, xlab = "Number of classes", ylab = "Information criterion",
    xaxt = "n", ...
  )
  graphics::axis(1, at = tab$nclass)
  graphics::legend("topright", legend = which, pch = seq_along(which),
                   lty = seq_along(which), bty = "n")
  invisible(x)
}

# -- tseLCA_measurement: prediction ---------------------------------------------

#' Class membership predictions from a measurement model
#'
#' Posterior class-membership probabilities P(X = t | Y) for the rows of
#' `newdata` (or the estimation sample), or their modal class. Indicators in
#' `newdata` are coded with the categories stored in the model; missing
#' indicator values are skipped (the posterior uses the observed ones), and
#' rows with no observed indicator get `NA`.
#'
#' @param object A `tseLCA_measurement` object.
#' @param newdata Optional data frame with the indicator columns. Omitted:
#'   the estimation sample.
#' @param type `"posterior"` (default) for an n x T matrix of probabilities,
#'   or `"class"` for the modal class of each row.
#' @param ... Unused.
#' @return A matrix (`type = "posterior"`) or integer vector (`"class"`).
#' @examples
#' d <- generate_data(300, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#' predict(m, newdata = d[1:5, ])
#' predict(m, newdata = d[1:5, ], type = "class")
#' @export
predict.tseLCA_measurement <- function(object, newdata = NULL, type = c("posterior", "class"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    post <- posterior(object)
  } else {
    s1 <- object$measurement_model
    Y.names <- s1$Y.names
    rec <- .recode_indicators(newdata, Y.names, s1$Y.levels)
    Y <- as.matrix(rec$data[, Y.names, drop = FALSE])
    observed <- rowSums(!is.na(Y)) > 0L
    Y_exp <- expand_Y(Y, s1$ivItemcat)
    mDesign <- (!is.na(Y_exp)) * 1L
    Y_exp[is.na(Y_exp)] <- 0
    post <- matrix(NA_real_, nrow(Y), object$n_classes,
                   dimnames = list(rownames(newdata), paste0("C", seq_len(object$n_classes))))
    if (any(observed)) {
      post[observed, ] <- step1_posteriors(
        Y_exp[observed, , drop = FALSE], mDesign[observed, , drop = FALSE],
        s1$fit0, s1$ivItemcat
      )
    }
  }
  if (type == "class") max.col(post, ties.method = "first") else post
}

#' @rdname predict.tseLCA_measurement
#' @export
fitted.tseLCA_measurement <- function(object, ...) posterior(object)

#' @export
formula.tseLCA <- function(x, ...) {
  if (is.null(x$formula)) {
    stop("This model was not fitted from a formula.", call. = FALSE)
  }
  x$formula
}

#' Use a measurement model with given parameters (Step 1)
#'
#' Creates a measurement model from given class sizes and item-response
#' probabilities, evaluated on `data`, without estimating it. This allows
#' Steps 2 and 3 to be based on a measurement model estimated elsewhere: in
#' another program, reported in a publication, or saved from an earlier
#' analysis.
#'
#' The Step-1 variance used for corrected standard errors in Step 3 is
#' computed on `data` at the given parameters, which is valid when they are
#' the maximum likelihood estimates for `data` (e.g. a model estimated on these
#' data and saved). For parameters estimated on another sample, use
#' `se = "robust"` in Step 3, or refit with [tse_lca()] using the parameters as
#' `start`.
#'
#' @param formula `cbind(Y1, Y2, ...) ~ 1`, as in [tse_lca()].
#' @param data A data frame.
#' @param class_sizes Class proportions, one per class (they are normalized
#'   to sum to one).
#' @param item_probs Item-response probabilities in the layout of
#'   [item_probs()]: one column per class, and one row per binary item
#'   (\eqn{P(Y = 1 \mid X = t)}, where 1 is the item's second category) or per
#'   category of a polytomous item (\eqn{P(Y = k \mid X = t)}), in the order
#'   of the indicators.
#' @param missing,control As in [tse_lca()].
#'
#' @return A `tseLCA_measurement` object, usable like one from [tse_lca()].
#' @examples
#' d <- generate_data(500, "high", "covariate", seed = 1)
#' m <- tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d, nclass = 3)
#'
#' # the same measurement model from its parameters
#' m2 <- as_tse_lca(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1, data = d,
#'                  class_sizes = class_sizes(m), item_probs = item_probs(m))
#' all.equal(logLik(m2), logLik(m), tolerance = 1e-6)
#' coef(tse_covariate(tse_classify(m2), ~ Zp))
#' @export
as_tse_lca <- function(
  formula,
  data,
  class_sizes,
  item_probs,
  missing = c("listwise", "fiml"),
  control = tse_control()
) {
  cl <- match.call()
  missing <- match.arg(missing)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  Y.names <- .indicators_from_formula(formula)
  opts <- .opts_from_control(control, incomplete = missing == "fiml",
                             include.intercept = TRUE)
  rec <- .recode_indicators(data, Y.names)
  dat <- .prepare_data(rec$data, Y.names, NULL, NULL, "gaussian", opts, rec$levels)
  ivItemcat <- dat$ivItemcat

  pi_t <- as.numeric(class_sizes)
  iT <- length(pi_t)
  if (iT < 2L || anyNA(pi_t) || any(pi_t <= 0)) {
    stop("`class_sizes` must be two or more positive proportions.", call. = FALSE)
  }
  pi_t <- pi_t / sum(pi_t)
  phi <- as.matrix(item_probs)
  n_rows <- sum(ifelse(ivItemcat == 2L, 1L, ivItemcat))
  if (!identical(dim(phi), c(n_rows, iT))) {
    stop(sprintf(
      "`item_probs` must be a %d x %d matrix (item rows as in item_probs(), one column per class).",
      n_rows, iT
    ), call. = FALSE)
  }
  if (anyNA(phi) || any(phi <= 0 | phi >= 1)) {
    stop("`item_probs` must be probabilities strictly between 0 and 1.", call. = FALSE)
  }
  full <- expand_Phi(phi, ivItemcat)
  item_of <- rep(seq_along(ivItemcat), ivItemcat)
  sums <- apply(full, 2L, function(p) tapply(p, item_of, sum))
  if (any(abs(sums - 1) > 1e-6)) {
    stop("The category probabilities of each polytomous item must sum to one in each class.",
         call. = FALSE)
  }

  labels <- unlist(lapply(seq_along(ivItemcat), function(h) {
    if (ivItemcat[h] == 2L) Y.names[h] else paste0(Y.names[h], ".", seq_len(ivItemcat[h]) - 1L)
  }))
  dimnames(phi) <- list(sprintf("P(%s|C)", labels), paste0("C", seq_len(iT)))

  log_joint <- sweep(log_lik_matrix(dat$Y.obs, full, dat$mDesign), 2L, log(pi_t), "+")
  row_max <- apply(log_joint, 1L, max)
  llik <- sum(row_max + log(rowSums(exp(log_joint - row_max))))
  npar <- (iT - 1L) + iT * sum(ivItemcat - 1L)
  fit0 <- list(
    vPi = matrix(pi_t, iT, 1L, dimnames = list(sprintf("P(C%d)", seq_len(iT)), "")),
    mPhi = phi,
    LLKSeries = matrix(llik),
    AIC = -2 * llik + 2 * npar,
    BIC = -2 * llik + npar * log(nrow(dat$Y.obs))
  )
  fit0$R2entr <- .entropy_R2(step1_posteriors(dat$Y.obs, dat$mDesign, fit0, ivItemcat))

  s1 <- .attach_step1_data(list(fit0 = fit0, fitZ = NULL), dat, 1L, fitted_here = TRUE)
  fit <- .new_measurement_fit(s1, dat, iT)
  fit$call <- cl
  fit$formula <- formula
  fit$missing <- missing
  fit$control <- control
  fit$data <- data
  fit
}
