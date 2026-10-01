#' Three-step LCA estimation with covariates and/or distal outcomes
#'
#' Fits a three-step latent class model through the following steps:
#' \enumerate{
#'   \item \strong{Measurement model}: estimates latent class parameters
#'     (\eqn{\pi}, \eqn{\phi}) using \pkg{multilevLCA}
#'     (Lyrvall et al., 2025).
#'   \item \strong{Classification-error matrix}: computes posterior class
#'     probabilities and the T x T misclassification probability matrix
#'     \eqn{P(W = s \mid X = t)}, with standard errors corrected for
#'     classification-error propagation (Bakk, Oberski & Vermunt, 2014).
#'   \item \strong{Structural model}: estimates covariate effects using
#'     two-step starting values (Bakk & Kuha, 2018) and/or distal outcome
#'     means following Bakk, Tekle & Vermunt (2013), with the ML correction (Vermunt, 2010) or BCH correction
#'     (Bolck, Croon & Hagenaars, 2004). See
#'     \code{vignette("tseLCA", package = "tseLCA")} for a worked example.
#' }
#'
#' @param data A data.frame containing all columns referenced by \code{Y.names},
#'   \code{Zp.names}, and \code{Zo.name}.
#' @param Y.names Character vector of indicator column names. Indicators
#'   may be factors, logicals, character, or numeric codes (see [tse_lca()]).
#' @param n_classes Integer. Number of latent classes.
#' @param Zp.names Character vector of covariate column names, or \code{NULL}
#'   for a measurement-only fit. Default \code{NULL}.
#' @param Zo.name Single character name of the distal outcome column, or
#'   \code{NULL}. Default \code{NULL}.
#' @param step1 Pre-fitted Step-1 object (output of [tseLCA::lca_step1()] or a
#'   prior \code{three_step()} call), or \code{NULL} to run Step 1 internally.
#'   Default \code{NULL}.
#' @param startval Optional starting classification for the Step-1
#'   measurement model, either an integer vector of length \code{nrow(data)}
#'   (a class assignment \code{1..n_classes} for every row) or a numeric
#'   matrix of conditional item-response probabilities
#'   \eqn{P(Y_h = k \mid X = t)} (one row per item-category pair in
#'   \code{Y.names} order, one column per class) from which a classification
#'   is derived internally. See [lca_step1_startval()] for the full
#'   description of both forms and typical sources (an external solver run
#'   with many random starts, or a published item-response table).
#'   \pkg{multilevLCA}'s default initialization (k-means on principal
#'   components) is deterministic given the data and can converge to a local
#'   optimum of the Step-1 log-likelihood; supplying \code{startval} bypasses
#'   it (\code{kmea = FALSE} with the classification injected as
#'   multilevLCA's \code{startval}). Mutually exclusive with \code{step1} and
#'   \code{n_init}. Default \code{NULL}.
#' @param n_init Optional positive integer. If supplied, fits the Step-1
#'   measurement model \code{n_init} times from independent uniform-random
#'   classifications (\code{kmea = FALSE}, not multilevLCA's k-means-on-PCA
#'   path) and keeps the fit with the highest log-likelihood -- the
#'   unconditional multi-start analog of \code{n_init} in \pkg{StepMix} or
#'   \code{nrep} in \pkg{poLCA}. This is distinct from
#'   \code{iter.measurement}, which reruns multilevLCA's own k-means
#'   initialization, and only when the entropy R\eqn{^2} of the default fit
#'   is below \code{R2.threshold}; \code{n_init} restarts always run.
#'   Mutually exclusive with \code{step1} and \code{startval}. Default
#'   \code{NULL}.
#' @param use.two.step Logical. Initialize Step-3 from two-step estimates.
#'   Default \code{TRUE}.
#' @param use.modal.assignment Logical. Use modal (hard) class assignments in
#'   Step 2 and 3. \code{FALSE} uses soft posterior weights. Default \code{TRUE}.
#' @param include.intercept Logical. Prepend an intercept column to the
#'   covariate design matrix. Default \code{TRUE}.
#' @param use.simple.cov Logical. Skip the Step-1 measurement-uncertainty
#'   correction and return only the robust sandwich variance. Faster but
#'   underestimates standard errors when class separation is low. Default
#'   \code{FALSE}.
#' @param incomplete Logical. FIML for partially missing indicators. See the
#'   \code{Missing Data} section of \code{vignette("tseLCA", package = "tseLCA")}.
#'   Default \code{FALSE}.
#' @param boundary.tol Scalar. Parameters within this tolerance of 0 or 1 are
#'   treated as fixed when computing the Step-1 variance matrix for numerical stability. Default
#'   \code{1e-2}.
#' @param maxIter.measurement Integer. Maximum EM iterations for Step 1.
#'   Default \code{5000L}.
#' @param measurement.tol Scalar. Convergence tolerance for the Step-1 EM
#'   algorithm. Default \code{1e-8}.
#' @param covariate.tol Scalar. Convergence tolerance for the Step-3
#'   Newton-Raphson or EM algorithm. Default \code{1e-6}.
#' @param iter.measurement Integer. Number of random restarts triggered when
#'   the Step-1 entropy R\eqn{^2} falls below \code{R2.threshold}. Default
#'   \code{10L}.
#' @param R2.threshold Scalar. Entropy R\eqn{^2} threshold below which Step-1
#'   random restarts are triggered. Default \code{0.70}.
#' @param use.bch Logical. Use the BCH estimator in Step 3 (default: the ML
#'   estimator). May error if BCH weights induce a non-positive semi-definite Hessian in the third step (common in cases of low separation). Default \code{FALSE}.
#' @param em.maxIter Integer. Maximum EM iterations for the Step-3 covariate
#'   or distal outcome model. Default \code{200L}.
#' @param get.twostep.vcov Logical. If \code{TRUE}, obtain \pkg{multilevLCA}'s
#'   bias-corrected variance-covariance matrix for the two-step gamma estimates
#'   and store it in \code{$two_step_vcov}. If the \code{fitZ} object passed
#'   through \code{step1} already contains a \code{Varmat_cor} (from a prior
#'   [fitZ_from_multiLCA()] or plain \code{multiLCA} call), it is attached
#'   automatically even when \code{get.twostep.vcov = FALSE}. Default
#'   \code{FALSE}.
#' @param rebase Character (e.g. \code{"C1"}, \code{"C2"}) or integer
#'   specifying which latent class to use as the reference category in the
#'   multinomial logit. The measurement model is permuted so this class becomes
#'   column 1 before any structural estimation. Default \code{"C1"}.
#' @param family Character. Distal outcome family: one of \code{"gaussian"}
#'   (class means), \code{"poisson"} (log-rates), \code{"binomial"}
#'   (logits), or \code{"multinomial"} (a saturated model for a nominal
#'   categorical outcome with 2 or more categories -- \code{Zo.name} may be
#'   a factor, character, or integer column; categories are taken from
#'   \code{sort(unique(data[[Zo.name]]))} with \code{factor()}). For
#'   \code{"multinomial"}, \code{coef()} returns a \code{T x C} matrix of
#'   class-conditional category probabilities
#'   \eqn{\hat\pi_{tc} = P(Zo = c \mid X = t)} (rows sum to 1), not a
#'   length-\code{T} vector, and \code{vcov()} returns its
#'   \code{(T*C) x (T*C)} sandwich covariance (necessarily singular, since
#'   each class's row sums to 1 -- see \code{\link{omnibus_test}()} for a
#'   Wald test that accounts for this). Unlike \code{"binomial"}, whose
#'   \code{coef()}/\code{vcov()} are on the logit scale, \code{"multinomial"}
#'   reports \code{coef()}/\code{vcov()} directly on the probability scale,
#'   so \code{Std.Error} is directly interpretable without a delta-method
#'   back-transform -- but a symmetric interval
#'   \code{Estimate +/- 1.96*Std.Error} can fall outside \eqn{[0, 1]} for a
#'   probability near a boundary, the same well-known limitation as a naive
#'   Wald interval for any sample proportion. The \code{z.value}/\code{p.value}
#'   columns \code{summary()}/\code{print()} show for this family test each
#'   probability against 0, which is rarely the question of interest;
#'   \code{\link{omnibus_test}()} is the intended, boundary-safe test of
#'   whether the outcome's distribution differs across classes. Combining
#'   \code{family = "multinomial"} with both \code{Zp.names} and
#'   \code{Zo.name} fully
#'   propagates both Step-1 measurement and Step-3 covariate uncertainty
#'   under \code{use.simple.cov = FALSE}, the same as the other families.
#'   Default \code{"gaussian"}.
#' @param correct.spec Logical. Estimate the Step-3 information matrix by the
#'   outer product of the case-wise scores, not the observed-data Hessian.
#'   Valid only when the Step-3 model is correctly specified; the default
#'   observed-Hessian sandwich is robust to misspecification. Default
#'   \code{FALSE}.
#' @param verbose Logical. Print convergence messages. Default \code{FALSE}.
#'
#' @return An S3 object of class \code{tseLCA}. The subclass depends on which
#'   models were estimated:
#'   \describe{
#'     \item{`tseLCA_measurement`}{Returned when neither \code{Zp.names} nor
#'       \code{Zo.name} is supplied. Contains the following elements:
#'       \describe{
#'         \item{`measurement_model`}{Step-1 output list from [tseLCA::lca_step1()].}
#'         \item{`llik`}{Final Step-1 log-likelihood.}
#'         \item{`AIC`, `BIC`}{Information criteria from the measurement model.}
#'         \item{`R2entr`}{Entropy R\eqn{^2} of the measurement model.}
#'         \item{`n_classes`}{Number of latent classes.}
#'         \item{`posteriors`}{n x T matrix of soft posterior class probabilities.}
#'         \item{`classifications`}{Length-n integer vector of modal class assignments.}
#'       }
#'     }
#'     \item{`tseLCA_covariate`}{Returned when \code{Zp.names} is supplied and
#'       \code{Zo.name} is \code{NULL}. Contains all elements of
#'       \code{tseLCA_measurement} plus:
#'       \describe{
#'         \item{`three_step`}{(Q+1) x (T-1) matrix of Step-3 gamma coefficients.}
#'         \item{`three_step_vcov`}{(Q+1)(T-1) x (Q+1)(T-1) variance-covariance matrix
#'           for \code{three_step}, with measurement-uncertainty correction
#'           unless \code{use.simple.cov = TRUE}.}
#'         \item{`two_step`}{(Q+1) x (T-1) matrix of two-step starting values, or
#'           \code{NULL} if \code{use.two.step = FALSE}.}
#'         \item{`two_step_vcov`}{\pkg{multilevLCA} bias-corrected vcov for the
#'           two-step estimates, or \code{NULL}.}
#'         \item{`estimator`}{Character: \code{"ML"} or \code{"BCH"}.}
#'         \item{`entropy.R2`}{Covariate-adjusted entropy R\eqn{^2}.}
#'         \item{`llik`}{Profile log-likelihood
#'           \eqn{\sum_i \log \sum_t P(X=t|Z_{p,i};\hat{\gamma}) P(Y_i|X=t;\hat{\phi})},
#'           with Step-1 parameters \eqn{\hat{\phi}} held fixed. By construction
#'           smaller than the equivalent one-step MLE likelihood.}
#'       }
#'     }
#'     \item{`tseLCA_distal`}{Returned when \code{Zo.name} is supplied and
#'       \code{Zp.names} is \code{NULL}. Contains:
#'       \describe{
#'         \item{`three_step`}{Named length-T vector of Step-3 distal outcome
#'           parameters (means, log-rates, or logits depending on
#'           \code{family}) -- or, for \code{family = "multinomial"}, a
#'           \code{T x C} matrix of class-conditional category probabilities
#'           (rows sum to 1).}
#'         \item{`three_step_vcov`}{T x T variance-covariance matrix for
#'           \code{three_step}, named \code{mu_C1} through \code{mu_CT} --
#'           or, for \code{family = "multinomial"}, a \code{(T*C) x (T*C)}
#'           (necessarily rank-deficient) matrix named \code{"C{t}:{category}"}.}
#'         \item{`three_step.llik`}{Step-3 distal log-likelihood
#'           \eqn{\log P(Z_o|X=t)} at converged estimates.}
#'         \item{`llik`}{Profile log-likelihood
#'           \eqn{\sum_i \log \sum_t P(X=t|\hat{\pi}) P(Z_{o,i}|X=t;\hat{\mu}) P(Y_i|X=t;\hat{\phi})},
#'           with Step-1 parameters \eqn{\hat{\pi}, \hat{\phi}} held fixed.
#'           By construction smaller than the equivalent one-step MLE likelihood.}
#'         \item{`AIC`}{Akaike information criterion based on \code{llik}.}
#'         \item{`BIC`}{Bayesian information criterion based on \code{llik},
#'           using the number of distal-complete observations.}
#'         \item{`family`}{Character. The distal outcome family used.}
#'         \item{`estimator`}{Character: \code{"ML"} or \code{"BCH"}.}
#'         \item{`posteriors`}{n x T soft posterior matrix.}
#'         \item{`classifications`}{Length-n modal class assignment vector.}
#'       }
#'     }
#'     \item{`tseLCA_both`}{Returned when both \code{Zp.names} and
#'       \code{Zo.name} are supplied. Contains:
#'       \describe{
#'         \item{`covariate`}{A \code{tseLCA_covariate}-structured sub-list
#'           (see above), including \code{llik}, \code{AIC}, \code{BIC},
#'           \code{entropy.R2}.}
#'         \item{`distal`}{A \code{tseLCA_distal}-structured sub-list
#'           (see above), including \code{llik}, \code{AIC}, \code{BIC},
#'           \code{three_step.llik}.}
#'         \item{`family`, `n_classes`, `estimator`}{Shared top-level fields.}
#'         \item{`posteriors`, `classifications`}{Shared n x T posterior
#'           matrix and length-n modal class vector.}
#'       }
#'     }
#'   }
#'
#' @references
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#'   between latent class membership and external variables using bias-adjusted
#'   three-step approaches. \emph{Sociological Methodology}, 43(1), 272--311.
#'   \doi{10.1177/0081175012470644}
#'
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#'   classes and external variables. \emph{Psychometrika}, 83(4), 871--892.
#'   \doi{10.1007/s11336-017-9592-7}
#'
#' Bakk, Z., Pohle, M. J., & Kuha, J. (2025). Bias-adjusted three-step
#'   estimation of structural models for latent classes. \emph{Multivariate
#'   Behavioral Research}. \doi{10.1080/00273171.2025.2473935}
#'
#' @seealso \code{vignette("tseLCA", package = "tseLCA")} for a full worked
#'   example; [tseLCA::lca_step1()] for standalone Step-1 estimation
#'   (including from an externally supplied starting classification, with its
#'   own `startval` argument); [fitZ_from_fit0()] and [fitZ_from_multiLCA()]
#'   for two-step covariate estimation.
#'
#' @examples
#' d <- generate_data(n = 200, separation = "high",
#'                    scenario = "covariate", seed = 1)
#'
#' # Measurement model only
#' fit_m <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' summary(fit_m)
#'
#' # ML three-step with simple SEs (fast)
#' fit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                   Zp.names = "Zp", use.simple.cov = TRUE)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#'
#' # Full measurement-uncertainty correction (see vignette for interpretation)
#' fit_cor <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.simple.cov = FALSE,
#'                       use.modal.assignment = FALSE)
#' summary(fit_cor)
#'
#' # BCH estimator
#' fit_bch <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zp.names = "Zp", use.bch = TRUE,
#'                       use.simple.cov = TRUE)
#' summary(fit_bch)
#'
#' # Change reference class
#' fit_c2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                      Zp.names = "Zp", use.simple.cov = TRUE,
#'                      rebase = "C2")
#' summary(fit_c2)
#'
#' # Gaussian distal outcome
#' d2 <- generate_data(200, "high", "distal", seed = 2)
#' fit_dis <- three_step(d2, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zo.name = "Zo", family = "gaussian",
#'                       use.simple.cov = TRUE)
#' summary(fit_dis)
#'
#' # Nominal categorical distal outcome (3+ categories): coef() returns a
#' # T x C matrix of class-conditional category probabilities; omnibus_test()
#' # gives a single Wald test of whether the category distribution differs
#' # across classes at all.
#' d2$Zcat <- factor(sample(c("low", "mid", "high"), nrow(d2), replace = TRUE))
#' fit_cat <- three_step(d2, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       Zo.name = "Zcat", family = "multinomial",
#'                       use.simple.cov = TRUE)
#' coef(fit_cat)
#' omnibus_test(fit_cat)
#'
#' # Pass a pre-fitted measurement model to skip Step 1
#' fit_step1 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3)
#' fit2 <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                    Zp.names = "Zp", step1 = fit_step1,
#'                    use.simple.cov = TRUE)
#' summary(fit2)
#'
#' # Supply an external starting classification for Step 1 (bypasses
#' # multilevLCA's k-means-on-PCA initialization; here we use the DGP's own
#' # true classes as a stand-in for e.g. a StepMix solution with many
#' # random starts)
#' fit_ext <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                       startval = d$X, use.simple.cov = TRUE)
#' summary(fit_ext)
#'
#' # Many random-classification restarts for Step 1, keeping the best
#' # (analogous to n_init in StepMix or nrep in poLCA)
#' fit_ninit <- three_step(d, Y.names = paste0("Y", 1:6), n_classes = 3,
#'                         n_init = 20L, use.simple.cov = TRUE)
#' summary(fit_ninit)
#'
#' # Plot item-response profiles from the measurement model
#' plot(fit)
#'
#' @section Deprecated:
#' Deprecated as of tseLCA 2.0.0. It keeps working (and gives the same
#' estimates) but warns once per session; set
#' `options(tseLCA.warn.deprecated = FALSE)` to silence the warning. Use
#' [tseLCA()] or the step-wise functions:
#'
#' | `three_step()` | tseLCA 2.0 |
#' |---|---|
#' | `Y.names`, `n_classes` | `tse_lca(cbind(...) ~ 1, nclass = )` |
#' | `Zp.names` | `tse_covariate(, ~ ...)` or `tseLCA(... ~ covariates)` |
#' | `Zo.name`, `family` | `tse_distal(, outcome ~ 1, family = )`, or [tseLCA()] with the outcome after the bar |
#' | `step1` (measurement model from another sample) | `tse_classify(, newdata = )` |
#' | `startval` | `tse_lca(start = )` |
#' | `use.modal.assignment` | `tse_classify(assignment = )` |
#' | `use.bch` | `method = "BCH"` |
#' | `use.simple.cov` | `se = "robust"` |
#' | `rebase` | `ref` argument, or `relevel()` |
#' | `incomplete` | `tse_lca(missing = "fiml")` |
#' | `n_init`, `maxIter.measurement`, `measurement.tol`, `iter.measurement`, `R2.threshold`, `em.maxIter`, `covariate.tol`, `boundary.tol`, `correct.spec` | [tse_control()] |
#' | `get.twostep.vcov` | `tse_twostep(se = TRUE)` |
#'
#' @export
three_step <- function(
  data,
  Y.names,
  n_classes,
  Zp.names = NULL,
  Zo.name = NULL,
  step1 = NULL,
  startval = NULL,
  n_init = NULL,
  use.two.step = TRUE,
  use.modal.assignment = TRUE,
  include.intercept = TRUE,
  use.simple.cov = FALSE,
  incomplete = FALSE,
  boundary.tol = 1e-2,
  maxIter.measurement = 5000,
  measurement.tol = 1e-8,
  covariate.tol = 1e-6,
  iter.measurement = 10L,
  R2.threshold = 0.70,
  use.bch = FALSE,
  em.maxIter = 200L,
  get.twostep.vcov = FALSE,
  rebase = "C1",
  family = "gaussian",
  correct.spec = FALSE,
  verbose = FALSE
) {
  .tse_deprecated(
    "three_step()",
    "tseLCA() or the step-wise tse_lca(), tse_classify(), tse_covariate(), and tse_distal()"
  )
  n_step1_inputs <- sum(!is.null(step1), !is.null(startval), !is.null(n_init))
  if (n_step1_inputs > 1L) {
    stop(
      "`step1`, `startval`, and `n_init` are mutually exclusive ways of ",
      "controlling Step 1: supply a pre-fitted measurement model with ",
      "`step1`, a starting classification (or item-response probability ",
      "matrix) to fit one with `startval`, or a number of random restarts ",
      "with `n_init`, not more than one of these.",
      call. = FALSE
    )
  }

  opts <- list(
    use.two.step = use.two.step,
    use.modal.assignment = use.modal.assignment,
    include.intercept = include.intercept,
    use.simple.cov = use.simple.cov,
    incomplete = incomplete,
    boundary.tol = boundary.tol,
    maxIter.measurement = maxIter.measurement,
    measurement.tol = measurement.tol,
    covariate.tol = covariate.tol,
    iter.measurement = iter.measurement,
    R2.threshold = R2.threshold,
    use.bch = use.bch,
    em.maxIter = em.maxIter,
    get.twostep.vcov = get.twostep.vcov,
    rebase = rebase,
    correct.spec = correct.spec,
    startval = startval,
    n_init = n_init,
    verbose = verbose
  )
  ref_idx <- parse_rebase(rebase, n_classes)

  # Indicators as 0-based codes; a reused measurement model brings its own
  # categories.
  rec <- .recode_indicators(data, Y.names, .measurement_levels(step1))
  data <- rec$data

  # -- Step 1: measurement model ----------------------------------------------
  s1 <- .fit_step1(data, Y.names, n_classes, Zp.names, step1, ref_idx, opts,
                   Y.levels = rec$levels)
  dat <- .prepare_data(data, Y.names, Zp.names, Zo.name, family, opts, rec$levels)
  s1 <- .attach_step1_data(s1, dat, ref_idx, fitted_here = is.null(step1))

  if (is.null(Zp.names) && is.null(Zo.name)) {
    return(.new_measurement_fit(s1, dat, n_classes))
  }

  # -- Step 2: classification --------------------------------------------------
  s2 <- .step2(dat, s1$fit0, n_classes, opts)

  # -- Step 3: structural models -----------------------------------------------
  # Step-1 variance for the uncertainty correction (not needed for BCH or
  # robust-only standard errors)
  Sigma.1 <- if (use.simple.cov || use.bch) {
    NULL
  } else {
    .step1_varmat(s1, dat, ref_idx, boundary.tol)
  }

  cov <- if (!is.null(dat$Z_mat)) {
    .fit_covariate(dat, s1, s2, Sigma.1, n_classes, opts)
  }
  if (!is.null(cov)) {
    s1 <- cov$s1
  }
  dis <- if (!is.null(dat$Zo_mat)) {
    .fit_distal(dat, s1, s2, Sigma.1, cov, n_classes, family, opts)
  }

  .new_structural_fit(s1, s2, cov$fit, dis, n_classes, family, use.bch)
}

#' Measurement-only tseLCA object
#' @noRd
.new_measurement_fit <- function(s1, dat, n_classes) {
  posts <- step1_posteriors(dat$Y.obs, dat$mDesign, s1$fit0, dat$ivItemcat)
  structure(
    list(
      measurement_model = s1,
      llik = s1$fit0$LLKSeries[nrow(s1$fit0$LLKSeries)],
      AIC = s1$fit0$AIC,
      BIC = s1$fit0$BIC,
      R2entr = s1$fit0$R2entr,
      n_classes = n_classes,
      npar = (n_classes - 1L) + n_classes * sum(dat$ivItemcat - 1L),
      nobs = nrow(dat$Y.obs),
      posteriors = posts,
      classifications = max.col(posts)
    ),
    class = c("tseLCA_measurement", "tseLCA")
  )
}

#' Structural tseLCA object: covariate, distal, or both
#'
#' `cov_fit` is the tseLCA_covariate object from .fit_covariate() and `dis`
#' the distal component from .fit_distal(); either may be NULL.
#' @noRd
.new_structural_fit <- function(s1, s2, cov_fit, dis, n_classes, family, use.bch,
                                estimator = if (use.bch) "BCH" else "ML") {
  if (is.null(dis)) {
    return(cov_fit)
  }
  shared <- list(
    family = family,
    n_classes = n_classes,
    estimator = estimator,
    posteriors = s2$all$p.xy,
    classifications = max.col(s2$all$p.xy)
  )
  if (is.null(cov_fit)) {
    out <- c(dis, list(measurement_model = s1), shared)
    class(out) <- c("tseLCA_distal", "tseLCA_structural", "tseLCA")
    return(out)
  }
  out <- c(list(measurement_model = s1, covariate = cov_fit, distal = dis), shared)
  class(out) <- c("tseLCA_both", "tseLCA_structural", "tseLCA")
  out
}

#' Label of the Step-3 estimator
#' @noRd
.estimator_label <- function(opts) {
  if (isTRUE(opts$uncorrected)) "uncorrected" else if (opts$use.bch) "BCH" else "ML"
}
