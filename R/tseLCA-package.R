#' tseLCA: Three-Step Estimation for Latent Class Analysis
#'
#' @description
#' \pkg{tseLCA} relates latent classes to covariates and distal outcomes by
#' bias-adjusted three-step estimation. The latent class measurement model is
#' estimated first and held fixed, so the structural variables cannot change
#' the meaning of the classes; the structural estimates are corrected for the
#' classification error of the class assignments (BCH and ML estimators), and
#' their standard errors account for the uncertainty of the measurement model.
#' Measurement models are estimated with \pkg{multilevLCA} (Lyrvall et al.,
#' 2025).
#'
#' @section The three steps:
#' \enumerate{
#'   \item \strong{Measurement model} ([tse_lca()]): class sizes and
#'     class-conditional item-response probabilities, estimated from the
#'     indicators alone. With several numbers of classes, a class-enumeration
#'     table (AIC, BIC, SABIC, entropy) for choosing the number of classes.
#'   \item \strong{Classification} ([tse_classify()]): posterior class
#'     probabilities, modal or proportional class assignments, and the
#'     classification-error probabilities \eqn{P(W = s \mid X = t)}.
#'   \item \strong{Structural model}: a multinomial logit of class membership
#'     on covariates ([tse_covariate()]), and/or class-specific distributions
#'     of a distal outcome ([tse_distal()]), with the ML (Vermunt 2010; Bakk,
#'     Tekle & Vermunt 2013) or BCH (Bolck, Croon & Hagenaars 2004) correction.
#' }
#' [tseLCA()] runs all three steps from one formula,
#' `indicators ~ covariates | distal outcome`; [measurement()],
#' [classification()], [covariate()], and [distal()] extract the components of
#' a fitted model.
#'
#' @section Estimators and standard errors:
#' \describe{
#'   \item{`method = "ML"` (default)}{Vermunt's (2010) maximum likelihood
#'     correction, treating the assigned class as an indicator of the true
#'     class with known classification-error probabilities.}
#'   \item{`method = "BCH"`}{The Bolck-Croon-Hagenaars correction, reweighting
#'     the assignments by the inverse classification-error matrix. Reliable
#'     when classes are well separated.}
#'   \item{`method = "none"`}{The uncorrected three-step estimator, for
#'     comparison.}
#'   \item{`se = "corrected"` (default)}{Sandwich standard errors plus the
#'     propagated uncertainty of the Step-1 measurement model (Bakk, Oberski &
#'     Vermunt 2014), and, for combined models, of the covariate model.}
#'   \item{`se = "robust"`}{Sandwich standard errors of Step 3 only.}
#' }
#' The two-step estimator of Bakk & Kuha (2018) is available with
#' [tse_twostep()].
#'
#' @section Features:
#' \itemize{
#'   \item Binary and polytomous indicators, coded as factors, logicals,
#'     characters, or numbers; full-information maximum likelihood for missing
#'     indicator values (`missing = "fiml"`).
#'   \item Covariate formulas with factors, interactions, and transformations;
#'     Wald tests by term ([anova.tseLCA_covariate()]); predicted class
#'     probabilities ([predict.tseLCA_covariate()]); any reference class.
#'   \item Gaussian, Poisson, binomial, and multinomial distal outcomes, with an
#'     omnibus test of equality across classes ([omnibus_test()]).
#'   \item Measurement models estimated on one sample and applied to another
#'     ([tse_classify()] with `newdata`).
#'   \item Standard methods for fitted models: `print()`, `summary()`,
#'     `coef()`, `vcov()`, `confint()`, `logLik()`, `AIC()`, `BIC()`,
#'     `nobs()`, `predict()`, `plot()`, and `update()`.
#'   \item Simulation from the design of Bakk & Kuha (2018)
#'     ([generate_data()]).
#' }
#' The 1.x function [three_step()] is deprecated; its help page maps each of
#' its arguments to the current interface.
#'
#' @section Getting started:
#' ```r
#' vignette("tseLCA-workflow", package = "tseLCA")
#' ```
#'
#' @references
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
#' association between latent class membership and external variables using
#' bias-adjusted three-step approaches. \emph{Sociological Methodology},
#' 43(1), 272--311. \doi{10.1177/0081175012470644}
#'
#' Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
#' assignments to external variables: Standard errors for correct inference.
#' \emph{Political Analysis}, 22(4), 520--540.
#' \url{https://www.jstor.org/stable/24573086}
#'
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#' classes and external variables. \emph{Psychometrika}, 83(4), 871--892.
#' \doi{10.1007/s11336-017-9592-7}
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' \emph{Political Analysis}, 12(1), 3--27. \doi{10.1093/pan/mph001}
#'
#' Lyrvall, J., Di Mari, R., Bakk, Z., Oser, J., & Kuha, J. (2025).
#' Multilevel latent class analysis: State-of-the-art methodologies and their
#' implementation in the R package \pkg{multilevLCA}. \emph{Multivariate
#' Behavioral Research}, 60(4), 731--747. \doi{10.1080/00273171.2025.2473935}
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis}, 18(4), 450--469.
#' \doi{10.1093/pan/mpq025}
#'
#' @author Sam Lee \email{samlee@@arizona.edu}, Jay Goodliffe \email{goodliffe@@byu.edu}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
