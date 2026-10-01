# tseLCA 2.0.0

## Step-wise interface

- New `tse_lca()` fits the Step-1 measurement model from a formula,
  `cbind(Y1, Y2, ...) ~ 1`, with a `control = tse_control()` argument.
  With several numbers of classes (e.g. `nclass = 1:6`) it returns a
  `tseLCA_select` class-enumeration table: log-likelihood, parameters, AIC,
  BIC, SABIC, entropy R^2, and smallest class. The table has `print()`,
  `plot()`, `as.data.frame()`, `best_model()`, and `[[` methods. The
  one-class (independence) model is fitted in closed form.
- New `tse_classify()` (Step 2) assigns observations to classes with a
  fixed measurement model (`assignment = "modal"` or `"proportional"`). It
  reports the classification-error probabilities P(W = s | X = t) and
  entropy R^2. With `newdata` it classifies another sample, replacing the
  `step1 =` route for using a measurement model from one sample on another.
  `tse_lca()` models now keep their `data` for this purpose.
- New Step-3 functions:
  - `tse_covariate()` fits the multinomial logit of class membership on a
    covariate formula (factors, interactions, transformations). Arguments:
    `method = "ML"`, `"BCH"`, or `"none"` (the uncorrected three-step
    estimator, for comparison); `se = "corrected"` or `"robust"`; `ref`
    (reference class); `start`. Methods: `predict()` (class probabilities
    given covariates), `relevel()`, and `anova()` (Wald tests per term).
  - `tse_distal()` fits a distal outcome (`Zo ~ 1`, any supported
    `family`, as a string or family object). Given a `tse_covariate()`
    model it fits the combined model.
  - `tse_twostep()` gives two-step estimates (Bakk & Kuha 2018). With
    `se = TRUE` they come with multilevLCA's corrected standard errors.
- New `tseLCA()` fits a whole model in one call:
  `cbind(indicators) ~ covariates | distal outcome` (Formula package). It
  chains the step-wise functions. Accessors `measurement()`,
  `classification()`, `covariate()`, and `distal()` extract the components
  of any fitted model.
- `tse_covariate()` and `tse_distal()` accept `data`: the classified data
  (same rows) with additional columns, e.g. variables created after
  classification.
- `omnibus_test()` returns a standard `htest` object.
- The vignette, README, package help page, and pkgdown reference index are
  rewritten for the step-wise interface. The vignette ends with a table
  mapping `three_step()` arguments to the new functions.
- `predict()` for measurement models gives posterior class probabilities
  (or modal classes) for new data, coded with the model's stored
  categories. `fitted()`, `formula()`, and `update()` also work.
- `item_probs(se = TRUE)` and `class_sizes(se = TRUE)` also return standard
  errors (delta method from the measurement model's variance).
- Fitted objects no longer keep multilevLCA's unused observation-level
  score matrix, which made up most of their size (about 80% smaller).
- `start` in `tse_lca()` and `tseLCA()` accepts `item_probs()` of a fitted
  model directly, e.g. to refit a chosen model with covariates.
- New `as_tse_lca()` builds a measurement model from given class sizes and
  item-response probabilities (e.g. estimates reported elsewhere, or saved
  from an earlier fit) without re-estimating it. The result can be passed to
  `tse_classify()` like any `tse_lca()` model.
- The replication and simulation scripts of the accompanying manuscript are
  no longer installed with the package; they are part of the manuscript's
  replication materials.

## Classes and methods

- Fitted objects now share a class hierarchy. `tseLCA` is the parent class of
  every fit; `tseLCA_structural` is the parent of `tseLCA_covariate`,
  `tseLCA_distal`, and `tseLCA_both`. Methods are defined once on the parent
  classes, not separately for each subclass.
- New methods on all fits: `logLik()` (with `df` and `nobs`, so `AIC()`
  and `BIC()` work), `nobs()`, `posterior()`, `classes()`,
  `class_sizes()`, and `item_probs()`.
- `summary()` now returns a `summary.tseLCA_structural` or
  `summary.tseLCA_measurement` object. Its tables print with
  `printCoefmat()` and are extracted with `coef(summary(fit))`.
- `plot()` now works on distal-outcome fits. These objects previously did
  not store the measurement model, so the plot failed.

## Data input

- Indicators can be factors, logicals, character, or numeric codes in any
  coding (e.g. 1..K). They are recoded to 0..K-1 internally, and their
  categories are stored with the measurement model
  (`$measurement_model$Y.levels`). Previously they had to be integers
  starting at 0. Codes such as 1..K were accepted but silently misread: the
  highest category was never matched.
- A measurement model reused through `step1` applies its own categories to
  new data. The new data may omit categories, and values outside them are
  an error.
- Covariate designs are built with `model.frame()`/`model.matrix()`, so
  factor covariates are dummy coded and unused factor levels are dropped.
  The intercept is named `"(Intercept)"` (was `"Intercept"`).
- Distal outcomes are validated for their family. Binomial outcomes may also
  be logical or two-category factors/characters.
- New `tse_control()` collects the numerical estimation settings.

## Bug fixes

- `tse_lca(missing = "fiml")` with random starts (`n_init`) failed with a
  cryptic error from multilevLCA ("sort_index(): detected NaN") when
  indicator values were missing: multilevLCA cannot start from a given
  classification on incomplete data. It now warns and uses multilevLCA's
  default initialization; a user-supplied `start` gives a clear error.
- Corrected (ML) standard errors: the Jacobian of the classification-error
  matrix with respect to the Step-1 parameters, which carries the Step-1
  uncertainty into Step 3, had its item-parameter columns ordered item by
  item, while the Step-1 variance is ordered class by class; with modal
  assignment it also differentiated the assignment as if it were the
  posterior probabilities. Both are fixed (checked against numerical
  derivatives). Point estimates, robust, and BCH standard errors are
  unchanged; corrected covariate standard errors typically increase by a
  few percent, and corrected standard errors no longer depend on the
  choice of reference class.
- Structural models with a reference class other than the first
  (`ref =`): `posterior()`, `classes()`, `item_probs()`, `class_sizes()`, and
  `plot()` reported the classes in the rebased order under the original
  labels, and in combined covariate + distal models the distal parameters
  of class t belonged to another class. All now use the measurement
  model's class order.
- Step-1 variance of polytomous indicators when the first (reference)
  category of an item has a boundary probability in some class: the free
  parameters of that item are log-ratios against this category, and their
  common shift, log P(Y = first | class), is not informed by the data. It is
  now treated as fixed, like other boundary parameters. Previously the
  variance of these log-ratios was huge (the information matrix was nearly
  singular) or all `NA` (singular), which made the corrected standard errors
  of the Step-3 models `NA`.
- If the Step-1 information matrix is still singular, Step-3 models warn
  and report robust standard errors (was `NA`).
- `tse_twostep(se = TRUE)` with a reference class other than the first
  returned the class-1-reference estimates and variance under the new class
  labels. They are now transformed to the requested reference class. It also
  no longer warns, wrongly, that multilevLCA's measurement model differs.
- For measurement-only fits, `posterior()` / `$posteriors` and `classes()` /
  `$classifications` were not in data-row order. They were read from
  multilevLCA's `fit0$mU`, which is sorted by response pattern. They are now
  computed from the Step-1 data in data-row order. Code that used
  `$classifications` from a measurement model, e.g. as `startval`, received
  mismatched rows in 1.1.x.
- Polytomous indicators were decoded from `fit0$mU` as category 0 whenever
  the measurement model was fitted with listwise deletion
  (`incomplete = FALSE`). The Step-1 information matrix was then singular, so
  `vcov()` of a polytomous measurement model and the Step-1-corrected Step-3
  variance (`use.simple.cov = FALSE`) were all `NA`. The measurement
  model now stores its Step-1 sample in data-row order and uses it for all
  of these computations. The fallback decoder for objects without stored data
  handles both of multilevLCA's codings.
- Parameter counts behind the reported AIC/BIC were wrong. Covariate models
  counted one-hot indicator columns, not free item parameters (e.g.
  40, not 22, for six binary items, three classes, one covariate).
  Combined covariate + distal models counted class sizes, not the
  covariate coefficients. Reported `AIC`/`BIC` for these models change;
  estimates and standard errors do not.

- Gaussian distal outcomes: the within-class variance was fixed at
  sigma2 = 1 in the Step-3 likelihood. Three-step ML now estimates a common
  sigma2 jointly with the class means, as in Bakk, Tekle & Vermunt (2013).
  Unlike ordinary regression, sigma2 does not factor out of the mean
  estimates here: it enters the posterior weights P(X = t | W, Zo). With the
  variance fixed at 1, ML class means were biased whenever the outcome's
  variance differed from 1 (e.g. roughly 3.5x too far apart for a residual SD
  of 3) and were not equivariant to rescaling the outcome. The Hessian,
  score, and Step-1/Step-2 uncertainty propagation now include sigma2. It is
  reported as `$sigma2` (estimate and, under ML, standard error), and counted
  in the log-likelihood's degrees of freedom. BCH class means were
  unaffected; their log-likelihood now uses the estimated variance.

- Three-step ML for distal outcomes with proportional assignment
  (`use.modal.assignment = FALSE`) maximized the wrong likelihood. It used
  log sum_t P(t) f(z|t) sum_s w_s P(W=s|t), with the assignment weights inside
  the log. The likelihood of Bakk, Tekle & Vermunt (2013), for the expanded
  data file weighted by the posterior assignment probabilities, is
  sum_s w_s log sum_t P(t) f(z|t) P(W=s|t). The covariate model already used
  this form. The E-step, score, Hessian/Jacobian, and Step-1/Step-2
  uncertainty propagation now use the correct likelihood for all distal
  families. Proportional-assignment ML distal estimates were biased away from
  zero (e.g. class means of -1.07 and 1.06, not -1 and 1, at mid
  separation). Modal-assignment fits, whose single assignment makes the two
  forms identical, BCH, and covariate models are unchanged.

- Covariate models with `include.intercept = FALSE` failed with a
  "length of 'dimnames'" error. The coefficient rows were always labeled
  with an `Intercept` row.

- In combined covariate + distal models, rows with a missing covariate but
  an observed distal outcome were kept in the distal model with `NA` class
  priors. The distal model now uses rows with complete covariates.

- When the Step-3 Hessian of an ML covariate model could not be inverted
  without an error being raised (non-finite entries), standard errors fell
  back to the outer product of the scores without a warning. Every fallback
  now warns. The documentation of `hessian = "opg"` (`correct.spec` in
  `three_step()`) described the outer-product estimator as model-robust; it
  is valid only under a correctly specified Step-3 model.

## Deprecated

- `three_step()` is deprecated in favor of `tseLCA()` and the step-wise
  functions. It keeps working, with the same estimates, and warns once per
  session. Its help page maps each argument to the 2.0 interface.
- `lca_step1()`, `lca_step1_startval()`, `fitZ_from_fit0()`, and
  `fitZ_from_multiLCA()` are deprecated in favor of `tse_lca()` and
  `tse_twostep()`. They warn once per session when called directly.
- `options(tseLCA.warn.deprecated = FALSE)` silences these warnings.

## Breaking changes

- `coef()` returns a named vector whose names match `vcov()`, so
  `confint()` works. `coef(fit, matrix = TRUE)` gives the previous
  Q x (T-1) (covariate) or T x C (multinomial distal) layout.
- `coef()` on a measurement model returns the log-ratio parameters that
  `vcov()` describes. The previous list is now available as
  `class_sizes()` and `item_probs()`.
- Covariate coefficient names use `"(Intercept)"` (was `"Intercept"`),
  e.g. `"(Intercept):C2"`.
- The `which` argument of `coef()` and `vcov()` is replaced by
  `component` (`"covariate"` or `"distal"`) and `step` (`"two_step"`).
  Passing `which` gives an error, so old code cannot silently return a
  different quantity.
- `vcov()` on a `tseLCA_both` object returns one matrix. The covariate
  and distal blocks are on the diagonal; the cross-covariances are not
  computed and are `NA`.

# tseLCA 1.1.1

- Corrected Jay Goodliffe's role in `Authors@R` from contributor (`"ctb"`) to
  author and copyright holder (`c("aut", "cph")`), reflecting his contribution
  to the package. No code changes.

# tseLCA 1.1.0

## Externally supplied Step-1 starting values

- Added `lca_step1_startval()`, a wrapper around `multilevLCA::multiLCA()`
  that injects a user-supplied Step-1 starting value and sets `kmea = FALSE`,
  bypassing `multilevLCA`'s default deterministic k-means-on-principal-
  components initialization. This addresses feedback that the default
  initialization can consistently converge to a local optimum of the Step-1
  log-likelihood on some datasets. `startval` accepts either:
  - an integer classification vector (`1..n_classes`, one entry per row of
    `data`), e.g. the modal class from an external solver run with many
    random starts ('StepMix', 'poLCA'); or
  - a numeric matrix of conditional item-response probabilities
    `P(Y_h = k | X = t)`, from which a classification is derived internally
    (naive-Bayes argmax under a flat class prior). This is the natural
    format for an externally estimated Step-1 solution that isn't tied to
    the current sample, e.g. `poLCA`'s `probs` output or a published
    item-response table.
  Either way, users can pass their external Step-1 solution straight through
  instead of hand-assembling a `tseLCA` measurement object.
- Added a matching `startval` argument to `lca_step1()` and
  `fitZ_from_multiLCA()`, so every place a measurement model is fit from raw
  data can use an external starting value instead of `multilevLCA`'s k-means
  initialization. Supplying `startval` skips the
  `iter.measurement`/`R2.threshold` random-restart logic, since restarting
  from fresh k-means seeds would defeat the purpose of a user-vetted start.
- Added a `startval` argument to `three_step()`, forwarded to `lca_step1()`
  (and to `fitZ_from_multiLCA()` when `get.twostep.vcov = TRUE`) so the
  measurement model can be fit from an external starting value in a single
  call.

## Multiple random-start Step-1 initialization (`n_init`)

- Added an `n_init` argument to `lca_step1()`, `fitZ_from_multiLCA()`, and
  `three_step()`: the unconditional multi-random-start analog of `n_init` in
  `StepMix` or `nrep` in `poLCA`. When supplied, the measurement model is
  fit `n_init` times from independent uniform-random classifications
  (`kmea = FALSE`, not `multilevLCA`'s deterministic k-means-on-PCA path),
  and the highest-log-likelihood fit is kept. This is a separate mechanism
  from the existing `iter.measurement`/`R2.threshold` restart logic, which
  reruns `multilevLCA`'s own k-means initialization and only when entropy
  R^2 is low; `n_init` restarts always run and never use k-means.
- `step1`, `startval`, and `n_init` are mutually exclusive ways of
  controlling Step 1 in `three_step()` (and `startval`/`n_init` in
  `lca_step1()` and `fitZ_from_multiLCA()`); supplying more than one errors.

## Bug fixes

- Fixed the orientation of the BCH weight matrix used throughout
  `use.bch = TRUE` estimation (both covariate and distal-outcome models).
  `pwx[s, t] = P(W = s | X = t)` is column-stochastic; the correct BCH
  weight matrix is `w.is %*% t(pwx)^-1` (Mplus Web Note 21: the row of
  `H^-1` for each case's most likely class, where `H = t(pwx)` is
  row-stochastic), not `w.is %*% pwx^-1`, which the code had been computing.
  The two orientations only agree when `pwx` is symmetric, so the bug was
  largely invisible in well-separated, balanced test cases; with genuine
  classification-error asymmetry it biased BCH point estimates and could
  produce the "negative column sums" error the package warns about (users
  were advised to fall back to `use.bch = FALSE`). With the corrected
  orientation, every row of the weight matrix sums to 1 and the class
  totals it implies exactly match the posterior class sizes. The
  weight-matrix computation is now consolidated into a single internal
  helper (`bch_weight_matrix()`) used by all four call sites that
  previously duplicated it, with a regression test asserting
  `rowSums(w.it) == 1` under a deliberately asymmetric classification-error
  matrix.

## Multinomial distal outcomes and an omnibus class-equality test

- Added `family = "multinomial"` to `three_step()`'s distal-outcome
  estimation, for a nominal categorical outcome with 2 or more categories
  (`Zo.name` may be a factor, character, or integer column). Estimates a
  saturated model -- the `T x C` matrix of class-conditional category
  probabilities `pi_hat[t, c] = P(Zo = c | X = t)` -- with the closed-form
  weighted-proportion estimator `pi_hat[t, c] = sum_i w_it * 1(y_i = c) /
  sum_i w_it`, for both `use.bch = TRUE` (BCH weights) and `use.bch =
  FALSE` (ML, through EM with the same closed-form M-step and
  responsibility-weighted E-step). This replaces the previous workaround of
  fitting one `family = "binomial"` model per category and renormalizing
  the resulting probabilities by hand, which wasn't constrained to the
  simplex before renormalizing and had no joint covariance across
  categories.
  - `coef()` returns the `T x C` probability matrix (rows sum to 1) instead
    of a length-`T` vector; `vcov()` returns its `(T*C) x (T*C)` sandwich
    covariance, which is necessarily rank-deficient (each class's row sums
    to 1).
  - For `use.bch = FALSE`, `use.simple.cov = FALSE` (the default) fully
    propagates both Step-1 measurement uncertainty and, when `Zp.names` is
    also supplied, Step-3 covariate uncertainty into the SEs, matching the
    existing gaussian/poisson/binomial ML paths exactly (the T x C
    generalization of each chain-rule term -- `C1_mat` for Step 1, `C_mat`
    for Step 2 -- only required expanding the "unit score" `g_it` across
    categories, since neither term's derivation otherwise depends on the
    outcome's dimensionality). The ML bread (`multinomial_ml_jacobian()`)
    is a closed-form Jacobian of the estimating equation, generalizing
    `ml_hessian_distal()`'s "observed = complete - missing information"
    correction to the T x C case. Both the bread and the full propagation
    (including the Step-2 term) were cross-validated against independent,
    from-scratch numerical differentiation (explicit per-case loops sharing
    no code with the package, `optim()` from multiple starting points for
    the point estimate) and matched to machine precision.
- Added `omnibus_test()`, a generalized Wald test of
  `H0: theta_1 = theta_2 = ... = theta_T` (the distal outcome's class-t
  parameter vector is the same for every class) for any `tseLCA_distal` or
  `tseLCA_both` object, regardless of family. Uses a Moore-Penrose
  pseudo-inverse of the contrast covariance so it remains valid when that
  covariance is singular, as it always is for `family = "multinomial"`; the
  resulting degrees of freedom recover the textbook `(T-1)*(C-1)` for a
  `T x C` chi-squared test of homogeneity in that case, and `T-1` for the
  scalar-parameter families.
- Unlike `family = "binomial"`, whose `coef()`/`vcov()` are on the logit
  scale, `family = "multinomial"` reports `coef()`/`vcov()` directly on the
  probability scale -- `Std.Error` is directly interpretable without a
  delta-method back-transform, but a symmetric interval
  `Estimate +/- 1.96*Std.Error` can fall outside the unit interval for a
  probability near a boundary, the same known limitation as a naive Wald
  interval for a sample proportion, and the per-cell `z.value`/`p.value`
  (testing each probability against 0) are rarely the question of interest.
  `print()`/`summary()` now print a one-line reminder of this after a
  multinomial distal-outcome table, pointing to `omnibus_test()` for the
  intended, boundary-safe test of whether the distribution differs across
  classes.

# tseLCA 1.0.0

-   Initial submission to CRAN.

## Core Estimation Framework

-   Implemented BCH and ML bias-adjusted three-step estimators for latent class analysis (LCA).
-   Added support for structural models containing covariates ($Z_p$), distal outcomes ($Z_o$), and combined models (estimating the relationship between $Z_p$ and the latent class first, followed by the distal outcome adjusting for covariate-adjusted posteriors).
-   Implemented analytic sandwich variance estimation to correctly propagate measurement uncertainty from the first-step LCA through classification-error correction in the final step.
-   Added a robust standard error option (`use.simple.cov = TRUE`) that bypasses the measurement-uncertainty correction for faster computation in large, well-separated samples.

## Measurement Model (Step 1) Integration

-   Integrated with the 'multilevLCA' package for efficient Step-1 measurement model estimation.
-   Added support for polytomous indicator items (0-based integer coding).
-   Implemented Full Information Maximum Likelihood (FIML) to handle missing data in the measurement model with the `incomplete = TRUE` argument (using a two-pass row-filtering strategy).
-   Added the ability to pass a pre-fitted measurement model (with the `step1` argument) to reuse across multiple structural models or apply to different sample subsets.
-   Implemented automated random restarts for the measurement model triggered when entropy $R^2$ falls below a user-specified threshold.

## Algorithmic Flexibility & Structural Models

-   Added support for both modal and proportional (soft) posterior class assignment (`use.modal.assignment`).
-   Integrated Gaussian, Poisson, and binomial families for distal outcome estimation.
-   Added the `rebase` argument to allow users to easily change the reference latent class for the multinomial logit parameterization while maintaining invariant log-likelihoods.
-   Implemented two-step EM estimation (`fitZ_from_fit0()`) to generate stable starting values for the three-step structural model.

## Utilities and Methods

-   Included standard S3 methods for `tseLCA` objects: `summary()`, `coef()`, `vcov()`, and `plot()` (which delegates to 'multilevLCA' for item-profile visualization).
-   Built a data-generating process (`generate_data()`) that replicates the Bakk & Kuha (2018) simulation study design for both covariates and distal outcomes under varying separation conditions.

# tseLCA 1.0.1

## CRAN resubmission
- Removed single quotes around acronyms in DESCRIPTION; added explanations
  of BCH, ML, and LCA.
- Replaced `T`/`F` with `TRUE`/`FALSE` throughout internal codebase
- Added `\value` tags to all exported functions missing them, including `bk2018_params`.
- `inst/examples`: examples now write to `tempdir()` instead of the home
  filespace.
- `inst/examples`: commented out `rm(list = ls())` calls.
- `inst/examples`: commented out `install.packages()` calls.

# tseLCA 1.0.2

## CRAN resubmission
- Title: corrected hyphenated compound capitalization to "Three-Step" per
  reviewer request.
