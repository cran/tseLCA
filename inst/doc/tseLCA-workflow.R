## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----setup, message = FALSE, warning=FALSE------------------------------------
library(tseLCA)

## ----generate-data------------------------------------------------------------
# High separation: P(Y_h = 1 | class) = 0.9 / 0.1
d <- generate_data(
  n = 500,
  separation = "high",
  scenario = "covariate",
  seed = 1
)
head(d)

## ----generate-data-low--------------------------------------------------------
# Low separation: P(Y_h = 1 | class) = 0.7 / 0.3
# Zp and X are identical to 'd' because seed = 1
d.low <- generate_data(
  n = 500,
  separation = "low",
  scenario = "covariate",
  seed = 1
)
head(d.low)

## ----measurement--------------------------------------------------------------
d.measurement <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  measurement.tol = 1e-8
)
summary(d.measurement)

## ----measurement-low----------------------------------------------------------
d.low.measurement <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  iter.measurement = 10,
  R2.threshold = 0.9
)
summary(d.low.measurement)

## ----plot-measurement, fig.width=6--------------------------------------------
plot(d.measurement)

## ----two-step-----------------------------------------------------------------
d.fitZ <- fitZ_from_fit0(
  fit0 = d.measurement$measurement_model$fit0,
  data = d,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp"
)
# True slopes: -1 (C2) and +1 (C3) relative to C1
d.fitZ$mGamma

## ----two-step-low-------------------------------------------------------------
d.low.fitZ <- fitZ_from_fit0(
  fit0 = d.low.measurement$measurement_model$fit0,
  data = d.low,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp",
  starting_val = d.fitZ$mGamma
)
d.low.fitZ$mGamma

## ----three-step-default-------------------------------------------------------
d.three_step <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp"
)
summary(d.three_step)

## ----s3-methods---------------------------------------------------------------
coef(d.three_step)
vcov(d.three_step)

## ----three-step-prop----------------------------------------------------------
d.three_step.prop <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE
)
summary(d.three_step.prop)

## ----three-step-simple--------------------------------------------------------
d.three_step.simple <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE
)
summary(d.three_step.simple)

## ----three-step-bch-----------------------------------------------------------
d.three_step.bch <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.bch = TRUE
)
summary(d.three_step.bch)

## ----bch-low, eval = FALSE----------------------------------------------------
# # Not run in vignette build (slow and and produces warnings)
# bch.fail <- three_step(
#   data = d.low,
#   Y.names = paste0("Y", 1:6),
#   n_classes = 3,
#   Zp.names = "Zp",
#   use.bch = TRUE,
#   maxIter.measurement = 2000,
#   iter.measurement = 10
# )

## ----three-step-low-----------------------------------------------------------
# Preferred approach for low separation
d.low.three_step.prop <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE
)
summary(d.low.three_step.prop)

## ----rebase-c1----------------------------------------------------------------
# Default: C1 as reference
summary(d.three_step.simple)

## ----rebase-c2----------------------------------------------------------------
d.three_step.simpleC2 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE,
  rebase = "C2"
)
summary(d.three_step.simpleC2)

## ----rebase-c3----------------------------------------------------------------
d.three_step.simpleC3 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE,
  rebase = "C3"
)
summary(d.three_step.simpleC3)

## ----step1-same-data----------------------------------------------------------
# Reuse the measurement model estimated above
d.three_step.prop2 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.measurement$measurement_model
)
summary(d.three_step.prop2)

## ----step1-different-sample---------------------------------------------------
# Measurement model from a larger low-separation sample
d.low2000 <- generate_data(
  n = 2000,
  separation = "low",
  scenario = "covariate",
  seed = 2
)
d.low.measurement2000 <- three_step(
  data = d.low2000,
  Y.names = paste0("Y", 1:6),
  n_classes = 3
)

# Apply to the smaller sample; get.twostep.vcov returns multilevLCA's
# bias-corrected vcov for the two-step estimates
d.low.three_step.prop2 <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.low.measurement2000$measurement_model,
  get.twostep.vcov = TRUE
)
summary(d.low.three_step.prop2)

## ----step1-inject-fitz--------------------------------------------------------
d.low.fitZ2 <- fitZ_from_fit0(
  fit0 = d.low.measurement2000$measurement_model$fit0,
  data = d.low,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp"
)
d.low.measurement2000$measurement_model$fitZ <- d.low.fitZ2

d.low.three_step.prop3 <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.low.measurement2000$measurement_model
)
summary(d.low.three_step.prop3)

## ----startval-vector----------------------------------------------------------
startval.vec <- d.measurement$classifications
d.three_step.startval <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  startval = startval.vec
)
summary(d.three_step.startval)

## ----startval-phi-------------------------------------------------------------
phi.compact <- d.measurement$measurement_model$fit0$mPhi # 6 binary items x 3 classes
phi.full <- do.call(
  rbind,
  lapply(1:6, \(h) rbind(1 - phi.compact[h, ], phi.compact[h, ]))
)
d.three_step.startval.phi <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  startval = phi.full
)
# Same solution as startval.vec above (both anchor Step 1 to the same fit)
summary(d.three_step.startval.phi)

## ----missing-data-setup-------------------------------------------------------
set.seed(42)
d.new <- generate_data(500, separation = "high", seed = 3)
sparsity <- 0.1
missing <- 1 -
  matrix(
    rbinom(prod(dim(d.new)), size = 1, prob = sparsity),
    nrow = nrow(d.new),
    ncol = ncol(d.new)
  )
missing[missing == 0] <- NA_real_
d.sparse <- d.new * missing
head(d.sparse)

## ----missing-listwise---------------------------------------------------------
d.sparse.measurement <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  incomplete = FALSE,
  verbose = TRUE
)
# Rows dropped = number of rows with at least one missing Y
sum(apply(d.sparse[, paste0("Y", 1:6)], 1, \(x) any(is.na(x))))
summary(d.sparse.measurement)

## ----missing-fiml-------------------------------------------------------------
d.sparse.measurement2 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  incomplete = TRUE,
  verbose = TRUE
)
summary(d.sparse.measurement2)

## ----missing-covariate--------------------------------------------------------
d.sparse.three_step <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = TRUE,
  verbose = TRUE
)
# Additional rows dropped from Step 3 due to missing Zp
sum(is.na(d.sparse$Zp))
summary(d.sparse.three_step)

## ----missing-reuse------------------------------------------------------------
d.sparse.three_step2 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = TRUE,
  step1 = d.sparse.measurement2$measurement_model,
  verbose = TRUE
)
summary(d.sparse.three_step2)

## ----polytomous-setup, message = FALSE----------------------------------------
data(election, package = "poLCA")
elec <- election
elec.items <- colnames(election)[1:12]

# Recode to 0-based integers as required by multilevLCA
elec[, elec.items] <- lapply(elec[, elec.items], \(x) as.integer(x) - 1L)

## ----polytomous-fit-----------------------------------------------------------
elec.measurement <- three_step(
  data = elec,
  Y.names = elec.items,
  n_classes = 3,
  #The poLCA example drops any row with a missing cell
  incomplete = FALSE
)

elec.three_step <- three_step(
  data = elec,
  Y.names = elec.items,
  n_classes = 3,
  Zp.names = c("PARTY"),
  step1 = elec.measurement$measurement_model,
  incomplete = FALSE,
  #With the neutral group as the base-category
  rebase = "C3"
)
summary(elec.three_step)

## ----elec-example-------------------------------------------------------------

party.x <- seq(from = 1, to = 7, length.out = 101)
pidmat <- cbind(1, party.x)
exb <- exp(pidmat %*% coef(elec.three_step))

matplot(
  party.x,
  (cbind(1, exb)) / (1 + rowSums(exb)),
  ylim = c(0, 1),
  type = "l",
  lwd = 3,
  col = 1,
  xlab = "Party ID: strong Democratic (1) to strong Republican (7)",
  ylab = "Probability of latent class membership",
  main = "Party ID as a predictor of candidate affinity class",
)
text(3.9, 0.60, "Other")
text(6.2, 0.6, "Bush affinity")
text(2.0, 0.65, "Gore affinity")


## ----distal-data--------------------------------------------------------------
d.distal <- generate_data(
  n = 500,
  separation = "high",
  scenario = "distal",
  seed = 4
)
# True class means: mu = (0, 1, -1) for C1, C2, C3

## ----distal-fit---------------------------------------------------------------
d.distal.measurement <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3
)

# ML estimator
d.distal.three_step.ml <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  step1 = d.distal.measurement$measurement_model,
  use.modal.assignment = FALSE,
  family = "gaussian"
)

# BCH estimator: closed-form M-step for distal outcomes
d.distal.three_step.bch <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  step1 = d.distal.measurement$measurement_model,
  use.modal.assignment = FALSE,
  use.bch = TRUE,
  family = "gaussian"
)

summary(d.distal.three_step.ml)
summary(d.distal.three_step.bch)

## ----multinomial-data---------------------------------------------------------
cat.probs <- matrix(c(.7, .15, .15, .15, .7, .15, .15, .15, .7), 3, 3, byrow = TRUE)
d.distal$Zcat <- factor(apply(
  cat.probs[d.distal$X, ],
  1,
  \(p) sample(c("low", "mid", "high"), 1, prob = p)
))

## ----multinomial-fit----------------------------------------------------------
d.distal.three_step.multi <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zcat",
  step1 = d.distal.measurement$measurement_model,
  family = "multinomial"
)
# coef() returns a T x C matrix of category probabilities (rows sum to 1)
coef(d.distal.three_step.multi)
summary(d.distal.three_step.multi)

## ----omnibus------------------------------------------------------------------
omnibus_test(d.distal.three_step.ml)
omnibus_test(d.distal.three_step.multi)

## ----covariate-distal-fit-----------------------------------------------------
d.covariate <- generate_data(
  n = 500,
  separation = "high",
  scenario = "covariate",
  seed = 4
)
d.covariate$Zo <- draw_Zo(d.covariate$X, bk2018_params$distal_params)
head(d.covariate)

d.covariate.three_step <- three_step(
  data = d.covariate,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  Zo.name = "Zo",
  use.modal.assignment = FALSE
)
summary(d.covariate.three_step)

## -----------------------------------------------------------------------------
three_step(
  data = d.covariate,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  use.modal.assignment = FALSE
) |>
  vcov() |>
  diag() |>
  sqrt()

## ----session-info-------------------------------------------------------------
sessionInfo()

