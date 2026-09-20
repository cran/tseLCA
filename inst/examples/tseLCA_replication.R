#####################################################
### Replication Script for the tseLCA package:
### tseLCA: Three-Step Estimation for Latent Class Analysis
### -------------------------------------------------
### By: Sam Lee
### E-Mail: samlee@arizona.edu
#####################################################

# Run this script with:
# source(system.file("examples", "tseLCA_replication.R", package = "tseLCA"))

###################################################
### preliminaries
###################################################

# rm(list = ls())
# gc()
# r_opts <- options(
#   prompt = "R> ",
#   continue = "+  ",
#   width = 77,
#   digits = 4,
#   useFancyQuotes = FALSE,
#   warn = 1
# )

# Loading libraries and installing if unavailable

# Install the development version from GitHub
# if (!require("tseLCA")) {
#   if (!require("pak")) {
#     install.packages("pak")
#   }
#   pak::pak("SamLeeBYU/tseLCA")
# }

# #For access to the polytomous data 'election'
# if (!require("poLCA")) {
#   install.packages("poLCA")
# }

library(tseLCA)

###################################################
### generate synthetic data
###################################################
message(
  "\n---- Generating synthetic data ------------------------------------------------------------------------------------"
)

#Data-generating process described in

# Bakk, Z., & Kuha, J. (2018). Two-Step Estimation of Models Between Latent Classes and External Variables. Psychometrika, 83(4), 871–892. doi:10.1007/s11336-017-9592-7

#Data with six dichotomous outcome variables where
#probability of a Yes (1) response is 0.9 (separation="high")
#Generated with covariate, Zp~Unif{1,...5}, predicting latent
#indicator X with multinomial logit
d <- generate_data(
  n = 500,
  separation = "high",
  scenario = "covariate",
  seed = 1
)
head(d)

#Data generated with low separation:
#probability of a Yes (1) response is 0.7 (separation="low")
d.low <- generate_data(
  n = 500,
  separation = "low",
  scenario = "covariate",
  seed = 1
)
#Note that the covariate Zp and latent indicator X will
#be the same as 'd' since seed=1
head(d.low)

###################################################
### estimate measurement model
###################################################
message(
  "\n---- Step 1: Measurement model ------------------------------------------------------------------------------------"
)

#Calls multilev::multiLCA under the hood
#Uses k-means to initialize starting parameter estimates
d.measurement <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  measurement.tol = 1e-8
)
#S3 method summary for tseLCA objects
summary(d.measurement)

#In cases of low-separation, fit multiple measurement models
# and pick the best one
d.low.measurement <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  maxIter.measurement = 5000, #Run for a bit longer to get convergence
  iter.measurement = 10, #How many models we want to try
  R2.threshold = 0.9 #If R^2 entropy is below this value => trigger iter.measurement # of restarts
)
summary(d.low.measurement)

## tseLCA also adopts the S3 plot method, from multilevLCA, so you can call plot() on a tseLCA object
plot(d.measurement)

###################################################
### obtain two-step estimates
###################################################
message(
  "\n---- Two-step estimates (fitZ_from_fit0) ----------------------------------------------------------------"
)

#Two-step LCA model converges smoothly with EM algorithm
#These estimates can be used as starting parameter values
# for three-step procedure

d.fitZ <- fitZ_from_fit0(
  fit0 = d.measurement$measurement_model$fit0,
  data = d,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp"
)
#Reference class is C1 by default
#True slope parameters are -1 and 1 for C2 and C3, respectively
d.fitZ$mGamma

d.low.fitZ <- fitZ_from_fit0(
  fit0 = d.low.measurement$measurement_model$fit0,
  data = d.low,
  Y.names = paste0("Y", 1:6),
  Zp.names = "Zp",
  starting_val = d.fitZ$mGamma #Starting parameters are 0 by default
)
d.low.fitZ$mGamma

#####################################################
### obtain three-step estimates
#####################################################
message(
  "\n---- Three-step estimation --------------------------------------------------------------------------------------------"
)

#Three-step estimation can be done with one tseLCA function call
#This estimates the measurement model if one is not provided
# (two-step estimation for initial values for three-step estimation can be turned off with use.two.step=FALSE)
d.three_step <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp"
)
summary(d.three_step)
#Other useful S3 methods
coef(d.three_step)
vcov(d.three_step)
#Note that fitZ_from_fit0 does not return a tseLCA object so S3 methods will not work with those objects

#Default is use.modal.assignment=TRUE
#Recommended to change to proportional assignment in cases of poor separation
#(This ensures a proper Jacobian in the uncertainty propogation computation from the measurement model)
d.three_step.prop <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE
)
summary(d.three_step.prop)

#The argument use.simple.cov=TRUE will ignore uncertainty from the first stage and use robust third-stage SEs
#In cases of high separation, uncertainty from the first stage may be negligible
#Thus, this is more computationally efficient for large sample sizes
d.three_step.simple <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE
)
#Standard errors are smaller
summary(d.three_step.simple)

#By default, three_step() uses the ML correction proposed by Vermunt (2010)
#But the BCH approach proposed by Bolck, Croon, and Hagenaars (2004) can be used as well
d.three_step.bch <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.bch = TRUE
)
summary(d.three_step.bch)

#BCH works well in cases of high separation
#When separation is low, then the BCH weights may cause the hessian in the Newton-Rhapson algorithm to be ill-defined (not positive semi-definite)

#The BCH algorithm will run for a certain amount of iterations and if the Hessian is still not PSD, then the function will exit with an error

#This will take longer to run than usual because both the
#measurement model estimation (re-estimated here) and BCH estimation procedures struggle to converge
message(
  "  [Note: the following BCH call on low-separation data may take 1-2 minutes]"
)
bch.fail <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.bch = TRUE,
  maxIter.measurement = 5000,
  iter.measurement = 10
)
#Do not trust the estimates from the BCH model if the hessian is not proper

#Better to use the ML method (and proportional assignment) in cases of low-separation
d.low.three_step.prop <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE
)
#Much better than the BCH approach
summary(d.low.three_step.prop)

###############################################################
### choose different class as reference category in regression
###############################################################
message(
  "\n---- Reference class selection (rebase) ------------------------------------------------------------------"
)

#The default reference class for the multinomial logistic parameterization
# in the three-step/two-step estimation procedure is class one ("C1")
# however, this can be changed with the argument 'rebase'

summary(d.three_step.simple) #C1 is reference

d.three_step.simpleC2 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE,
  rebase = "C2"
)
summary(d.three_step.simpleC2)

d.three_step.simpleC3 <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.simple.cov = TRUE,
  rebase = "C3"
)
summary(d.three_step.simpleC3)


##############################################################
### pass in measurement model as argument
##############################################################
message(
  "\n---- Passing pre-fitted measurement model with step1 ------------------------------------------"
)

#The first stage measurement can be passed in as an argument
# even if computed on a different sample (uncertainty is properly calibrated)
#This is the biggest advantage of using tseLCA as opposed to one-step estimation (poLCA) or using mutilevLCA

d.three_step.prop2 <- three_step(
  data = d,
  #Y still needs to be in the data for covariate estimation for step 2
  #(which computes posterior probabilities for the current sample)
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  #Pass a measurement model here
  step1 = d.measurement$measurement_model
)
#Same as summary(d.three_step.prop)!
summary(d.three_step.prop2)

#A measurement model from a different sample with lower uncerainty (due to higher sample size)
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
#But covariates are only available in d.low
d.low.three_step.prop2 <- three_step(
  data = d.low,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  use.modal.assignment = FALSE,
  step1 = d.low.measurement2000$measurement_model,
  #You can also get the covariance matrix of the two-step esitmates returned by multiLCA (set FALSE by default)
  get.twostep.vcov = TRUE
)
summary(d.low.three_step.prop2)

#However, if we just want the two-stage estimates for starting values, we can call fitZ_from_fit0 and pass that into step1 as well
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
  #Now updated with fitZ
  step1 = d.low.measurement2000$measurement_model
)
summary(d.low.three_step.prop3)

##############################################################
### fixed Step-1 initialization (startval)
##############################################################
message(
  "\n---- Step 1 from a fixed starting classification (startval) --------------------------------"
)

#multiLCA's default k-means initialization is deterministic and can land on a
# local optimum; startval fixes the Step-1 EM start instead (disables k-means)

#startval can be a classification vector (one predicted class per row)...
startval.vec <- d.measurement$classifications
d.three_step.startval <- three_step(
  data = d,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  startval = startval.vec
)
summary(d.three_step.startval)

#...or a T x C item-response probability matrix (e.g. from an external solver
# such as poLCA, or, as here, a first-pass tseLCA fit's own mPhi expanded to
# both categories per item)
phi.compact <- d.measurement$measurement_model$fit0$mPhi #6 binary items x 3 classes
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
#Same solution as startval.vec above (both anchor Step 1 to the same fit)
summary(d.three_step.startval.phi)

##############################################################
### missing data
##############################################################
message(
  "\n---- Missing data handling --------------------------------------------------------------------------------------------"
)

#tseLCA handles missing data similar to multilevLCA

#Example
d.new <- generate_data(500, separation = "high", seed = 3)
sparsity <- 0.1
missing <- 1 -
  matrix(
    rbinom(prod(dim(d.new)), size = 1, prob = sparsity),
    nrow = dim(d.new)[1],
    ncol = dim(d.new)[2]
  )
missing[missing == 0] <- NA_real_
d.sparse <- d.new * missing
head(d.sparse)

#If incomplete = FALSE, the measurement model will drop *any* row in the Y matrix with a missing cell
d.sparse.measurement <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  incomplete = FALSE,
  verbose = TRUE
)
#This is how many we should drop
sum(apply(d.sparse[, 1:6], 1, \(x) any(is.na(x))))

summary(d.sparse.measurement)

#If incomplete = TRUE, we only drop rows where *all* of the data is missing
#We instead use FIML and retain incomplete item-responses
d.sparse.measurement2 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  incomplete = TRUE,
  verbose = TRUE
)

#Parameter estimates should be unbiased with data with MCAR-type missingness
summary(d.sparse.measurement2)

#Similar to multiLCA, we drop rows in the Z matrix if *any* covariate is missing
d.sparse.three_step <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = FALSE,
  verbose = TRUE
)
#In the above call we first drop rows where *any* Y is missing, then estimate the measurement model
#The Z that are dropped after that are the remaining rows where there's a missing Z (but not missing Y in that row)
sum(is.na(d.sparse$Zp[!apply(d.sparse[, 1:6], 1, \(x) any(is.na(x)))]))

#Similar, even if incomplete=TRUE, we still drop rows in Z with any missing covariate cell
# incomplete=TRUE only affects the way the measurement model is estimated
d.sparse.three_step2 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = TRUE,
  verbose = TRUE
)
#The above call drops rows with *all* missing data first and then proceeds to drop missing Z data
#Since there are no rows where all Y are missing, the number of Z dropped is simply the number of missing Zp here
sum(is.na(d.sparse$Zp))

summary(d.sparse.three_step2)

#Again any measurement model (with same structure) can be passed into the three.step function
d.sparse.three_step3 <- three_step(
  data = d.sparse,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zp.names = "Zp",
  incomplete = TRUE,
  verbose = TRUE,
  step1 = d.sparse.measurement2$measurement_model #Use measurement model from FIML
)
#Same as above
summary(d.sparse.three_step3)

##############################################################
### polytomous outcomes
##############################################################
message(
  "\n---- Polytomous items (election data) ----------------------------------------------------------------------"
)

#Like other LCA software, we can accomodate polytomous outcomes
#Here we demonstrate tseLCA on the election data, following the example from poLCA

#Data preparation
data(election, package = "poLCA")
elec <- election
elec.items <- colnames(election)[1:12]
#Like multiLCA, we require that all variables in Y are coded as sequential integers with base level coded as 0
elec[, elec.items] <- lapply(elec[, elec.items], function(x) as.integer(x) - 1L)
elec.measurement <- three_step(
  data = elec,
  Y.names = elec.items,
  n_classes = 3,
  incomplete = FALSE
)

elec.three_step <- three_step(
  data = elec,
  Y.names = elec.items,
  n_classes = 3,
  Zp.names = "PARTY",
  use.simple.cov = TRUE,
  step1 = elec.measurement$measurement_model,
  incomplete = FALSE,
  rebase = "C3"
)

party.x <- seq(from = 1, to = 7, length.out = 101)
pidmat <- cbind(1, party.x)
exb.tse <- exp(pidmat %*% coef(elec.three_step))
probs.tse <- (cbind(1, exb.tse)) / (1 + rowSums(exb.tse))

f.party <- cbind(
  MORALG,
  CARESG,
  KNOWG,
  LEADG,
  DISHONG,
  INTELG,
  MORALB,
  CARESB,
  KNOWB,
  LEADB,
  DISHONB,
  INTELB
) ~ PARTY
nes.party <- poLCA::poLCA(f.party, election, nclass = 3, verbose = FALSE)

exb.polca <- exp(pidmat %*% nes.party$coeff)
probs.polca <- (cbind(1, exb.polca)) / (1 + rowSums(exb.polca))

matplot(
  party.x,
  probs.tse,
  ylim = c(0, 1),
  type = "l",
  lwd = 3,
  col = 1,
  lty = 1,
  xlab = "Party ID: strong Democratic (1) to strong Republican (7)",
  ylab = "Probability of latent class membership",
  main = "Party ID as a predictor of candidate affinity class"
)
matlines(party.x, probs.polca, lwd = 2, col = 1, lty = 2)

text(3.9, 0.60, "Other")
text(6.2, 0.6, "Bush affinity")
text(2.0, 0.65, "Gore affinity")
legend(
  "topright",
  legend = c("tseLCA (3-step)", "poLCA (1-step)"),
  lty = c(1, 2),
  lwd = c(3, 2),
  col = 1,
  bty = "n"
)

###################################################
### distal outcomes
###################################################
message(
  "\n---- Distal outcomes --------------------------------------------------------------------------------------------------------"
)

#Unlike multilevLCA, we can accomodate distal outcomes

#DGP Reminder:
# Zp -> X -> Y (scenario="covariate")
# Zo <- X -> Y (scenario="distal")

d.distal <- generate_data(
  n = 500,
  separation = "high",
  scenario = "distal",
  seed = 4
)

#For estimation, however, we need a distributional assumption on P(Z|X)
#Current families available in tseLCA are:
# - family="gaussian": for continuous Zo
# - family="poisson": for positive discrete Zo
# - family="binomial": for binary Zo

d.distal.measurement <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3
)

#Both BCH and ML approaches are available here, too

d.distal.three_step.ml <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  step1 = d.distal.measurement$measurement_model,
  #Recommended but won't matter much here
  use.modal.assignment = FALSE,
  family = "gaussian" #default
)

#BCH with distal outcomes also has a propensity to fail here too due to problems with the hessian, as before
d.distal.three_step.bch <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zo",
  step1 = d.distal.measurement$measurement_model,
  #For comparison, although again, won't matter too much
  use.modal.assignment = FALSE,
  use.bch = TRUE,
  family = "gaussian"
)
#The advantage with BCH with distal outcomes is that there exists a closed-form solution when family="gaussian"
# for both parameter estimates and hessian (though again, not guarunteed to be PSD)

#Note that the mean (mu) parameters for each class were generated as -1, 1, and 0 for classes C1, C2, and C3, respectively
summary(d.distal.three_step.ml)
summary(d.distal.three_step.bch)

###################################################
### nominal categorical distal outcomes
###################################################
message(
  "\n---- Distal outcomes: nominal categorical (family=\"multinomial\") -------------------------------"
)

#For a Zo with 2+ categories and no natural ordering, family="multinomial"
# estimates a saturated model: the T x C matrix of class-conditional
# category probabilities pi_hat[t, c] = P(Zo = c | X = t)
cat.probs <- matrix(c(.7, .15, .15, .15, .7, .15, .15, .15, .7), 3, 3, byrow = TRUE)
d.distal$Zcat <- factor(apply(
  cat.probs[d.distal$X, ],
  1,
  \(p) sample(c("low", "mid", "high"), 1, prob = p)
))

d.distal.three_step.multi <- three_step(
  data = d.distal,
  Y.names = paste0("Y", 1:6),
  n_classes = 3,
  Zo.name = "Zcat",
  step1 = d.distal.measurement$measurement_model,
  family = "multinomial"
)
#coef() returns a T x C matrix of category probabilities (rows sum to 1)
coef(d.distal.three_step.multi)
summary(d.distal.three_step.multi)

###################################################
### omnibus test of class equality (distal outcomes)
###################################################
message(
  "\n---- Omnibus test: does the distal outcome differ across classes? -------------------------------"
)

#omnibus_test() runs a generalized Wald test of H0: the distal outcome's
# distribution is the same for every class -- works for any family
omnibus_test(d.distal.three_step.ml)
omnibus_test(d.distal.three_step.multi)

#########################################################################
### three-step estimation with both covariates (Zp) and distal outcomes
#######################################################################

message(
  "\n---- Covariate estimation and distal outcomes ------------------------------------------------------"
)

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

###################################################
### reset options
###################################################
message(
  "\n---- Done ------------------------------------------------------------------------------------------------------------------------------"
)
# options(r_opts)
# rm(list = ls())
# gc()

###################################################
### Print Session Information
###################################################

sessionInfo()
