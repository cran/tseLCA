## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----setup, message = FALSE---------------------------------------------------
library(tseLCA)

## ----data---------------------------------------------------------------------
d <- generate_data(n = 1000, separation = "high", scenario = "covariate", seed = 1)
d$Zo <- draw_Zo(d$X, bk2018_params$distal_params) # add a distal outcome
head(d)

## ----enumeration--------------------------------------------------------------
f_items <- cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ 1
sel <- tse_lca(f_items, data = d, nclass = 1:4)
sel

## ----enumeration-plot, fig.height = 4-----------------------------------------
plot(sel)

## ----best---------------------------------------------------------------------
m <- best_model(sel, criterion = "BIC")

## ----measurement--------------------------------------------------------------
summary(m)

## ----measurement-plot---------------------------------------------------------
plot(m)

## ----measurement-accessors----------------------------------------------------
class_sizes(m)
item_probs(m)

## ----classify-----------------------------------------------------------------
cl <- tse_classify(m)
cl

## ----covariate----------------------------------------------------------------
fc <- tse_covariate(cl, ~ Zp)
summary(fc)

## ----covariate-tools----------------------------------------------------------
confint(fc)
AIC(fc)
anova(fc) # Wald test of each covariate term, across all classes

## ----covariate-predict--------------------------------------------------------
predict(fc, newdata = data.frame(Zp = 1:5))

## ----estimators---------------------------------------------------------------
fc_bch <- tse_covariate(cl, ~ Zp, method = "BCH")
fc_raw <- tse_covariate(cl, ~ Zp, method = "none")
ft <- tse_twostep(m, ~ Zp)
round(cbind(
  ML = coef(fc), BCH = coef(fc_bch), uncorrected = coef(fc_raw), two.step = coef(ft)
), 3)

## ----se-----------------------------------------------------------------------
round(cbind(
  corrected = sqrt(diag(vcov(fc))),
  robust = sqrt(diag(vcov(tse_covariate(cl, ~ Zp, se = "robust"))))
), 4)

## ----ref----------------------------------------------------------------------
coef(relevel(fc, ref = "C3"), matrix = TRUE)

d$group <- factor(ifelse(d$Zp > 3, "high", "low"))
anova(tse_covariate(cl, ~ Zp + group, data = d))

## ----distal-------------------------------------------------------------------
fd <- tse_distal(cl, Zo ~ 1)
summary(fd)

## ----omnibus------------------------------------------------------------------
omnibus_test(fd)

## ----multinomial--------------------------------------------------------------
d$Zcat <- cut(d$Zo, c(-Inf, -0.5, 0.5, Inf), labels = c("low", "mid", "high"))
fm <- tse_distal(cl, Zcat ~ 1, family = "multinomial", data = d)
round(coef(fm, matrix = TRUE), 3)
omnibus_test(fm)

## ----combined-----------------------------------------------------------------
fb <- tse_distal(fc, Zo ~ 1)
fb

## ----one-call-----------------------------------------------------------------
fit <- tseLCA(cbind(Y1, Y2, Y3, Y4, Y5, Y6) ~ Zp | Zo, data = d, nclass = 3)
all.equal(coef(fit), coef(fb))
round(classification(fit)$D, 3)

## ----multisample--------------------------------------------------------------
sub <- d[1:300, ]
fc_sub <- tse_covariate(tse_classify(m, newdata = sub), ~ Zp)
coef(fc_sub)
nobs(fc_sub)

## ----missing------------------------------------------------------------------
d_miss <- d
set.seed(2)
d_miss$Y1[sample(nrow(d), 100)] <- NA
d_miss$Y2 <- factor(d_miss$Y2, labels = c("no", "yes"))
m_fiml <- tse_lca(f_items, data = d_miss, nclass = 3, missing = "fiml")
nobs(m_fiml)
nobs(tse_lca(f_items, data = d_miss, nclass = 3))

## ----control------------------------------------------------------------------
tse_control(step1.maxit = 10000, n_init = 20)

