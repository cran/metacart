## ----load----------------------------------------------------------------
library(metacart)
?SimData
summary(SimData)
set.seed(1)

## ----sim-----------------------------------------------------------------
res.simRE <- REmrt(formula = efk ~ m1 + m2 + m3 + m4 +m5, data = SimData, vi = vark, c = 0.5)
res.simRE

## ----sim2----------------------------------------------------------------
plot(res.simRE)
summary(res.simRE)


## ----FEsim---------------------------------------------------------------
res.simFE <- FEmrt(formula = efk ~ m1 + m2 + m3 + m4 +m5, data = SimData, vi = vark, c = 0.5)
res.simFE
plot(res.simFE)
summary(res.simFE)


## ----library-------------------------------------------------------------
data("dat.BCT2009")
summary(dat.BCT2009)

## ----RE1-----------------------------------------------------------------
set.seed(2017)
REres1 <- REmrt(formula = g ~ T1 + T2 + T4 + T25, vi = vi, data = dat.BCT2009, c = 1)


## ----summary-------------------------------------------------------------
summary(REres1)

## ----RE2-----------------------------------------------------------------
REres0 <- REmrt(formula = g ~ T1 + T2 + T4 + T25, vi = vi, data = dat.BCT2009, c = 0)
REres0
summary(REres0)
plot(REres0)


## ----sumRE2--------------------------------------------------------------
summary(REres1)

## ----FE------------------------------------------------------------------
FEres <- FEmrt(formula = g ~ T1 + T2 + T4 + T25, vi = vi, data = dat.BCT2009, c = 0.5)
FEres
summary(FEres)
plot(FEres)

