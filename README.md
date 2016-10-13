
# Supplemental codes

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

## Generate data
```R
library(deSolve)
source("modelConditions.R")
tn1 <- c(1, 2, 3, 5, 7, 9)
set.seed(123)
cDa <- ChooseData(tn1, 5, 0.5)
```

## Run DE
```R
library(DEoptim)
source("cost.R")
DEoptions <- DEoptim.control(parallelType = 1, packages = c("deSolve"), 
                parVar = c("UIV","state", "lower", "upper", "times", "cDa","W"), 
                trace = 1, itermax = 10000, steptol = 30, reltol = 1e-8, F = 0.8, CR = 0.9, NP = 50)
DEargs <- list(RMSE, lower, upper, DEoptions)
fitDE <- do.call("DEoptim", DEargs)
```
## Run Metropolis Hasting
```R
source("MH.R")
Y <- t(cDa[,-1  ])
TimePoints <- t(cDa[,-2])
# List of two parameter distribution (normal or gamma)
nPriors <- data.frame(firstM = c(-5, 0, 1, 0.5, 0), #first moment
                      secondM = c(0.5, 0.5, 0.5, 0.5, 1), #second moment
                      family = c('norm','norm','norm','norm','uni'), 
                      stringsAsFactors = FALSE) #family
realAP  <- fastMH(Y, tData = TimePoints, vPriors = nPriors, nTunes = 10, 
                 rMonitor = 500, nPost = 1000, thinning = 10, optim = TRUE)
```

## Run robust
```R
library("sn")
source("robust.R")
set.seed(123)
t24 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
data1000 <- replicate(1000, ChooseData2(t24, 10, 0.5, skew = +10), FALSE)
bootRb <- function(iData) {
    environment(LS) <- environment(robustODE) <- environment(resODE) <- environment()
    cDa <- iData
    fit <- robustODE(LS = LS, resODE = resODE, optim.ctrl = optim.ctrl)
    return(fit)
}
outRb <- do.call(c, parallel::mclapply(data1000, bootRb, mc.cores = 4) )
```
