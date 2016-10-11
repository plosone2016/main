# Supplemental codes

Examples

## Generate data
```R
sources("modelConditions.R")
tn1 <- c(1, 2, 3, 5, 7, 9)
set.seed(123)
cDa <- ChooseData(tn1, 5, 0.5)
```

## Run DE
```R
library(DEoptim)
sources("cost.R")
DEoptions <- DEoptim.control(parallelType = 1, packages = c("deSolve"), 
                parVar = c("UIV","state", "lower", "upper", "times", "cDa","W"), 
                trace = 1, itermax = 10000, steptol = 30, reltol = 1e-8, F = 0.8, CR = 0.9, NP = 50)
DEargs <- list(RMSE, lower, upper, DEoptions)
fitDE <- do.call("DEoptim", DEargs)
```
## Run Metropolis Hasting
```R
sources("MH.R")
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
sources("robust.R")
set.seed(123)
t24 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
data1000 <- replicate(1000, ChooseData2(smidR:::t24, 10, 0.5, skew = +10), FALSE)
bootRb <- function(iData) {
    environment(LS) <- environment(robustODE) <- environment(resODE) <- environment()
    cDa <- iData
    fit <- robustODE(LS = LS, resODE = resODE, optim.ctrl = optim.ctrl)
    return(fit)
}
outRb <- do.call(c, parallel::mclapply(data1000, bootRb, mc.cores = 4) )
```

The MIT License (MIT)
Copyright (c) 2016 Van Kinh Nguyen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
