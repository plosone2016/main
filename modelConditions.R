Threestates = function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    dU = - pBeta*U*V
    dI = pBeta*U*V - pDelta*I
    dV = pP*I - pC*V
    list(c(dU,dI,dV))
  })
}
UIV <- function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    dU = - pB*U*V
    dI =   pB*U*V - pD*I
    dV =   pP*I   - pC*V
    list(c(dU,dI,dV))
  })
}
times   <- seq(0, 12, by = 0.01)
state   <- c(U = 10^6, I = 0, V = 10)
paras   <- c(pBeta = 1e-05, pDelta = 1.6, pP = 5, pC = 3.7)
lower   = log10(c(1e-7, 1e-2, 1e+0, 1e-1))
upper   = log10(c(1e-3, 1e+2, 2e+2, 1e+2))

# The "correct" data
Data0   <- ode(state, times, Threestates, paras)
Data0   <- as.data.frame(Data0)
## ------------------------------------------------------------------------
## Data generate
## ------------------------------------------------------------------------
## Generate using adjusted parameter set
## Adding noise
## Randomly select time points and number of data points
## Adding noise as random, assume log-normal distribution of
## viral load measurements
## ------------------------------------------------------------------------
ChooseData <- function(timepoints, n.datapoints, std, undetect = FALSE) {
  ## ------------------------------------------------------
  ## std is the desire standard deviation
  ## extract only data in chosen time points
  ## Using as.character to couple the floating point issue when using %in%
  ## include undetectable data?
  ## ------------------------------------------------------
  sel <- Data0[which(as.character(Data0[, "time"]) %in% as.character(timepoints)), "V"]
  tmp <- NULL
  for (i in 1:length(sel)) {
    tmp[[i]] <- log10(sel[i]) + rnorm(1000, 0, std) # Normal on log scale
    tmp[[i]] <- sample(tmp[[i]], n.datapoints, replace = TRUE)
  }
  # export from list to data frame and give column names
  tmp <- cbind(rep(timepoints, each = n.datapoints), unlist(tmp))
  tmp <- data.frame(tmp)
  colnames(tmp) <- c("time","V")
  # remove Undetectable
  if (!undetect) tmp <- tmp[tmp[, "V"] >= log10(50), ]
  return(logV = tmp)
}
