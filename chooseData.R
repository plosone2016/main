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
