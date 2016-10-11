robustODE <- function(LS, init = initPars, crit.Tukey = 0.01, crit.huber = 0.05, optim.ctrl, resODE, maxit = 10) {
	environment(LS) <- environment()
	if (missing(resODE)) stop("Provide function to calculate residual")
	## =========================================================================
	## IRLS and robust
	## =========================================================================
	# 1. Obtain starting values. 
	## -------------------------------------------------------------------------
	# Use OLS to compute residuals and 
	message("Step 1. Estimating initial scale and weights")
	weighted <- FALSE
	fit <-  optim(init, LS, optim.ctrl)
	e <- resODE(fit)
	s <- mad(e)  # MAD scale median(abs(resid))/0.6745
	
	# 2. Huber estimation: 
	## -------------------------------------------------------------------------
	# Use Huber function with C = 1.345. 
	es <- abs(e) / s
	old.Huber <- w <- first.w <- getWeight(es, "Huber")
	# iteration until the maximum change in weights value is less than 0.05.
	hubercv <- FALSE
	weighted <- TRUE
	CheckConv <- function(old, current, crit) {
		change <- abs(current - old)  # absolute changes
		return(max(change) < crit)
	}
	message("Step 2. Starting Huber estimation...")
	iter <- 0
	while(!hubercv) {
		iter <- iter + 1
		message("H-iter: ", iter)
		current.fit <- optim(init, LS, optim.ctrl)
		current.e 	<- resODE(current.fit)
		current.s 	<- mad(current.e)  # MAD scale median(abs(resid))/0.6745
		current.es 	<- abs(current.e) / current.s
		current.Huber <- getWeight(current.es, "Huber")
		convergedH <- CheckConv(old.Huber, current.Huber, crit.huber)
		if (convergedH) {
			hubercv <- TRUE
			message("Huber weights converged!")
		} else if (iter == maxit) {
			hubercv <- TRUE
			message("Huber weights are not converged after ", maxit, " iterations!")
			message("Continuing Tukey step")
		} else {
			old.Huber <- w <- current.Huber
		}
	}
	# 3. Tukey estimation: 
	## -------------------------------------------------------------------------
	# Calculate case weight by using the biweight function with c = 4.685.
	tukeycv <- FALSE
	iter <- 0
	old.Tukey <- w <- getWeight(current.es, "biweight")
	message("Step 3. Starting Tukey estimation...")
	while(!tukeycv) {
		iter <- iter + 1
		message("T-iter: ", iter)
		current.fit <-  optim(init, LS, optim.ctrl)
		current.e 	<- resODE(current.fit)
		current.s 	<- mad(current.e)  # MAD scale median(abs(resid))/0.6745
		current.es 	<- abs(current.e) / current.s
		current.Tukey <- getWeight(current.es, "biweight")
		convergedT <- CheckConv(old.Tukey, current.Tukey, crit.Tukey)
		if (convergedT) {
			tukeycv <- TRUE
			message("Tukey weights converged!")
		} else if (iter == maxit) {
			tukeycv <- TRUE
			message("Tukey weights are not converged after ", maxit, " iterations!")
		} else {
			old.Tukey <- w <- current.Tukey
		}
	}
	message("Done")
	# return(list(L2 = fit$par, robust = current.fit$par) )
	return(list(L2 = fit$par, robust = current.fit$par, 
	       first.w = first.w, first.r = e, 
	       final.w = current.Tukey, final.r = current.e))
}
robustODE <- compiler::cmpfun(robustODE)
# 1. define RMS function (with weight options)
LS <- function(x) {
    pB      <- x[1]
    pD      <- x[2]
    pP      <- x[3]
    pC      <- x[4]
    paras   <- 10^c(pB = pB, pD = pD, pP = pP, pC = pC)
    tryCatch( 
    {
        out   <- deSolve::ode(state, times, UIV, paras)
        yhat  <- out[as.character(out[,1]) %in% as.character(cDa$time), "V"]
        # make sure viral load does not go negative
        if (any(yhat < 0)) { 
            message("Predict negative viral load, ABORT!")
            rms <-  1e+8
            return(rms)
        } else {
            # Number of replicates
            nm  <- rle(as.character(cDa$time))$lengths
            # Residuals
            x   <- cDa$V - rep(log10(yhat), times = nm)
            if (weighted) rms <- sqrt(mean( w * x^2, na.rm = TRUE))
            else rms <- sqrt(mean(x^2, na.rm = TRUE))
            return(rms)
        }
    }, error = function(e) {
            message(e)
            rms <- 1e+8
            return(rms)
        }
    )
}
# function to extract residuals
resODE <- function(fit) {
	names(fit$par) <- c("pB", "pD", "pP", "pC")
	out   <- deSolve::ode(state, times, UIV, 10^(fit$par))
	yhat  <- out[as.character(out[,1]) %in% as.character(cDa$time), "V"]
	nm  <- rle(as.character(cDa$time))$lengths
	x   <- cDa$V - rep(log10(yhat), times = nm)
	return(x)
}
# 1. Options for optim
initPars <- (lower + upper) / 2
optim.ctrl <- list(method = "L-BFGS-B", lower = lower, upper = upper, 
                   control = list(trace=1, factr = 1))  # hessian = TRUE?

# A skew data right
## Skew normal: Giam phuong sai, thay doi mean
## ----------------------------------------------------------------------------
## alpha <- 0
## mu <- 0
## omega <- 0.5
## delta <- alpha / sqrt(1 + alpha^2)
## mu + omega * delta * sqrt(2 / pi)
## ----------------------------------------------------------------------------
exMean <- function(omega = 0.5, alpha =-2) {
	delta <- alpha / sqrt(1 + alpha^2)
	return(omega * delta * sqrt(2 / pi))
}
## omega^2 * (1 - (2 * delta^2) / pi)
## ============================================================================

ChooseData2 <- function(timepoints, n.datapoints, std, skew = 0, undetect = FALSE) {
  ## ------------------------------------------------------
  ## std is the desire standard deviation
  ## extract only data in chosen time points
  ## Using as.character to couple the floating point issue when using %in%
  ## include undetectable data?
  ## ------------------------------------------------------
  sel <- Data0[which(as.character(Data0[, "time"]) %in% as.character(timepoints)), "V"]
  tmp <- NULL
  for (i in 1:length(sel)) {
    # tmp[[i]] <- log10(sel[i]) + rnorm(1000, 0, std) # Normal on log scale
    tmp[[i]] <- log10(sel[i]) + (sn::rsn(n=1000, xi=0, omega=std, alpha = skew) - exMean(std, skew))
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
