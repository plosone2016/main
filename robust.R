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
