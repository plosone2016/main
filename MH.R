MetropolisAP <- function(
					   Y, 		# the data
                       tData, 	# corresponding time points 
                       Inits 	= smidR:::state,# Vector of initial values
                       model 	= smidR::UIV, 	# Model for deSolve
                       tODE 	= smidR:::times,# Time evaluating with deSolve
					   Sidm 	= c(3), # indices of the observed state in model
					   Sidd 	= c(1), # indices of the observed state in data
					   nPars 	= 5, 	# number of parameters
					   nStates 	= 3, 	# number of equations in the model
					   vPriors,			# require
					   uniSteps = 0.1, gui = FALSE, nPost = 1000, rMonitor = 500, 
					   nTunes = 10, n.chains = 3, 
					   Burnin = nTunes * rMonitor, thinning = 1, unitarget = 0.45, 
					   target = 0.35, largetarget = 0.234,
					   acceptTol = 0.075, optim = FALSE, Scale = 2.38)
{
	startt 		<- proc.time()

	tData <- as.character(tData)  # avoid rounding error
	nReps <- rle(tData)$length  # save time computing LL
	tData <- unique(tData)  

	minusLLL <- function(x) {
	    pB  <- x[1]
	    pD  <- x[2]
	    pP  <- x[3]
	    pC  <- x[4]
	    noise <- x[5]
	    paras <- c(pB = pB, pD = pD, pP = pP, pC = pC)
	    tryCatch(
	    {	
	        out   <- deSolve::ode(Inits, tODE, model, 10^(paras))
	        yhat  <- out[as.character(out[,1]) %in% tData, "V"]
	        if (any(yhat < 0)) { 
	            # message("Predict negative viral load, ABORT!")
	            mll <-  1e+8
	            return(mll)
	        } else {
	            x   <- Y - rep(log10(yhat), times = nReps)
	            ll  <- dnorm(x, mean = 0, sd = noise, log = TRUE)
	            mll <- -sum(ll)
	            return(mll)
	        }
	    }, error = function(e) {
	            message(e)
	            mll <- 1e+8
	            return(mll)
	        }
	    )
	}

	findLL <- function(params, noise) {
		tryCatch(
		{
	    	# Run ODE, get UIV and sensitivities
	    	XData <- deSolve::ode(Inits, tODE, model, 10^(params))  # par in raw scale

	    	XData <- XData[as.character(XData[,1]) %in% tData, "V"]

		    Xpred <- rep(XData, times = nReps)
		    error <- Y - log10(Xpred)

		    # Calculate the log-likelihoods of the parameters
		    # # Calculate the current likelihoods of the current X estimates
	        # Noting that the estimates are in raw scale
		    tempLL <- sum(dnorm(error, mean = 0, sd = noise, log=TRUE))
	 		return(tempLL)
	    }, warning = function(warn) {
	    	# message(warn)
	    	tempLL <- -1e300
	 		return(tempLL)
	    }, error=function(e) {
	    	message(e)
	    	tempLL <- -1e300
	 		return(tempLL)
		})
	}

	naiveMu 	<- vPriors$firstM[-5]
	naiveSigma 	<- InitSigma()
	naivePars 	<- MASS::mvrnorm(n = 1, naiveMu, naiveSigma)
	naiveNoise 	<- Priors('random', vPriors$firstM[5], vPriors$secondM[5], 
	                     vPriors$family[5])

	lower <- c(qnorm(1e-4, vPriors$firstM[1:4], vPriors$secondM[1:4]), 
	          qunif(1e-4, vPriors$firstM[5], vPriors$secondM[5]))
	upper <- c(qnorm(1-1e-4, vPriors$firstM[1:4], vPriors$secondM[1:4]), 
	          qunif(1-1e-4, vPriors$firstM[5], vPriors$secondM[5]))

	fit <- NULL
	if (optim) {
		message("\nEstimating asymtotic posterior covariance matrix...")
		# Generating starting values for finding covariance matrix
		Starts 	<- c(naivePars, naiveNoise)
		fit 	<- tryCatch({
			fit <- optim(Starts, minusLLL, method = "L-BFGS-B", 
			  lower = lower, upper = upper, control = list(trace=1), hessian = TRUE)
		}, error = function(e) {
			message("Failed to estimates the posterior asymtotic covariance matrix,
			        use naive sigma")
			message(e)
			return(NULL)
		})
	}
	# This part can add option robust to make positive definite
	# and handling error of not able to inverse
	if (!is.null(fit)) {
		fisher_info <- solve(fit$hessian[1:4, 1:4])
		prop_sigma 	<- sqrt(diag(fisher_info))
		Sigma 		<- diag(prop_sigma) * (Scale^2) / 4
	} else {
		Sigma <- InitSigma(std = 0.5)
	}

	# print(Sigma)

	message("\nInitialisation Completed. Computing chains...")

	# =====================================================================
	chainFn <- function(...) 
	{
	# Still starting at random, but propose using the estimated covariance
	if ( optim & (!is.null(fit)) ) {
		Pars 	<- MASS::mvrnorm(n = 1, fit$par[1:4], Sigma)
		Noise 	<- fit$par[5] + rnorm(1) * uniSteps
	} else {
		Pars 	<- MASS::mvrnorm(n = 1, naiveMu, Sigma)
		Noise 	<- naiveNoise + rnorm(1) * uniSteps
	}
	
	Pars 		<- setNames(Pars, c("pB", "pD", "pP", "pC"))
	CurrentLL 	<- findLL(Pars, Noise)

	# Set up proposal counters
	Accepted 	<- rep(0, 2)
	Mutation 	<- rep(0, 2)

	# Set up parameter history variable
	ParaHistory <- matrix(0, nPost, nPars)
	LLHistory   <- vector('numeric', nPost) #only one state

	# Store the burnIn for adjusting JUMP
	BurnResults <- matrix(0, Burnin, nPars)

	# Initialise iteration number
	nIter 		<- 0
	Continue 	<- TRUE
	Converged   <- FALSE

	# Main loop
	while (Continue) {

	    nIter = nIter + 1
	    if(nIter %% rMonitor == 0) cat('iter ', nIter, '\n')
		NewParas 	<- MASS::mvrnorm(n = 1, Pars, Sigma)
		Mutation[1] <- Mutation[1] + 1

		# checking boundary
		isInBound 	<- all(NewParas[-5] < upper[-5] & NewParas[-5] > lower[-5])
		if (isInBound == FALSE | any(is.na(NewParas)) ) {
			badMove <- TRUE 
		} else {
			badMove <- FALSE
		}

		# JumpBack 	<- Jump(Pars, NewParas, Sigma, Noise) # the same
		# JumpNext	<- Jump(NewParas, Pars, Sigma, Noise) # the same
		
		isAccept 	<- FALSE
		
		if (badMove == FALSE) {
			ProposedLL 	<- findLL(NewParas, Noise)
			Ratio 		<- ProposedLL - CurrentLL
			# Ratio 	<- ProposedLL - CurrentLL + JumpBack - JumpNext 
		    isAccept 	<- !is.na(Ratio) & ( Ratio > 0 | (Ratio > log(runif(1))) )
		}

	    if (isAccept) {
            Pars      	<- NewParas
            CurrentLL 	<- ProposedLL
            Accepted[1] <- Accepted[1] + 1
	    }

	    # Now mutating the noise given the new (if updated)/current par
	    newNoise 	<- Noise + rnorm(1) * uniSteps
	    Mutation[2] <- Mutation[2] + 1
	    ProposedLL 	<- findLL(Pars, newNoise)

	    isAccept 	<- FALSE
	    Ratio 		<- ProposedLL - CurrentLL
	    isAccept 	<- !is.na(Ratio) & ( Ratio > 0 | (Ratio > log(runif(1))) )
	    
	    if (isAccept) {
            Noise 		<- newNoise
            CurrentLL 	<- ProposedLL
            Accepted[2] <- Accepted[2] + 1
	    }

	    # Save parameters if converged and according to thinning rate
	    if (Converged & (nIter %% thinning)==0 ) {
	        ParaHistory[(nIter-nIterConverg)/thinning, 1:4] <- Pars
	        ParaHistory[(nIter-nIterConverg)/thinning, 5] 	<- Noise
	        LLHistory[(nIter-nIterConverg)/thinning] 		<- CurrentLL
	    }
	    
	    # If not yet converged
	    if (!Converged) {
	    	
		    # Save parameters in burn in for adjusting
	        BurnResults[nIter, 1:4] <- Pars
	        BurnResults[nIter, 5] 	<- Noise

			# Adjust parameter proposal widths
	        if (nIter > rMonitor & (nIter %% rMonitor) == 1) {

	            message("Iteration", nIter)
	            message(rep('=', 80))

	            arate 	<- Accepted/Mutation
	            cat(100*arate[1], '% mutation acceptance for parameter ', 1:4, '\n')
	            cat(100*arate[2], '% mutation acceptance for parameter ', 5, '\n')

	            # Adjust JUMP if out range
	            outBound <- (arate[1] < (target - acceptTol) | arate[1] > (target + acceptTol))
	            if (outBound) {
		            # Update default scale
		            Scale <- tuneScale(Scale, arate[1], targetRate = target)
	                # update cov
	                region 	<- (nIter - rMonitor):(nIter - 1)
	                covv 	<- cov(BurnResults[region, 1:4])
					Sigma 	<- tuneCov(covv, Sigma) * (Scale^2) /4
	            }
            	# Adjust for sd
                if (arate[2] < (unitarget - acceptTol)) {
                    uniSteps <- uniSteps * 0.9
                } else if (arate[2] > (unitarget + acceptTol)) {
                    uniSteps <- uniSteps * 1.1
                }

	            message('Current estimates: ')
	            cat(round(BurnResults[nIter,], 5), "\n")
	            # Reset counters
	            Accepted <- rep(0, 2)
	            Mutation <- rep(0, 2)
	        }
	        if (nIter >= Burnin) {
	            Converged    <- TRUE
	            nIterConverg <- nIter
	            message("\n", rep('=', 80))
	            message('Stop tuning at iteration number ', nIterConverg, "\n")
	            cat('Accepted rate ', arate, "\n")
	            message(rep('=', 80))
	        }
	    } else if (nIter == nIterConverg + nPost * thinning) {
             Continue <- FALSE
	    }
	}
	colnames(ParaHistory) <- colnames(BurnResults) <- c("pB", "pD", "pP", "pC","sd")
	attr(ParaHistory, "mcpar") <- c(thinning, nPost*thinning, thinning)
	attr(BurnResults, "mcpar") <- c(1, Burnin, thinning)
	attr(ParaHistory, "class") <- attr(BurnResults, "class") <- "mcmc"
	return(list(ParaHistory = ParaHistory, LLHistory = LLHistory, Burn = BurnResults))
	} 
	#end chains
	# =====================================================================

	OutParallel <- do.call( c, parallel::mclapply(seq_len(n.chains), chainFn, 
	                       mc.cores = n.chains, mc.silent = FALSE))
	
	MCMClist 	<- OutParallel[c(1, 4, 7)]
	LLlist 		<- OutParallel[c(2, 5, 8)]
	Burnlist 	<- OutParallel[c(3, 6, 9)]

	class(MCMClist) <- class(Burnlist) <- "mcmc.list"

	endt			<- proc.time() - startt

	return(list(MCMC = MCMClist, Burn = Burnlist, LL = LLlist, runtime = endt))
}
