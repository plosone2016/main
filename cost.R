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
    }, error = function(e) 
    {
        message(e)
        mll <- 1e+8
        return(mll)
    })
}
W <- 1
RMSE <- function(x) {
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
            rms <- sqrt(mean(W * x^2, na.rm = TRUE))
            return(rms)
        }
    }, error = function(e) {
            message(e)
            rms <- 1e+8
            return(rms)
        }
    )
}
