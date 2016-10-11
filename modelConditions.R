UIV = function(t,state,parameters) {
  with(as.list(c(state,parameters)),{
    dU = - pBeta*U*V
    dI = pBeta*U*V - pDelta*I
    dV = pP*I - pC*V
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
