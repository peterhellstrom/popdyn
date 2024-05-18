# Investigate effects of first-order environmental autocorrelation on extinction risk
#' @export
stoch.growth.autocorr <- function(mu,sig2,rho,Nc,Ne,tmax) {

	N <- numeric(tmax+1)
	N[1] <- Nc

	beta <- sqrt(1-rho^2)
	sigma <- sqrt(sig2)

	# Define starting deviate - Morris & Doak used a random deviate sampled from N(mean=0,sd=1):
	eold <- rnorm(1,0,1)
	# Multiplying the random deviate with sigma*beta has a rather large effect on the outcome
	# Which method is preferable???
	# eold <- sigma*beta*rnorm(1,0,1)

	e <- numeric(tmax+1)
	e[1] <- eold

	for (t in 1:tmax) {
		enew <- (rho*eold) + (sigma*beta*rnorm(1,0,1))
		N[t+1] <- exp(mu + enew) * N[t]
		e[t+1] <- enew
		eold <- enew
		if (N[t+1] <= Ne) break
	}

	cbind(e,N)
}
