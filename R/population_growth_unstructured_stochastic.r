# Function for simulation of population growth in a stochastic environment:
#' @export
geom.stoch <- function(N0,mu,prob=NULL,sigma2=NULL,tmax,plot=TRUE,type="b") {
# NOTE: The for-loops could in theory be replaced with N0*cumprod(lambda)
	# Set up vectors that holds output
	N <- numeric(tmax + 1)
	N[1] <- N0 # Initiate the vector N, initial population size

	if (length(mu) > 1) {

		if (is.null(sigma2) != TRUE) stop("sigma2 is not allowed if you enter several mu's")

			if (is.null(prob) == TRUE) {
				lambda <- sample(x=mu, size=tmax, replace=TRUE)
			}

			if (is.null(prob) != TRUE) {
				if (length(mu) != 2) stop("Number of growth rates (mu) must be 2 if you enter a probability")
				# Binomial sampling of two growth rates
				lambda <- ifelse(rbinom(tmax,1,prob)==1,mu[1],mu[2])
			}

			lambda <- exp(lambda)

			for (t in 1:tmax) N[t+1] <- lambda[t] * N[t]
		}

	if (length(mu) == 1) {

		if (is.null(sigma2) == TRUE) {
			lambda <- exp(mu)
			for (t in 1:tmax) N[t+1] <- lambda* N[t]
		}

		if (is.null(sigma2) != TRUE) {
			lambda <- rlnorm(tmax,mu,sqrt(sigma2))
			for (t in 1:tmax) N[t+1] <- lambda[t]* N[t]
		}
	}

	if (plot == TRUE) {
		plot(0:tmax, N, type=type, pch=16, xlab="Time", ylab="Population size", font.lab=2, las=1)
	}

	N
}

#' @export
geom.stoch.n <- function(nsim=1000,N0,mu,prob=NULL,sigma2=NULL,tmax,plot=TRUE) {

	test.sim <- replicate(geom.stoch(N0,mu,prob,sigma2,tmax,plot=FALSE),n=nsim)
	lambda.a <- lognormal(mu, sigma2, "exp")[1]


	if (plot == TRUE) {
		ylims <- c(0, max(log(test.sim)))
		matplot(log(test.sim), type="l",ylim=ylims,xlab="Time", ylab="log(N)", font.lab=2, las=1,
			main="Density indenpendent growth in a random environment",bty="l")
		# Observed long-term (average) growth rate
		points(0:tmax, rowMeans(log(test.sim)), col=2, lwd=3, type="l")
		# Expected long-term growth rate if sigma2 = 0
		points(0:tmax, log(N0*lambda.a^(0:tmax)), type="l", lwd=2)
		#abline(a=log(N0),b=mu, lwd=4, col=4)
		#abline(a=log(N0),b=log(lambda.a), lwd=4, col=3)
		# Expected long-term growth rate in a stochastic environment
		points(0:tmax, log(N0) + mu*(0:tmax), type="l", lty=2, lwd=2)
		legend("topleft", c("No environmental variance", "Expected stochastic growth", "Observed stochastic growth"),
			lwd=c(2,2,3),col=c(1,1,2),lty=c(1,2,2), bty="n")
	}
}
