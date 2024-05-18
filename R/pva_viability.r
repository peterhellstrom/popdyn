# Write function that calculates basic viability measures:
# Source: Dennis et al. 1991 (p. 124-130) and Morris et al 1999 (Handbook..., p. 19-21)
# Parameters:
# mu =
# sigma2 =
# N =
# Nc = Current population size
# Ne = Extinction threshold
#' @export
viability <- function(mu, sigma2, N, Nc, Ne) {

	t <- length(N) # length of time populations has been observed
	q <- length(diff(N)[!is.na(diff(N))]) # Number of transitions in data set

	# 1) Average population growth rate (r, lambda and CI's)
	r <- mu + (sigma2/2)
	# Confidence limits for r, Equation 68
	r.ci <- qnorm(0.025,0,1) * sqrt(sigma2*((1/t) + (sigma2/(2*q-1))))
	r.ll <- r + r.ci
	r.ul <- r - r.ci

	# 2) P for reaching a lower extinction threshold, Equation 84
	sigma2.hat <- (q-1)*(sigma2/q)
	pi.hat <- ifelse(mu<0, 1, (Ne/Nc)^(2*mu/sigma2.hat))

	# 3) The mean time to extinction, Equation 91
	xd <- log(Nc/Ne)
	theta <-  xd / abs(mu)
	var.theta <- (I(xd^2)*sigma2)/(I(mu^4)*t)
	theta.ci <- qnorm(0.025,0,1)*sqrt(var.theta)
	theta.ll <- ifelse(theta + theta.ci < 0, 0, theta + theta.ci)
	theta.ul <- theta - theta.ci

	out <- list(
		r = c(r = r, '2.5 %' = r.ll,'97.5 %' = r.ul),
		lambda = exp(c(lambda = r,'2.5 %' = r.ll,'97.5 %' = r.ul)),
		p.extinciton = pi.hat,
		extinction.time = c(theta = theta, '2.5 %' = theta.ll, '97.5 %'=theta.ul)
	)
	out
}
