# Function that simulates demography in a Markov environment
# A = list with Leslie matrices
# P = transition Matrix for Markov chain
# w = initial PROPORTIONS (stable stage distribution)
# n = length of Markov chain
#' @export
markov.sim.stoch <- function(A,P,w=NULL,n) {

	if (class(A) != "list") stop("matrices must be a list")
	if (length(A) != ncol(P)) stop("Number of matrices in A must match number of states in P")

	# Create a sequence of "environments", Markov chain
	z <- markov.sim(P,len=n)

	# Generate sequence of matrices
	A.markov <- lapply(1:n, function(i) A[[z[i]]])

	if (is.null(w) == FALSE) w <- w # stable stage distribution

	if (is.null(w) == TRUE) {
		# Calculate mean matrix and stable stage distribution for the mean matrix
		A.mean <- Reduce("+", A.markov)/length(A.markov)
		w <- eigenall(A.mean)$w
	}

	r <- numeric(n)
	n0 <- w

	for (t in 1:n) {
		A.tmp <- A.markov[[t]]
		n0 <- A.tmp %*% n0
		N <- sum(n0)
		r[t] <- log(N)
		n0 <- n0/N # Rescale to proportions to avoid numerical constraints
	}
	# Calculate mean of r
	loglsim <- mean(r)
	# Create list with output
	out <- list(
		z = z,
		A.markov = A.markov,
		r = r, # vector with growth rates
		loglsim = loglsim, # r = log(lambda)
		lsim = exp(loglsim) # lambda
	)

	cat(paste("Stochastic lambda =", round(exp(loglsim),4),"\n"))
	invisible(out)
}
