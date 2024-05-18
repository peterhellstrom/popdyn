# Here's a simple function that does the basic parts:
# matrices is a list with matrices to sample from
# maxt is the number of simulations
# w is a vector with proportions, rather than number (as we have used previously).
# Remember that w ideally should represent the stable stage distribution
#' @export
stoch.sim <- function(matrices,maxt,w) {
	if (class(matrices) != "list") stop("matrices must be a list")
	n <- length(matrices) # Number of matrices
	prob <- rep(1/n, n) # Probabilities of selecting each matrix in "matrices", equal probabilities
	maxt <- maxt # Number of iterations
	r <- numeric(maxt) # Vector that holds output
	w <- w # Should really be stable stage proportion
	n0 <- w # Initiate the population (with proportions!!!)

	for (t in 1:maxt) {
		ind <- sample(1:n, 1, prob = prob)
		A <- A.list[[ind]]
		n0 <- A %*% n0
		N <- sum(n0)
		r[t] <- log(N)
		n0 <- n0/N # Rescale to proportions!
	}
loglsim <- mean(r)
cat(paste("Stochastic lambda =", round(exp(loglsim),4),"\n"))
# exp(loglsim)
}

#' @export
proj.stoch <- function(n0, nsim=100, repl=1000, Alist, QET=500, plot=TRUE) {

	# Create mean matrix
	A <- Reduce("+",Alist)/length(Alist)
	# Set up matrix that stores final result of each replicate
	Nrepl <- matrix(nrow=repl, ncol=ncol(A))
	# Start outer loop, looping over replicates
	for (i in 1:repl) {
		# Create matrix that stores output for each time step
		N <- matrix(nrow=nsim+1, ncol=ncol(A))
		N[1,] <- n0 # Set first row to initial pop size
			# Start of inner loop, loop over time
			for (t in 1:nsim) {
				# Sample one matrix from list Alist
				Arand <- sample(Alist,1)[[1]]
				# Population projection
				N[t+1,] <- Arand %*% N[t,]
				# Final pop size
				Nfinal <- N[nrow(N),]
			} # End of inner loop
		# Send final pop size to the matrix Nrepl
		Nrepl[i,] <- Nfinal
	} # End of outer loop

	# Start calculations on the results from the simulations
	Nfinal <- rowSums(Nrepl) # Sum pop-size in final year

	# Calculate quasi-extinction risk:
	# The which() part evaluates and lists the position of the values in Nfinal
	# that are at or below the user-defined threshold QET
	# length() gives the number of values in a vector.
	# So the scalar pqext is the number of replicates that dropped below QET, divided by the number of replicates
	pqext <- length(which(Nfinal <= QET)) / repl

	# If plot==TRUE, draw histograms, with various extra details (such as lines at mean and median values).
	if (plot==TRUE) {
		par(mfrow=c(1,2))
		hist(Nfinal,breaks=20,col="lightgrey",
			xlab=paste("Final pop. size at t =",nsim),ylab="Frequency",
			font.lab=2,main=paste(repl, "simulations"))
		abline(v=mean(Nfinal),col=1,lty=2)
		abline(v=median(Nfinal),col=2,lty=3)
		legend("topright", c("Mean","Median"), col=c(1,2), lty=c(2,3), bty="n")

		# Also plot log of final population size for comparison
		hist(log(Nfinal),breaks=20,col="lightgrey",
			xlab=paste("log(Final pop. size at t) =",nsim),ylab="Frequency",
			font.lab=2,main=paste(repl, "simulations"))
		abline(v=mean(log(Nfinal)),col=1,lty=2)
		abline(v=median(log(Nfinal)),col=2,lty=3)
	par(mfrow=c(1,1))
	}
	# Create a list that holds all output
	out <- list(
		mean = mean(Nfinal),
		median = median(Nfinal),
		sd = sd(Nfinal),
		cv = sd(Nfinal)/mean(Nfinal),
		ci = quantile(Nfinal,c(0.025,0.5,0.975)),
		range = range(Nfinal),
		pqext = pqext
		)
	# Print list out as output
	out
}
