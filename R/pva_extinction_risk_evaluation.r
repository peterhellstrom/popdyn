# Used for setting pop sizes in simulated trajectories below a certain threshold to 0.
# x is a matrix with trajectories (rows = time, columns = realizations [replicates]) or
# a numeric vector with a single realization.
# Ne = extinction threshold
# IMPORTANT: in this function, populations cannot recover if they drop below Ne (they go extinct)
#' @export
find.first.0 <- function(x) {
  which(x==0)[1]
}

#' @export
eval.ext <- function(x, Ne=1) {

	# Evaluate if population size at time t is less than Ne
	ext.logic <- x < Ne
	# Convert logical matrix to binary (0=below Ne, 1=above Ne)
	ext.logic <- ifelse(ext.logic==TRUE,0,1)
	# Calculate a matrix with the cumulative product of each binary realization
	# This means that all (binary) numbers are set to 0 after quasi-extinction
	if (class(x) == "matrix") extinct <- apply(ext.logic, 2, cumprod)
	else if (class(x) == "numeric") extinct <- cumprod(ext.logic)
	else stop("Class not supported")

	# Create time variable
	time <- 0:(nrow(x)-1)

	# Calculate the number of extinct populations at each step (cumulative)
	N.G <- ncol(x) - rowSums(extinct)
	p.G <- N.G / ncol(x)

	# Calculate extinction times
	# Subtract 1, we start at t=0
	g <- apply(extinct,2, find.first.0) - 1
	N.g <- table(g)
	N.g <- merge(as.data.frame(time),data.frame(time=names(N.g), freq=as.numeric(N.g)),all.x=TRUE)

	# Matrix with counts (set extinct populations to zero)
	x <- x * extinct

	rownames(x) <- time
	colnames(x) <- 1:ncol(x)

	list(
		time = time,
		nG = N.G,
		pG = p.G,
		g = g,
		pg = N.g[,2] / ncol(x),
		x = x
	)
}

# Evaluate extinction risk, by comparing simulations with theoretical CDF
#' @export
eval.ext2 <- function(N0,mu,sigma2,Ne,tmax=50,nrep=1000,plot=TRUE) {

	Nc <- N0
	outmat <- replicate(nrep,geom.stoch(N0=N0,mu=mu,sigma2=sigma2,tmax=tmax,plot=FALSE))

	# Logically evaluate of population is at or below Ne-threshold
	eval.ext <- outmat <= Ne
	# Find extinction times (longest run of TRUE, all populations start as extant)
	# This actually gives t the year BEFORE the population hits Ne.
	# ext.t <- sapply(1:ncol(eval.ext), function(i) rle(eval.ext[,i])$lengths[1])
	# Number of time steps that population has spent at or below Ne
	ext.steps <- apply(eval.ext,2,cumsum)
	ext.mat <- ifelse(ext.steps>=1,NA,1)
	ext.pop <- outmat * ext.mat
	# Proportion extinct by future time t
	prop.ext <- 1 - (sapply(1:nrow(ext.pop), function(i) length(na.omit(ext.pop[i,]))) / nrep)

	# Extant populations
	# n.alive <- length(ext.t[ext.t==(tmax+1)])
	# Extinct populations
	# n.extinct <- length(ext.t[ext.t!=(tmax+1)])
	# Total number of populations
	# (nrep - n.alive) / nrep # Proportion extinct


	# Create plots
	if (plot==TRUE) {
		leg.text <- as.expression(bquote(paste(mu==.(round(mu,4)), ", ", sigma^2==.(round(sigma2,4)))))
		par(mfrow=c(2,2))

		ylimt <- quantile(outmat,c(0.01,0.99))
		plot(c(0,tmax),ylimt,xlab="Time",ylab="Population size",
			font.lab=2,type="n",main=paste(nrep, "simulations"),bty="l")
		title(sub=leg.text, cex=0.8)
		xv <- 0:tmax
			if (nrep > 100) {
				inds <- sample(size=100, x=1:nrep, replace=FALSE)
					for (i in 1:length(inds)) lines(xv,outmat[,inds[i]],lty=2)
					}
			if (nrep <= 100) for (i in 1:nrep) lines(xv,outmat[,i],lty=2)

		# Mean for each row
		lines(xv, rowMeans(outmat), lwd=2, col=2)

		hist(outmat[tmax+1,], col="lightgrey", breaks=30, xlim=quantile(outmat[tmax+1,],c(0,0.995)),
			xlab=paste("Population size at t =", tmax), ylab="Frequency", font.lab=2, main="")

		hist(log(outmat[tmax+1,]), col="lightgrey", breaks=30,
			xlab=paste("log(Population size) at t =", tmax), ylab="Frequency", font.lab=2, main="")
		abline(v=log(Ne), col=2, lty=2)

		ylim=c(0,xCDF(t=tmax,mu=mu,sigma2=sigma2,d=log(Nc/Ne)))
		plot(0:tmax,prop.ext,type="l",ylim=ylim,xlab="Time",ylab="Quasi-extinction prob. before time T",font.lab=2,bty="l",main="CDF")
		extCDF2(mu, sigma2, Nc=Nc, Ne=Ne, tmax=tmax, add=TRUE)

		par(mfrow=c(1,1))
	}

	invisible(outmat)
}
