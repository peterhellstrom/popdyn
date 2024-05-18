# markov.models
# Useful Markov chains (see Caswell p. 427-430)
# Mover-advancer model
# Mover-stayer model
# Random walk model

# N = number of states, q = probability of state change
#' @export
markov.models <- function(N,q,model="mover-stayer",plot=TRUE) {

	if (model == "mover-advancer") {
		P <- matrix(rep((1 - q)/N, N^2),ncol=N,nrow=N)
		P[row(P) == col(P)+1] <- q + ((1-q)/N)
		P[1,N] <- q + ((1-q)/N)
	}

	if (model == "mover-stayer") {
		P <- matrix(rep((1 - q)/N, N^2),ncol=N,nrow=N)
		diag(P) <- (q + (1 - q)/N)
	}

	if (model == "random-walk") {
		P <- matrix(rep(0, N^2),ncol=N,nrow=N)
		diag(P) <- q
		P[row(P) == col(P)+1] <- (1-q)/2
		P[row(P) == col(P)-1] <- (1-q)/2
		P[1,1] <- q + (1-q)/2
		P[N,N] <- q + (1-q)/2
	}
	P
}

# markov.2
# Construct transition matrix P for two-state Markovian environment
#' @export
markov.2 <- function(p,q) {
	P <- matrix(c(1-p,p,q,1-q), byrow=FALSE, nrow=2, ncol=2)
	P
}

# markov.plot
# z = a (time) series with states

# sample = TRUE/FALSE (if TRUE, take a sample of length samplen)
# start = start value of series z, and for sample
# samplen = if sample=TRUE, plot a sample of this length, with starting position given by start
#' @export
markov.plot <- function(z, start=0, sample=TRUE, samplen=100) {

	z.ts <- ts(z, start=start)

	if (sample == TRUE && length(z) >= samplen) {
		if (start > (length(z) - samplen)) stop("Invalid start position")
		z.ts <- ts(z.ts[start:(start+samplen)],start=start)
	}

	xv <- rep(time(z.ts),each=2)
	xv <- xv[-1]
	yv <- rep(z.ts, each=2)
	yv <- yv[-length(yv)]
	plot(xv,yv,type="l",xlab="Time",ylab="Environmental state",font.lab=2,main="")
}

# markov.sim
# Function that simulates transitions according to a Markov chain.
# Required input: a matrix P with transition probabilities, tmax = length of chain
# method = c("sample","caswell")
# Methods differ primarily in generation of inital environment
#' @export
markov.sim <- function(P, tmax=1000, method="sample") {
	if (unique(colSums(P)) != 1) stop("Columns do not sum to 1")

	if (method == "sample") {
		n <- nrow(P)
		# Initial environment (not entirely correct, approach compare with Caswell p. 420)
		x <- sample(1:n, 1) # Initiate chain by sampling a value between 1 and n (n being number of possible states)
		env <- numeric(tmax) # result is a vector that stores states
		env[1] <- x # Write start value to vector result
		# Sample a new value, dependent on previous state
		# Breakdown of code: result[i-1] returns the previous value (state) of the chain
		# P[result[i-1],] returns a column from the matrix P, and the values in this row are used as probabilities by function sample
			for (i in 2:tmax) {
				env[i] <- sample(1:n, 1, prob=P[,env[i-1]])
			}
	}

	if (method == "caswell") { # See Caswell 2001 p. 420
		P.dim <- nrow(P)
		# Generate initial environment
		eig <- eigen(P)
		d <- diag(eig$values, nrow=P.dim)
		w <- eig$vectors
		indmax <- which(d==max(abs(d)))
		w <- abs(w)
		w <- w[,indmax] / sum(w[,indmax])
		rand <- runif(n=1)
		x <- sum(rand >= cumsum(w))+1 # initial environment
		# Generate sequence of environments
		Pcum <- apply(P,2,cumsum)
		u <- runif(n=tmax)
		env <- numeric(tmax)
			for (t in 1:tmax) {
				x <- sum(u[t] >= Pcum[,x]) + 1
				env[t] <- x
			}
	}

	env
}

# Function for simulation of Markov environments:
#' @export
markov.sim2 <- function(P, len=1000) {
	if (unique(colSums(P)) != 1) stop("Rows do not sum to 1")
	n <- nrow(P)
	ini <- sample(1:n, 1) # Initiate chain by sampling a value between 1 and n (n being number of possible states)
	result <- numeric(len) # result is a vector that stores states
	result[1] <- ini # Write start value to vector result
	# Sample a new value, dependent on previous state
	# Breakdown of code: result[i-1] returns the previous value (state) of the chain
	# P[result[i-1],] returns a row from the matrix P, and the values in this row are used as probabilities by function sample
		for (i in 2:len) {
			result[i] <- sample(1:n, 1, prob=P[,result[i-1]]) }
	result
}
