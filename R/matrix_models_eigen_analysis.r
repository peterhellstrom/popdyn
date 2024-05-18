# EIGEN ANALYSIS OF A SQUARE MATRIX
#' @export
eigenall <- function(A) {

	lambda <- Re(eigen(A)$values[1]) # lambda, dominant eigenvalue

	# Stable stage distribution
	W <- eigen(A)$vectors # Right eigenvectors
	w <- W[,1] # Stable distribution
	w <- abs(w/sum(w)) # Rescale stable distribution to proportions

	# Reproductive value
	V <- Conj(solve(W)) # Left eigenvectors. Conj() returns complex conjugate, solve() solves a system of equations
						# Can also be extracted by using eigen-decomposition of the transposed matrix:
						# V <- eigen(t(A))$vectors; v <- abs(V[,1]/V[1,1])
	v <- abs(V[1,]/V[1,1]) # Reproductive values: rescale V relative to class 1

	# Sensitivity
	S <- v %o% w / sum(v*w) # %o% gives the outer product of v and w
	# Elasticity
	E <- (S*A) / lambda

	list(
		A = A,
		lambda=lambda,
		logLambda = log(lambda),
		eigenvals = Re(eigen(A)$values),
		w=w, v=v, S=S, E=E
	)
}

# The maximum proportional change in lambda that we can hope for by altering a vital rate, given the maximum possible rate (rmax) of that vital rate
#' @export
E.max <- function(A,i,j,rmax) {
	# From p. 343 in Morris & Doak, equation 9.13
	r <- A[i,j]
	E.r <- eigenall(A)$E[i,j]
	E.r*((rmax-r)/r)
}

# E.max(A=A, i=3, j=2, rmax=0.5)
