# Calculate stage transitions, assuming fixed stage duration,
# and given total survivorship for each class.
# lambda is an initial estimate of lambda
# px = survivorship for class x
# d = duration (yrs) of stage x
# gamma is the fraction of individuals in the last year of stage x.

# Use a root-finding approach to solve the problem.
# Search for the value of lambda that minimizes the difference between the initial estimate (x)
# and the final estimate, as obtained by eigen-analysis (dominant eigenvalue):

# NOTE: in the current code, the final stage is not truncated!

# gs.A - basic function, x = lambda
#' @export
gs.A <- function(px, mx, d, x) {
	z <- px/x
	gamma <- (z^d - z^(d-1)) / ((z^d) - 1)
	P <- px*(1 - gamma)
	P <- c(P[-length(P)],px[length(px)]) # Do not truncate final stage, keep original estimate!!!
	G <- px*gamma
	F <- mx*P + c(mx[-1],0)*G
	A <- matrix(0,ncol=length(px),nrow=length(px))
	A[1,] <- F
	A[row(A)==col(A)] <- P # diagonal
	A[row(A)==col(A)+1] <- G[-length(G)] # sub-diagonal
	A
}

# Second function, combines the use of gs.A and gs.root (defined within).
# Returns the estimated lambda and the corresponding matrix A
#' @export
stage.duration <- function(px, mx, d, interval=c(0.01,1.5), tol=1e-15, maxiter=500) {

	inp <- data.frame(px=px, mx=mx, d=d)

	gs.root <- function(x) {
		A <- with(inp, gs.A(px,mx,d,x))
		L1 <- Re(eigen(A)$values[1]) # Dominant eigenvalue (lambda)
		x - L1
	}

	with(inp, {
		lambda <- uniroot(f=gs.root, interval=interval,tol=tol,maxiter=maxiter)$root
		A <- gs.A(px,mx,d,x=lambda)
		list(lambda=lambda,A=A)
	})
}

#' @export
life.span <- function(x = 0.9) {
  1/(1-x)
}

# life.span(x=0.8)
# x = adult survival
# curve(life.span(x), from=0, to=1, col=2, lwd=2,log="y")
