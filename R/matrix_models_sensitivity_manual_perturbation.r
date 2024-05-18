# Evaluate the effect of A[i,j] on lambda, by manually changing A[i,j].
# from,to,by arguments defines over which interval the relation between matrix element A[i,j] and lambda is evaluated.
# For instance a survival rate could be evaluated with S.m.eval(A=A,i=6,j=6,from=0,to=1,by=0.01)
# Ouput is a plot with the numerical result, as well as the analytical solution for S, evaluated at A[i,j].
#' @export
S.m.eval <- function(A,i,j,from=0,to=1,by=0.01) {

	x <- seq(from,to,by) # Vector with values, evaluate lambda for all values in x
	# For each value of x, calculate the resulting lambda
	lambda.m <- sapply(1:length(x), function(k) {
		A.temp <- A
		A.temp[i,j] <- x[k]
		eigenall(A.temp)$lambda
		}
	)

	# Plot analytical solution, S at current value of A[i,j]
	S <- S.tangent(A,i,j)
	a <- as.numeric(S["a"])
	S.ij <- as.numeric(S["S.ij"])

	# Create a graph
	me <- paste(i,j,sep="")
	plot(x,lambda.m,type="l",lwd=2,
		xlab=substitute(paste("Value of matrix element ", A[me])),
		ylab=expression(lambda),las=1,
		main="Sensitivity")
	abline(v=A[i,j],lty=2,col=2)
	abline(a=a, b=S.ij, lty=1, lwd=1, col=1)

	legend("topleft", c("Current value","Analytical solution","Manual perturbation"),
		col=c(2,1,1), lty=c(2,1,1), lwd=c(1,1,2), bty="n", cex=1)

	# Collect output
	out <- c(a=a, b=S.ij)
	invisible(out) # Print output, in this case I don't want the output to be displayed(visible, unless I store the output as an object

}

# Predict lambda at values of "predict", i.e. a vector of values for A[i,j]
# Intended for use only together with S.tangent
#' @export
S.m.predict <- function(a,b,predict) a + b*predict

# S.tangent returns the intercept and slope (sensitivity) of the tangent to a curve, evaluated at A[i,j]
#' @export
S.tangent <- function(A,i,j) {
		eig <- eigenall(A) # Eigen-analysis of matrix A
		A.ij <- A[i,j]
		lambda.ij <- eig$lambda
		S.ij <- eig$S[i,j] # Slope of line
		a <- lambda.ij - (S.ij*A.ij) # Intercept of line
		c(a=a,S.ij=S.ij) # Parameters for a line given by the equation a + bx, where b=S.ij
}

# Calculate sensitivity with the manual perturbation method.
# Requires a transition matrix A,
# which element in the matrix that is changed A[i,j]
# where i indicates row in matrix A, and j indicates column in matrix A.
# delta is the magnitude of perturbation (i.e. A[i,j] + delta)
#' @export
S.m <- function(A,i,j,delta) {
	r.orig <- A[i,j] # Extract current value
	eig.orig <- eigenall(A) # Eigen-analysis
	A.new <- A # Create a copy of the original matrix
	r.new <- A[i,j] + delta # Calculate new rate for matrix element A[i,j] by adding a SMALL value, delta
	A.new[i,j] <- r.new # Add the new matrix element rate to the matrix A.new
	eig.new <- eigenall(A.new) # Eigen-analysis of the new matrix
	S <- (eig.new$lambda - eig.orig$lambda) / (r.new - r.orig) # Calculate S (sensitivity)
	E <- (r.orig / eig.orig$lambda) * S # Convert to E (elasticity)
	c(S=S,E=E) # Bind values for S and E together and print to screen
}
