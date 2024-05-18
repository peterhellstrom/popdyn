# Create a function for use together with uniroot():
# Inputs:
# x = the parameter to evaluate (necessary input argument for uniroot, should not be supplied by the user)
# expr = an expression containing the elements of the matrix. This expression is used to
# populate the transition matrix, and it is important to remember that the resulting matrix
# is filled by row, starting with the top (fecundities) row.
# pars = a named list with demographic parameters, excluding the parameter given by argument vitalrate
# vitalrate = the name of the demographic parameter that should be evaluated (written as character)
# lambda = the lambda value to maintain, defaults to 1. This value can be changed, e.g. if you want to
# calculate a fitness landscape
# Output: the value of the vital rate specified in "vitalrate" that maintains lambda, given the
# other vital rates in the transition matrix.

#' @export
maintain.lambda <- function(x, expr, pars, vitalrate=NULL, lambda=1) {
	if (length(vitalrate) > 1) stop("Can only replace one vitalrate")
	expr.tmp <- gsub(pattern=vitalrate, replacement="x", expr) # replace pattern/vitalrate with "x"
	expr.tmp <- as.expression(parse(text=expr.tmp)) # convert back to expression
	pars$x <- x # supply the parameter x in the parameter list
	v <- sapply(expr.tmp, eval, pars, NULL) # calculate the transition matrix elements
	A <- matrix(v, nrow=sqrt(length(v)), byrow=TRUE) # populate the transition matrix, this is done by row!
	lambda - as.numeric(abs(eigen(A)$values[1])) # the root. Calculate the dominant eigenvalue, and subtract this from the desired lambda
}
