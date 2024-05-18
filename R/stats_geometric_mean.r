# geometric mean
#' @export
gmean <- function(x, method = "N") {
	if (method == "N") {
		n <- length(x)
		q <- length(diff(x)[!is.na(diff(x))])
		lambdas <- x[-1]/x[-length(x)]
	}
	if (method == "lambda") {
		lambdas <- x
		q <- length(diff(x)[!is.na(diff(x))])
	}
	#prod(lambdas,na.rm=T)^(1/q)
	exp(sum(log(lambdas))/q)
}
