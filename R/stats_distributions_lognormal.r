# Convert from/to lognormal and normal distributions
#' @export
lognormal <- function(M,V,method=c("lnorm","exp")) {

	method <- match.arg(method)
	M <- as.numeric(M)
	V <- as.numeric(V)
	if (method == "exp") {
		# Use method="exp" if you have data on log-scale (mu and sigma2 in Morris & Doak)
		Mn <- M; Vn <- V
		Mln <- exp(Mn + 0.5*Vn)
		Vln <- exp(2*(Mn + 0.5*Vn)) * (exp(Vn) - 1)
		out <- c(M=Mln, V=Vln)
	} else if (method == "lnorm") {
		# Use method="lnorm" if you have data on original scale (e.g. lambdas)
		Mln <- M; Vln <- V
		Mn <- log(Mln) - 0.5*log((Vln/Mln^2) + 1)
		Vn <- log((Vln/Mln^2) + 1)
		out <- c(Mn=Mn, Vn=Vn)
	}
	out
}

# Function that calculates mean and variance of a simulated lognormal variable
#' @export
lognormal.sim <- function(n, mu, sigma) {
	x <- rlnorm(n, mu, sigma)
	c(mean=mean(x),var=var(x))
}
