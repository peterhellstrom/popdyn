# fderiv
# FUNCTION THAT CALCULATES AND PLOTS
# - first derivative
# - second derivative
# Required input:
# f = an expression that contains the function to be evaluated
# parms = a vector or list with defined parameters for fn
#' @export
fderiv <- function(f,parms=NULL,from=0,to=40,n=1001,tol=1e-15,tangent=TRUE) {

	if(!is.null(parms)) {
		for(x in unlist(attributes(parms)))
		assign(x, parms[[x]],envir = .GlobalEnv)
	}

	f0 <- deriv(body(f), "x", func=TRUE)
	# First derivative
	f.d1 <- D(body(f),"x")
	f1 <- deriv(f.d1,"x",func=TRUE)
	# Second derivative
	f.d2 <- D(f.d1, "x")
	f2 <- deriv(f.d2,"x",func=TRUE)

	# Turnpoints of first derivative
	tp <- uniroot.all(f=f1, interval=c(from,to), n=n, tol=tol)
	# Inflection points
	ip <- uniroot.all(f=f2, interval=c(from,to), n=n, tol=tol)

	tangent.d1 <- function(x) {
		slope <- f1(x)
		y.at.x <- f0(x)
		a <- y.at.x - slope*x
		b <- slope
		abline(a,b,lty=2)
	}

	par(mfrow=c(1,3))
		# Function
		do.call("curve", list(body(f), from=from, to=to, n=1001, font.lab=2, cex.lab=1.3, cex.axis=1.3, xlab='x', ylab="f (x)", main="Function"))
		points(ip, f0(ip), pch=16, col=2, cex=1.4)
		if (tangent == TRUE) if (length(ip) > 0) for (i in 1:length(ip)) tangent.d1(x=ip[i])

		# First derivative
		do.call("curve", list(f.d1, from=from, to=to, n=1001, font.lab=2, cex.lab=1.3, cex.axis=1.3, xlab='x', ylab="f '(x)", main="First derivative"))
		points(ip, f1(ip), pch=16, col=2, cex=1.4)

		# Second derivative
		do.call("curve", list(f.d2, from=from, to=to, n=1001, font.lab=2, cex.lab=1.3, cex.axis=1.3, xlab='x', ylab="f ''(x)", main="Second derivative"))
		abline(h=0, col=1, lty=2)
		points(ip, y = rep(0, length(ip)), pch = 16, col=2, cex=1.4)
	par(mfrow=c(1,1))

	# Print output
	n.ip <- length(ip)

	if (n.ip >= 1) {
		out <- list(
			n = n.ip,
			InflectionPoints = ip,
			MaxFirstDeriv = f1(ip),
			TpFirstDeriv = tp)
		invisible(out)
	}

}

# fderiv.simple
#' @export
fderiv.simple <- function(f,parms,from,to) {

  if (class(f) != "function") stop("The argument 'f' must be a function")

  f.d1 <- D(body(f),"x")
  f.d2 <- D(f.d1, "x")

  with(parms, {
    par(mfrow=c(1,3))
      do.call("curve", list(body(f), from=from, to=to, n=1001,
        xlab='x', font.lab=2, cex.lab=1.3, cex.axis=1.3, ylab="f(x)", main="Function"))
      do.call("curve", list(f.d1, from=from, to=to, n=1001,
        xlab='x', font.lab=2, cex.lab=1.3, cex.axis=1.3, ylab="f '(x)", main="First derivative"))
      do.call("curve", list(f.d2, from=from, to=to, n=1001,
        xlab='x', font.lab=2, cex.lab=1.3, cex.axis=1.3, ylab="f ''(x)", main="Second derivative"))
      abline(h=0, lty=2)
      par(mfrow=c(1,1))
  })
}
