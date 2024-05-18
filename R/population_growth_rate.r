#' @export
growth.rate <- function(x, lag.max=4, log=FALSE,
                        add.constant = c("never","if.zeros","always"),
                        constant=0, base=exp(1),
                        type=c("data.frame","list"),
                        na.omit=TRUE) {

	add.constant <- match.arg(add.constant)
	type <- match.arg(type)
	if (lag.max < 1) stop("number of lags must exceed 1") # What happens if I set lags to zero?

	if (class(x) == "ts") {
		start.x <- start(x)[1]
		x <- ts(x, start=start.x)
	} else {
		x <- ts(x)
	}

	if (add.constant == "never") {
		constant <- 0
	} else if (add.constant == "always") {
		constant <- constant
	} else if (add.constant == "if.zeros") {
		if (any(x == 0, na.rm=TRUE) == TRUE) {
			if (constant == 0) {
				stop("Set constant to a value different than 0!")
			} else if (constant != 0) {
				constant <- constant
			}
		}
	}

	x <- x + constant
	inds <- length(x)

	if (log==FALSE) {
		R <- ts(x[-1]/x[-inds], start=start.x+1)
	} else if (log==TRUE) 	{
		x <- log(x, base=base)
		R <- diff(x)
	} else if (log=="resp.only") {
		R <- diff(log(x, base=base))
	}

	R <- ts(R, start=(start(R)[1])-1)

	if (type == "list") {
		out <- lapply(0:lag.max, function(i) cbind(x=x,R=lag(R,i)))
		names(out) <- paste("Lag",0:lag.max,sep="")
	} else if (type == "data.frame") {
		E <- x
		for (i in 1:lag.max) {
			E <- cbind(E,lag(x,-i))
		}
		out <- cbind(E, R=R)
		colnames(out) <- c(paste("Lag",0:lag.max,sep=""),"R")
	}
	#if (na.omit == TRUE) return(na.omit(out))
	#if (na.omit == FALSE) return(out)
	out
}

# NOTE! This function does not account for irregular census intervals
# E.g. time series must have only a single deltat / frequency!
#' @export
gr <- function(N, method="r", plot=TRUE) {

	if (class(N) != "ts") stop("Object N must be a time series object")

	na.inds <- which(is.na(N))
	if (length(na.inds) > 0) message("Object N contains missing/NA values")

	g <- switch(method,
		lambda = N[-1]/N[-length(N)],
		r = diff(log(N))
	)

	# Time component of the variable g is wrong (it is associated with t instead of t-1)
	# and must therefore be changed
	g <- ts(g, start=start(N)[1], deltat=deltat(N))

	title.text <- switch(method,
		lambda = expression(bold(paste("Growth rate =", lambda))),
		r = expression(bold(paste("Growth rate = log(", lambda, ")")))
	)

	if (plot) {

		# LOESS for N
		loess.N <- loess(N ~ time(N))
		# LOESS for g
		loess.g <- loess(g ~ time(g))
		# LOESS for dd in g
		dd.N <- N[-length(N)]
		loess.dd <- loess(g ~ dd.N)
		pred.dd.N <- seq(min(N,na.rm=T),max(N,na.rm=T),length=101)
		pred.loess.dd <- predict(loess.dd, newdata=pred.dd.N)

		# Density dependence in growth rates
		lm.dd <- lm(g ~ N[-length(N)])
		coefs <- summary(lm.dd)$coefficients[,1]
		p <- summary(lm.dd)$coefficients[2,4]
		lm.text <- paste(round(coefs[1],6), " + ", round(coefs[2],6), "x, p =", round(p,6), sep="")

		par(mfrow=c(2,2))
		plot(N, las=1, type="n", main="Population counts", xlab="Time", ylab="Counts")
			lines(N, lty=2)
			points(N, pch=16, cex=1)
			points(time(N), predict(loess.N), type="l", col=2, lwd=2)
		plot(g,las=1,type="n",main=title.text,xlab="Time",ylab="Growth rate")
			lines(g,lty=2)
			points(g,pch=16,cex=1)
			points(time(g), predict(loess.g), type="l", col=2, lwd=2)
		plot(N[-length(N)],las=1,g,pch=16,cex=1,xlab="Population density",ylab="Growth rate",main="Density-dependence in growth rates")
			points(pred.dd.N,pred.loess.dd,col=2,lwd=2,type="l")
		acf(g, main="Auto-correlation in growth rates", na.action=na.pass)
		if (length(na.inds) > 0) title(sub="N contains NA values")
		par(mfrow=c(1,1))
	}
	# Print output
	invisible(list(TimeSeries=N, GrowthRate=g))
}

# Estimate growth rate with linear regression (Dennis approach)
# Required input is a time series object with population counts (N).
# Morris & Doak (2002), p. 66-79
# This function should be expanded to:
# a) allow for multiple regression (i.e. split the time variable in, and create a model.matrix)
# purpose: get multiple values of mu, but only a single estimate of sigma2
# b) investigate temporal trends in mu and sigma2 (p. 94 in Morris & Doak)
#' @export
gr.da <- function(N, plot=TRUE) {

	if (class(N) != "ts") stop("Object N must be a time series object")

	t <- as.numeric(time(N))
	na.inds <- which(is.na(N))

	if (length(na.inds) > 0) {
		N <- N[-na.inds]
		t <- t[-na.inds]
	}

	gr <- log(N[-1] / N[-length(N)]) # Growth rate

	x <- sqrt(diff(t))
	y <- gr/x

	# model.matrix(~x-1)

	fm.da <- lm(y ~ x-1)

	# MEAN (slope)
	mu <- as.numeric(coef(fm.da))
	# VAR (mean square in analysis of variance table)
	sigma2 <- anova(fm.da)[["Mean Sq"]][2]

	df1 <- length(gr)-1
	CI.mu <- confint(fm.da, level=0.95)
	CI.sigma2 <- matrix(df1*sigma2 /qchisq(c(.975, .025), df = df1), nrow=1, ncol=2)
	colnames(CI.sigma2) <- colnames(CI.mu)
	rownames(CI.sigma2) <- rownames(CI.mu)

	if (plot==TRUE) {
		# plot like Fig 3.7 in Morris & Doak
		plot(x,y, xlim=c(0,max(x)), ylim=c(min(y),max(y)), pch=16, las=1, bty="l",
		xlab=expression((t[t+1]-t[i])^{1/2}),
		ylab=expression(log(N[t+1]/N[t]) / (t[t+1]-t[i])^{1/2}) ,
		main=expression(paste("Estimating ", mu, " and ", sigma^2)))
		abline(fm.da, lwd=2, col=2, lty=2)
		abline(h=0, lwd=1, lty=3)
		legend("topleft", as.expression(bquote(paste(mu==.(mu), ", ", sigma^2==.(sigma2)))),
			col=2, lwd=2, lty=2, cex=0.8, bty="n")
	}

	out <- list(
		mu = c(mu = mu , '2.5 %' = CI.mu[1], '97.5 %' = CI.mu[2]),
		sigma2 = c(sigma2 = sigma2, '2.5 %' = CI.sigma2[1], '97.5 %' = CI.sigma2[2]),
		n = length(y),
		model = fm.da
	)
	out
}
