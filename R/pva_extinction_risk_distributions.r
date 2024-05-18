# Cumulative density function for extinction probability
# Box 3.3 and equation 3.5 in Morris & Doak
# Taken from popbio package
# mu = mean growth rate, sig2 = growth rate variance,
# Nc = current population size, Ne = Extinction threshold,
# tmax = projection interval
#' @export
erf <- function(y) {
  2 * pnorm(y * sqrt(2)) - 1
}

#' @export
phi <- function(z) {
  0.5 * (1 + (erf(z/sqrt(2))))
}

#' @export
xCDF <- function(t,mu,sigma2,d) {
  phi((-d - mu * t)/sqrt(sigma2 * t)) + exp(-2 * mu * d/sigma2) * phi((-d + mu * t)/sqrt(sigma2 * t))
}

#' @export
extCDF2 <- function (mu, sigma2, Nc, Ne, sq=1,
                     tmax=50, plot=TRUE,
                     col=4, lwd=2, lty=1,
                     add=FALSE, output=TRUE, log=NULL) {

	sig2 <- sigma2
	d <- log(Nc/Ne)
	if (d < 0) stop("Nc is already below Ne, nonsensical CDF")

	if (output == TRUE) {
		t <- seq(1,tmax,by=sq)
		G <- sapply(1:length(t), function(i) xCDF(t=t[i],mu=mu,sigma2=sig2,d=d))
		}

	if (plot == TRUE) {
			if (add == TRUE) {
				curve(xCDF(t=x,mu=mu,sigma2=sig2,d=d),
					from=0, to=tmax, n=1001, col=col, lwd=lwd, lty=lty, add=TRUE)
				}
			if (add == FALSE) {
				curve(xCDF(t=x,mu=mu,sigma2=sig2,d=d),
					from=0, to=tmax, n=1001, font.lab=2, col=col, lwd=lwd, lty=lty, bty="l",
					xlab="Years", ylab="Quasi-extinction probability before time T",
					main=paste("Quasi-extinction CDF: Nc =", Nc, ", Ne =", Ne),log=log)
				legend("topleft", as.expression(bquote(paste(mu==.(mu), ", ", sigma^2==.(sigma2)))),
					col=col, lwd=lwd, cex=1.1, bty="n")
			}
		}

	if (output == TRUE) {
		out <- data.frame(
			Time = t,
			P = G
		)
		invisible(out)
	}
}

# Extinction probability density function
# Returns a data frame with Time in column 1 and Probability (P) of quasi-extinction in column 2
# Additional parameters to mu and sig2 are:
# Nc = Current population size
# Ne = Quasi-extinction threshold
# sq = sequence interval (used for generating the time interval for which we want to plot the PDF)
# tmax = latest time to calculate extinction probability, default 50

#' @export
xPDF <- function(t, mu, sigma2, d) {
  (d / (sqrt(2*pi*sigma2*t^3))) * exp((-(d + mu*t)^2) / (2*sigma2*t))
}

#' @export
extPDF2 <- function (mu, sigma2, Nc, Ne, sq=1, tmax=50,
                     plot=TRUE, col=4, lwd=2, lty=1,
                     add=FALSE, output=FALSE) {

	sig2 <- sigma2
	d <- log(Nc/Ne)

	if (d <= 0) stop("Nc is already below or equal to Ne, nonsensical PDF")

	if (output == TRUE) {
		t <- seq(0,tmax,by=sq) # Differs from popbio (where t is 1:tmax)
		P <- (d / (sqrt(2*pi*sig2*t^3))) * exp((-(d + mu*t)^2) / (2*sig2*t))
		}

		if (plot == TRUE) {
			if (add == TRUE) {
				curve(xPDF(t=x,mu=mu,sigma2=sig2,d=d),
					from=0, to=tmax, n=1001, col=col, lwd=lwd, lty=lty, add=TRUE)
				}
			if (add == FALSE) {
				curve(xPDF(t=x,mu=mu,sigma2=sig2,d=d),
					from=0, to=tmax, n=1001, font.lab=2, col=col, lwd=lwd, lty=lty, bty="l",
					xlab="Years", ylab="Quasi-extinction probability",
					main=paste("Quasi-extinction PDF: Nc =", Nc, ", Ne =", Ne))
				legend("topright", as.expression(bquote(paste(mu==.(mu), ", ", sigma^2==.(sig2)))),
					col=col, lwd=lwd, lty=lty, cex=1.1, bty="n")
			}
		}

	if (output == TRUE) {
		out <- data.frame(
			Time = t,
			P = P)
		invisible(out)
	}
}

# Estimate uncertainty - bootstrapping
# Input variables:
# mu = mean log growth rate
# sig2 = variance of mean log growth rate
# nt = number of transitions in data set (= length of time series - 1)
# Nc = Current density
# Ne = Quasi-extinction threshold
# tq = length of the census in years
# tmax = latest time to calculate extinction probability, default 50
# Nboot = number of bootstrap samples

# Slightly modified function (made it compatible with my own extCDF2)
#' @export
countCDFxt2 <- function (mu, sigma2, nt, Nc, Ne, tq=nt, sq=1, tmax=50, Nboot=500,
    plot = TRUE) {

	# Make sure that inputs of mu and sig2 are numeric
	mu <- as.numeric(mu)
	sig2 <- as.numeric(sigma2)

    SEmu <- sqrt(sig2/tq) # Standard error of mu
    # Set up matrix for confidence interval estimation (bootstrapping)
	Glo <- matrix(1, tmax, 1) # Lower limit, fill with 1's
    Gup <- matrix(0, tmax, 1) # Upper limit, fill with 0's

	# Confidence intervals for mu and sig2
	CI_mu <- c(mu - qt(0.975, nt - 1) * SEmu, mu + qt(0.975,
        nt - 1) * SEmu)
    CI_sig2 <- c((nt - 1) * sig2/qchisq(0.975, nt - 1), (nt -
        1) * sig2/qchisq(0.025, nt - 1))

	###################
	# Changes made here
	# Note that sq must equal 1
	# "Best" estimate of G
	Gbest <- extCDF2(mu, sig2, Nc, Ne, sq, tmax)$P
    # Convert to matrix
	# Gbest <- Gbest[-1] # Remove year zero (t=0), previous version of extCDF2 evaluated from to = 0, now t = 1
	###################

	Gbest <- matrix(Gbest, tmax, 1)

	# Bootstrapping (repeat Nboot times)
	for (i in 1:Nboot) {

		# First set mu and sig2 to infinity, then while the random number (infinite, in this case) is outside the calculated confidence limits
		# create a random number (mu and sig2) within the confidence limits.
		murnd <- Inf
        while (murnd < CI_mu[1] | murnd > CI_mu[2]) {
            murnd <- mu + SEmu * rnorm(1)
        }
        sig2rnd <- Inf
        while (sig2rnd < CI_sig2[1] | sig2rnd > CI_sig2[2]) {
            sig2rnd <- sig2 * rchisq(1, nt - 1)/(nt - 1)
        }
        # Calculate the CDF for the randomly selected mu and sig2
		###################
		G <- extCDF2(murnd, sig2rnd, Nc, Ne, sq, tmax)$P # Updated code here as well
		# G <- G[-1] # Remove year zero
		###################

		# Output estimated confidence limits, write to Glo and Gup
		# Remember that Glo is a matrix filled with 1's and Gup a matrix filled with 0's
		for (x in 1:tmax) {
            if (G[x] < Glo[x]) {
                Glo[x] <- G[x]
            }
            if (G[x] > Gup[x]) {
                Gup[x] <- G[x]
            }
        }
    }
	# Plot
    if (plot) {
        plot(Gbest, log = "y", type = "l", pch = 16, col = "blue", font.lab=2, bty="l",
            ylim = c(min(Glo[Glo != 0]), max(Gup)), main = "Extinction CDF",
            xlab = "Years into the future", ylab = "Cumulative probability of quasi-extinction")
        # Add bootstrap confidence limits
		lines(Glo, col = "red", lty = 2)
        lines(Gup, col = "red", lty = 2)
    }
	# Print output data frame
    data.frame(Gbest, Glo, Gup)
}
