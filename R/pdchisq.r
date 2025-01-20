#Calculates pairwise chi-squared statistic for observations x given sd's s
#if s is a scalar function, it is applied to x to obtain an estimate of s
#
# Reference: Steele et al.
#
# This version (2018-09-10) is based on 
# msd.r
#

# No quantile calculations are included as Steele et al recommended simulation without critical values.
# Approximate crit vals are produced by the barplot method using qchisq

#
# NB: No separate file for the class plotting and print objects this time.
#

#
#Index-based calculation
#Still O(n^2) for distance calculation, but a lot faster and lighter on RAM than outer()

pdchisq <- function(x, s=sd, cov=NULL, cor = NULL, na.rm=FALSE, ...) {
        if(is.null(cov) ) {
		#No covariance matrix
		ss <- if(is.function(s)) {
			rep(s(x, na.rm=na.rm, ...), length(x))
		} else {
			if(length(s) == 1) {
				rep(s, length(x))
			} else {
				s
			}
		}
	        Cov <- if(!is.null(cor)) {
			outer(ss, ss, "*") * cor	        	
	        } else {
	        	diag(ss*ss)
	        }
        }
	
	pwch.local <- function(x, Cov, na.rm) {
		N <- length(x)
		I <- 1:N
        	#Correlated form
        	sapply(I, function(n) sum( (x[n] - x[-n])^2 /(Cov[n,n] + diag(Cov[-n,-n]) -2 * Cov[n, -n] ), na.rm=na.rm ) ) /
        		( N - 1 )
	
	}
	
        structure( 
        	pwch.local(x, Cov, na.rm),
        	names=if(is.null(names(x))) paste(1:length(x)) else names(x),
        	x=x,
        	s=sqrt(diag(Cov)),
        	cov=cov,
        	cor=cor,
        	class=c("PDchisq", "mtr.paircomp", "numeric")
        )
}

#Simple print and (bar)plot methods for PDchisq object

print.PDchisq <- function(x, ...) {
	print(c(x), ...)
	invisible(x)
}

plot.PDchisq <- function(x, type='h', ylab="Pair-difference chi-squared", ylim=NULL, ...) {
	plot.local<-function(x, y, axes, ...) plot(x, y, axes=FALSE, ...)
	ldots<-list(...)
	axes <- if(is.null(ldots$axes)) TRUE else ldots$ann 
	frame.plot <- if(is.null(ldots$frame.plot)) axes else ldots$frame.plot 
	
	if(is.null(ylim) ) {
		ylim <- range(pretty(c(0, x)))
	}
	
	plot.local(1:length(x), c(x), type=type, ylab=ylab, ylim=ylim, ...)
	if(frame.plot) box()
	if(axes) {
		axis(1, at=1:length(x), labels=names(x), ...)
		axis(2, ...)
	}
	return(invisible(NULL))
}

barplot.PDchisq <- function(height, ylab="Pair-difference chi-squared", names.arg=names(height), 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), 
	probs=c(0.95, 0.99), n=length(height), ylim=NULL, ... ) {
	
	if(is.null(names.arg)) names.arg <- paste(1:length(height))
	
	if(is.null(ylim) ) {
		ylim <- range(pretty(c(0, height)))
	}
	
	mids <- barplot(as.vector(height), ylab=ylab, names.arg=names.arg, ...)
	
	if(crit.vals && length(na.omit(probs)) > 0) {
		abline(h=qchisq(probs, n-1)/(n-1), lty=lty.crit, lwd=lwd.crit, col=col.crit)
	}
	return(invisible(mids))
}

bootPDchisq <- function(x, B=3000, probs=c(0.95, 0.99), 
	method=c("rnorm", "lhs"), keep=FALSE, labels=names(x), ...) {

	boot.mtr.pairwise(x=attr(x, "x"), s=attr(x, "s"), B=B, 
		cov=attr(x, "cov"), cor=attr(x, "cor"), probs=probs,
		stat="PDchisq", method=method, keep=keep, labels=names(x), ...)
}

