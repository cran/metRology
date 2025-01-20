#
# metrology\R\MSD-class.r
#

# Author: S L R Ellison
# Created: 2017-09-04


#Simple print and (bar)plot methods for MSD object

print.MSD <- function(x, ...) {
	print(c(x), ...)
	invisible(x)
}

plot.MSD <- function(x, type='h', ylab="MSD", ylim=NULL, ...) {
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

barplot.MSD <- function(height, ylab="MSD", names.arg=names(height), 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), 
	probs=c(0.95, 0.99), n=length(height), ylim=NULL, ... ) {
	
	if(is.null(names.arg)) names.arg <- paste(1:length(height))
	
	if(is.null(ylim) ) {
		ylim <- range(pretty(c(0, height)))
	}

	mids <- barplot(as.vector(height), ylab=ylab, names.arg=names.arg, ylim=ylim, ...)
	
	if(crit.vals && length(na.omit(probs)) > 0) {
		abline(h=qmsd(probs, n), lty=lty.crit, lwd=lwd.crit, col=col.crit)
	}
	return(invisible(mids))
}

