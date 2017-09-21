#
# R\bootMSD.R
#

#
# Author: S L R Ellison
# Created: 2016-08-30
# This version: 2017-09-02
# This version omits latin hypercube sampling from simulation methods
# to reduce unnecessary package loads

# lhs code is retained, commented against future development

# Calculates parametric bootstrapped median scaled difference for observations x given sd's s
# if s is a scalar function, it is applied to x to obtain an estimate of s

#Define generic
bootMSD <- function(x, ...) {
	UseMethod("bootMSD")
}

#Define variant for class mtR.msd (metRology-msd)

bootMSD.MSD <- function(x, B=3000, probs=c(0.95, 0.99), 
	method=c("rnorm", "lhs"), keep=FALSE, labels=names(x), ...) {

	bootMSD.default(attr(x, "x"), attr(x, "s"), B=B, probs=probs, 
		method=method, keep=keep, labels=names(x), ...)
}

bootMSD.default<-function(x, s=mad , B=3000, probs=c(0.95, 0.99), 
	method=c("rnorm", "lhs"), keep=FALSE, labels=names(x), ...) {

 	method <- match.arg(method)

        #Get standard deviations if not a vector
        ss <- if(is.function(s)) {
                rep(s(x, ...), length(x))
        } else {
                if(length(s) == 1) {
                        rep(s, length(x))
                } else {
                        s
                }
        }
        
        #Get the raw values of msd from the data
        t0 <- c(msd(x, ss)) #c() strips unwanted attributes from msd object
        
	N <- length(x)
	
	#Generate simulated data
	if(method=="rnorm") {
		m <- matrix(rnorm(B*N, mean=0, sd=ss), byrow=TRUE, ncol=N)
	# lhs code omitted in this version
	#} else if(method=="lhs") {
	#	require(lhs) ## Shouold be Imported for pkg if used
	#	m <- sweep(qnorm(randomLHS(B, N)), MARGIN=2, STATS=ss, FUN='*')
	} else {
		stop(paste("Method", method, "not implemented"))
	}
	
	#Generate the msd simulation
	t <- t(apply(m, 1, function(x, s) as.vector(msd(x, s)), s=ss))

	#Quantiles (critical values)
	ncrit <- length(na.omit(probs))
	if(ncrit > 0) q <- apply(t, 2, quantile, probs=probs)
	if(ncrit==1) q <- matrix(q, nrow=1)
	
	p <- apply(sweep(t, MARGIN=2, STATS=t0, FUN="-"), 2, function(x) sum(x > 0) )
	p <- p/B
	
	structure(list(msd=t0, labels=labels, probs=probs, critical.values=q, 
			pvals=p, B=B, method=method, t=if(keep) t else NA),
			class="bootMSD")
        
}

print.bootMSD <- function(x, ...) {
	print(c(x$msd), ...)
}

summary.bootMSD <- function(object, p.adjust="none", ...) { 
	#'digits=NULL' removed 2017-09-02 - SLRE
	p.adjust <- match.arg(p.adjust, p.adjust.methods)
	p.adj <- p.adjust(object$pvals, method=p.adjust)
	structure(c(object[c("msd", "labels", "probs", "critical.values", "pvals")], 
		list(p.adjust=p.adjust, p.adj=p.adj),
		object[c('B', "method")]),
		class="summary.bootMSD")
}

print.summary.bootMSD <- function(x, digits=3, ..., signif.stars = getOption("show.signif.stars"), 
		signif.legend=signif.stars) 
{
	cat("Median Scaled Difference parametric bootstrap\n")
	cat(sprintf("%s replicates\n", format(x$B)))
	cat(sprintf("Sampling method: %s\n", x$method))
	cat(sprintf("P-value adjustment: %s\n", x$p.adjust))
	df.crit <- as.data.frame(t(x$critical.values), check.names=FALSE)
	names(df.crit) <- paste("Upper", names(df.crit))
	m <- format( cbind(data.frame(MSD=x$msd, row.names=x$labels), 
		     df.crit, "P(>MSD)"=x$p.adj) )
	which.0 <- which(x$pvals < 1/x$B)
	if(length(which.0) > 0) {
		#Recalculate adjusted p '+1' for zeros
		p.plus <- p.adjust( pmax(1 + x$B * x$pvals)/x$B,
					method=x$p.adjust )
		m[["P(>MSD)"]][which.0] <- sprintf(" < %7.1e", p.plus[which.0] )
	}
	if (signif.stars && any( x$p.adj < 0.1 )) {
		Signif <- symnum(x$p.adj, corr = FALSE, na = FALSE, 
		cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
			symbols = c("***", "**", "*", ".", " "))
                m <- cbind(m, " "=format(Signif))
	} else {
	       #Nothing significant so no legend
	       signif.legend <- FALSE
	}
 	print(m, digits=digits, ...)
	if (signif.legend) {
		#Using code borrowed from print.Coefmat ...
		if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, 
		    "legend"))) 
		    sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
		cat("---\nSignif. codes:  ", sleg, sep = "", fill = w + 
		    4 + max(nchar(sleg, "bytes") - nchar(sleg)))
	}
 	invisible(x)
}

plot.bootMSD <- function(x, ...) {
	UseMethod("barplot")
}

barplot.bootMSD <- function(height, ylab="MSD", names.arg=height$labels, 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), ... ) {
	
	mids <- barplot(height$msd, ylab=ylab, names.arg=names.arg, ...)
	
	if(crit.vals && length(na.omit(height$probs)) > 0) {
		ncrit <- length(height$probs)
		dmid <- diff(mids)[1]
		cw <- 0.98*dmid/2
		lty.crit <- rep(lty.crit, ncrit)
		lwd.crit <- rep(lwd.crit, ncrit)
		col.crit <- rep(col.crit, ncrit)
		for(i in 1:ncrit) segments(mids-cw, height$critical.values[i,], mids+cw, height$critical.values[i,],
			lty=lty.crit[i], lwd=lwd.crit[i], col=col.crit[i], lend=2)
	}
	return(invisible(mids))
}

