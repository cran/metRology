#
# New boot function for pairwise stats
#

boot.mtr.pairwise<-function(x, s=mad , B=3000, probs=c(0.95, 0.99), cov=NULL, cor = NULL, 
	 stat=c("MSD", "PDchisq"), method=c("rnorm", "lhs"), keep=FALSE, 
	 labels=names(x), ...) {

 	method <- match.arg(method)
	
	stat <- match.arg(stat)
	
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
        #c() strips unwanted attributes from msd object
        t0 <- switch(stat, 
			MSD=c(msd(x, ss)),
			PDchisq=c(pdchisq(x, ss, cov=cov, cor = cor))
		     )

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
	
	#Generate the simulation
	t <- switch(stat, 
		MSD=t(apply(m, 1, function(x, s) as.vector(msd(x, s)), s=ss)),
		PDchisq=t(apply(m, 1, function(x, s, ...) as.vector(pdchisq(x, s, ...)), 
			s=ss, cov=cov, cor=cor))
	     )

	#Quantiles (critical values)
	ncrit <- length(na.omit(probs))
	if(ncrit > 0) q <- apply(t, 2, quantile, probs=probs)
	if(ncrit==1) q <- matrix(q, nrow=1)
	
	p <- apply(sweep(t, MARGIN=2, STATS=t0, FUN="-"), 2, function(x) sum(x > 0) )
	p <- p/B
	
	structure(list(t0=t0, labels=labels, probs=probs, critical.values=q, 
			pvals=p, B=B, method=method, stat=stat, t=if(keep) t else NA),
			class="bootMtrPairs")
        
}


print.bootMtrPairs <- function(x, ...) {
	print(c(x$t0), ...)
}

summary.bootMtrPairs <- function(object, p.adjust="none", ...) { 
	#'digits=NULL' removed 2017-09-02 - SLRE
	p.adjust <- match.arg(p.adjust, p.adjust.methods)
	p.adj <- p.adjust(object$pvals, method=p.adjust)
	structure(c(object[c("t0", "labels", "probs", "critical.values", "pvals")], 
		list(p.adjust=p.adjust, p.adj=p.adj),
		object[c('B', "method", "stat")]),
		class="summary.bootMtrPairs")
}

print.summary.bootMtrPairs <- function(x, digits=3, ..., signif.stars = getOption("show.signif.stars"), 
		signif.legend=signif.stars) 
{
	cat(switch(x$stat,
		MSD="Median Scaled Difference parametric bootstrap\n",
		PDchisq="Pair-difference chi-square parametric bootstrap\n"))
	cat(sprintf("%s replicates\n", format(x$B)))
	cat(sprintf("Sampling method: %s\n", x$method))
	cat(sprintf("P-value adjustment: %s\n", x$p.adjust))
	df.crit <- as.data.frame(t(x$critical.values), check.names=FALSE)
	names(df.crit) <- paste("Upper", names(df.crit))
	m <- format( cbind(data.frame(t0=x$t0, row.names=x$labels), 
		     df.crit, P=x$p.adj) )
	which.0 <- which(x$pvals < 1/x$B)
	if(length(which.0) > 0) {
		#Recalculate adjusted p '+1' for zeros
		p.plus <- p.adjust( pmax(1 + x$B * x$pvals)/x$B,
					method=x$p.adjust )
		m$P[which.0] <- sprintf(" < %7.1e", p.plus[which.0] )
	}
	
	names(m)[c(1, length(m))] <- switch(x$stat, 
		MSD=c("MSD", "P(>MSD)"), 
		PDchisq=c("PDchisq", "P(>PDchisq)"))
	
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

plot.bootMtrPairs <- function(x, ...) {
	UseMethod("barplot")
}

barplot.bootMtrPairs <- function(height, ylab=NULL, names.arg=height$labels, 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), ylim=NULL, ... ) {
	
	if(is.null(names.arg)) names.arg <- paste(1:length(height$t0))
	if(is.null(ylab)) ylab <- switch(height$stat, MSD="MSD", PDchisq="Pair-difference chi-square")
	
	if(is.null(ylim) ) {
		if( crit.vals ) {
			# If plotting crit vals, 
			# make sure range includes critical values too
			ylim <- range(pretty(c(0, height$t0, height$critical.values)))
			
		} else {
			ylim <- range(pretty(c(0, height$t0)))
		}
	}
	
	mids <- barplot(height$t0, ylab=ylab, names.arg=names.arg, ylim=ylim, ...)
	
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

