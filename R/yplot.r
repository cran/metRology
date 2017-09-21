#
# yplot_v2.r
# 
# Functions to provide and support Youden plots
#
# Created: 2015-08-22
# Author: S Ellison
#
# Amended (location/center) 2016-05-27
#


youden.plot <- function(x, ...) {
	UseMethod("youden.plot")
}

yplot <- function(x, ...) {
	UseMethod("youden.plot")
}

#youden.plot.data.frame <- function(x, ..., main=paste("Youden plot: ", deparse(substitute(x)))) {
#	youden.plot(as.matrix(x), ..., main=main)	
#}

youden.plot.default <- function(x, y=NULL, type=c("points", "labels", "both", "outliers"), 
	labels, probs=c(0.95, 0.99), 
	x0, y0, pch=par("pch"), cex=par("cex"), col=par("col"), bg=par("bg"),
	main, xlab, ylab, xlim=c("data", "ellipse", "all"), ylim=c("data", "ellipse", "all"), 
	col.axes=2, lwd.axes=1, lty.axes=1, cex.lab=0.7, pos=3, out.method=c("F", "chisq", "n"), n.out, p.out=0.99, 
	add=FALSE, ...) {
	
	#setting col.axes=NA suppresses plotting of youden plot axes.
	
	type <- match.arg(type)
	
	#Check dimensions
	if( is.null(y) ) {
		X <- x
	} else {
		X <- cbind(x, y)
		colnames(X) <- c(deparse(substitute(x)), deparse(substitute(y)))
		if(!is.null(names(x))) rownames(X) <- names(x)
	}
	
	if(is.null(dim(X)) || dim(X)[2] != 2)
		stop("x must be a 2-colum matrix and y omitted, or x and y must both be present and vectors")

	cov.dellipse.local <- function(x, y=NULL, cov.method=NULL, 
		scalefn=NULL, locfn=NULL, cov.control=list(), ...) {
	
		cov.dellipse(x=x, y=y, cov.method=cov.method, 
			scalefn=scalefn, locfn=locfn, cov.control=cov.control)
	}
	
	cov.data <- cov.dellipse.local(X, ...)

	if(missing(x0)) x0 <- cov.data$center[1]
	if(missing(y0)) y0 <- cov.data$center[2]
		#lets x0, y0 default to calculated values while allowing user control
	n <- cov.data$n
		
	ellipses <- data.ellipse(cov.data, probs=probs, plot=FALSE)
	e.xrange <- range(sapply( ellipses$ellipses, function(e) range(e[,1]) ))
	e.yrange <- range(sapply( ellipses$ellipses, function(e) range(e[,2]) ))
	
	if(is.character(xlim)) {
		xltype <- match.arg(xlim)
		if(xltype=="data") xlim <- range(pretty(X[,1]))
		else if(xltype=="ellipse") xlim <- range(pretty(e.xrange))
		else xlim <- range(pretty(c(X[,1], e.xrange)))
	}
	
	if(is.character(ylim)) {
		yltype <- match.arg(ylim)
		if(yltype=="data") ylim <- range(pretty(X[,2]))
		else if(yltype=="ellipse") ylim <- range(pretty(e.yrange))
		else ylim <- range(pretty(c(X[,2], e.yrange)))
	}
	
	if(missing(main)) main <- paste("Youden plot:", deparse(substitute(x)))
	
	if(missing(xlab)) xlab <- if(!is.null(colnames(X))) 
					colnames(X)[1]
				  else deparse(substitute(x))
	if(missing(ylab)) ylab <- if(!is.null(colnames(X))) 
					colnames(X)[2]
			          else deparse(substitute(y))
	
	# dot-masking local plot function 
	lplot <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, 
	    log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
	    ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL, 
	    panel.last = NULL, asp = NA, ...) {
	    
	    plot(x = x, y = y, type = type, xlim = xlim, ylim = ylim, 
	    	    log = log, main = main, sub = sub, xlab = xlab, ylab = ylab, 
	    	    ann = ann, axes = axes, frame.plot = frame.plot, panel.first = panel.first, 
	            panel.last = panel.last, asp = asp) #note absence of ...
	}

	#set up plot window unless add=TRUE
	if(!add) lplot(mean(xlim), mean(ylim), type="n", xlim=xlim, ylim=ylim, 
		xlab=xlab, ylab=ylab, main=main, ...) 
	plot(ellipses, xlim=xlim, ylim=ylim, add=TRUE,  ...)

	abline( h=y0 , v=x0, col=col.axes, lty=lty.axes, lwd=lwd.axes )
	
	if(missing(labels)) {
		labels <- if(!is.null(row.names(X))) row.names(X) else paste(1:nrow(X))
	}

	if(type %in% c("points", "outliers", "both" ) )
		points(X[,1], X[,2],  pch=pch, col=col, cex=cex, bg=bg)
	if(type=="labels") text(X[,1], X[,2],  labels, cex=cex.lab)
	if(type=="both") text(X[,1], X[,2],  labels, cex=cex.lab, pos=pos)
	if(type=="outliers") {
		out.method <- match.arg(out.method)
		md <- mahalanobis(X, center=c(x0, y0), cov=cov.data$cov)
		if(missing(n.out)) n.out <- min(n, max(5, floor(n/10)))
		out.md <- if(out.method=="n") 
			which(rank(md) > n - n.out)
		     else if(out.method=="chisq")
			which(md > qchisq(p.out, 2) ) #Chi^2, 2df
		     else
			which( md > 2 * (n-1) * qf(p.out, 2, n-1) / (n-2) ) #F dist
		     
		ellipses$outliers <- list(method=out.method, n.out=n.out, p.out=p.out, 
			which=out.md, labels=labels[out.md], coords=X[out.md,] )
		if(length(out.md)) 
			text(X[out.md,1], X[out.md,2],  labels[out.md], cex=cex.lab, pos=pos)
	}
	return(invisible(ellipses))
}

data.ellipse <- function(cov, probs=0.95, plot=TRUE, npoints=100, ...) {
        ##cov is an object of class cov.ellipse
        #probs is a vector of confidence levels
                
	
	Cov <- cov$cov
	n <- cov$n.obs
	x0 <- cov$center[1]
	y0 <- cov$center[2]
	sx <-  cov$scale[1]
	sy <-  cov$scale[2]

        rho <- Cov[1,2]/sqrt(Cov[1,1]*Cov[2,2])
        
        T2 <- if(is.na(n)||is.infinite(n)) 
                2 * qf(probs, 2, Inf)
              else
                2 * (n-1) * qf(probs, 2, n-1) / (n-2)
        
        L <- length(probs)
        ellipses <- as.list(rep(NA, L))
        names(ellipses) <- paste("p=", probs, sep="")
        for(i in 1:L) {

              zx.1 <- sqrt(T2[i]) * cos( seq(pi, 0, length.out=npoints))  
              zx.2 <- rev(zx.1[-c(1,length(zx.1))])
              zy <- rho * zx.1 + sqrt( pmax( (1 - rho*rho) * ( T2[i] - zx.1*zx.1), 0 ) )
              zy <- c(zy, rho * zx.2 - sqrt( pmax( (1 - rho*rho) * ( T2[i] - zx.2*zx.2), 0 ) ))
              
              x <- c(zx.1, zx.2) * sx + x0
              y <- zy * sy + y0
              
              ellipses[[i]] <- cbind(x,y)
        }
	
	scale <- c(sx,sy)
	centre <- c(x0, y0)
	names(centre) <- names(scale) <- c("x","y")

	rv <- structure(list(ellipses=ellipses, probs=probs, cov=cov), #cov contains center, scale etc
		class="d.ellipse" )
	if(plot) {
	      plot(rv, ...)
	}
	
	return(invisible(rv))
}

plot.d.ellipse <- function(x, col.ellipse=1, lty.ellipse=1, lwd.ellipse=1, fill=NA, density=NULL, angle=45, 
                add=FALSE, npoints=100, xlim=NA, ylim=NA, 
                prinax=FALSE, col.prinax=1, lty.prinax=1, lwd.prinax=1,  
                xlab=NULL, ylab=NULL, ...) {
	
        L <- length(x$ellipses)
        col <- rep(col.ellipse, length.out=L)
        lty <- rep(lty.ellipse, length.out=L)
        lwd <- rep(lwd.ellipse, length.out=L)
        fill <- rep(fill, length.out=L)
        if(!is.null(density[1])) density <- rep(density, length.out=L)
        angle <- rep(angle, length.out=L)
	
	if(!add) {
		if(is.na(xlim[1])) xlim<-c(min(sapply(x$ellipses, function(x) range(x[,1]))), max(sapply(x$ellipses, function(x) range(x[,1])))) 
		if(is.na(ylim[1])) ylim<-c(min(sapply(x$ellipses, function(x) range(x[,2]))), max(sapply(x$ellipses, function(x) range(x[,2])))) 
		if(is.null(xlab)) {
			xlab <-dimnames(x$cov$cov)[[1]][1]
			if(is.null(xlab)) xlab <- "X"
		}
		if(is.null(ylab)) {
			ylab <-dimnames(x$cov$cov)[[1]][2]
			if(is.null(ylab)) ylab <- "Y"
		}
		
		plot(mean(xlim), mean(ylim), type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...) 

	} 

	for(i in 1:L ) {
		polygon(x$ellipses[[i]], border=col[i], col=fill[i], lty=lty[i], lwd=lwd[i],
			density=density[i], angle=angle[i])
	}

	centre<-x$cov$center
	Cov <- x$cov$cov
	if(prinax) {
		do.prinax <- function(vec, centre, ...) {
			b <- vec[2]/vec[1]
			a <- centre[2] - b * centre[1]
			abline(a,b, col=col.prinax, lty=lty.prinax, lwd=lwd.prinax)
		}
		ev <- eigen(Cov, symmetric=TRUE)$vectors
		apply(ev, 2, do.prinax, centre=centre, col=col.prinax, lwd=col.prinax, lty=lty.prinax) 	
	}
	return(invisible(NULL))
}

summary.d.ellipse <- function(object, ...) {
	lapply( object, summary, ...	)
}


print.d.ellipse <- function(x,  ...) {
	hprint.list <- function(xx, name, ...) {
		cat(sprintf("\n$%s", name))
		lapply(names(xx), function(nn, name) {
			cat(sprintf("\n$%s$%s\n", name, nn))
			print(xx[[nn]])
		}, name=name)
	}
	cat("$ellipses\n")
	print(summary(x$ellipses), ...)
	cat("\n$probs\n")
	cat(paste(c(format(x$probs, ...), "\n"), collapse="  "))
	hprint.list(x$cov, 'cov', ...)
	if(!is.null(x$outliers)) hprint.list(x$outliers, 'outliers')
	return(invisible(x))
}


cov.dellipse <- function(x, y=NULL, cov.method=c("spearman", "kendall","pearson","MCD","OGK","GK","gk","rgk","mcd", "mve"), 
	scalefn=NULL, locfn=NULL, cov.control=list()) {
	#Returns an object of class cov.dellipse, which is a list with (at least) components
	#    method	Character string describing method
	#    cov	2x2 covariance matrix
	#    cor	2x2 correlation matrix
	#    center	vector (length 2) specifying centre of ellipse
	#    scale	vector, length 2, specifying scale estimates for each variable
	#    n.obs	number of points (rows) used in the covariance estimate
	#
	# This list is intended to be consistent with that returned by cov.wt.

	# Based on bivariate normal simulations, OGK or GK are recommended robust cov estimates 
	# for smaller data sets (10-30). "MCD" is strongly positively biased for N~10 and highly variable.
	# "mve" and "mcd" are more precise but are consistently biased low by 20-30% from N=10-60. 

	# Note that 'pearson' still allows specification of a scale function, making it possible 
	# (if not very sensible) to use a combination of robust scale function with the 
	# pearson correlation coefficient.

	cov.method <- match.arg(cov.method)
	X <- if(is.null( y )) x else cbind(x, y)
	X <- na.omit( X )
	n.obs <- nrow( X )
	
	if( cov.method %in% eval(formals(cor)$method) ) {
		if(is.null(scalefn)) scalefn <- 
			if(cov.method=="pearson") sd else mad
		if(is.null(locfn)) locfn <- 
			if(cov.method=="pearson") mean else median
		cov.control$method <- cov.method
		cor.mat <- do.call( "cor", c( list(x=X), cov.control ) )
		center <- apply( X, 2, locfn )
		scale <- apply( X, 2, scalefn )
		cov.mat <- cor.mat * outer( scale, scale, "*")
		rv <- structure( list( method=cov.method, cov=cov.mat, cor=cor.mat,  
					center=center, scale=scale, n.obs=n.obs),
					class="cov.dellipse" )
	}

	if( cov.method == "MCD" ) {
		cov.control$cor=TRUE
		m <- do.call( covMcd, c( list( x=X ), cov.control) )
		rv <- structure( list( method=cov.method, cov=m$cov, cor=m$cor,  
				center=m$center, scale=sqrt(diag(m$cov)), n.obs=n.obs),
				class="cov.dellipse" )
	}

	if( cov.method == "OGK" ) {
		if(is.null(cov.control$sigmamu)) {
		    null.fn <- c(is.null(scalefn), is.null(locfn) )
		    cov.control$sigmamu <- if(!any(null.fn)) {
				function(x, mu.too=FALSE, ...) {
					if(mu.too) c(locfn(x), scalefn(x))
					else scalefn(x)
				}
			} else {
				#Consider checking whether scalefun is 
				#present and has a mu.too argument, as for
				#"GK" below
				if(sum(null.fn)==1) warning(
					"One or more of scalefn and locfn NULL: Using scaleTau2 as OGK sigmamu",
					call. = FALSE)
				scaleTau2
			}
		} 
		m <- do.call( covOGK, c( list( X=X ), cov.control) )
		rv <- structure( list( method=cov.method, cov=m$cov, cor=cov2cor(m$cov),  
				center=m$center, scale=sqrt(diag(m$cov)), n.obs=n.obs),
				class="cov.dellipse" )
	}

	if( cov.method %in% c("GK", "gk", "rgk") ) {
		#Variants from (or based on) Gnanadesikan and Kettenring's
		#early suggestions, using covGK from robustbase
		#or (for rgk) a local implementation.
		#GK relies on a scalefn but does not use a mu.too argument
		if(is.null(cov.control$scalefn)) {
		    #scalefn missing from cov.control: 
		    cov.control$scalefn <- if(is.null(scalefn) ) {
		    	scaleTau2  #default
		    } else {
		    	scalefn
		    } 
		} #now have guaranteed control$scalefn
		
		#
		scale <- apply(X, 2, cov.control$scalefn)
		Cov <- if( cov.method == "GK" ) {
			do.call( covGK, c( list( x=X[,1], y=X[,2] ), cov.control) )
		} else if( cov.method == "gk" ) {
			#Cov based on scaled X to avoid overwhelming
			#a small variable with a large one.
			prod(scale) * do.call( covGK, c( list( x=X[,1]/scale[1], y=X[,2]/scale[2] ), cov.control) )
		} else if(cov.method == "rgk") { 
			#This is intended to implement Gnanadesikan and Kettenring's 
			#second covariance estimate based on scaled variables and a 
			#ratio calculation for robust correlation rho.
			#Advantage over "gk" is guaranteed \rho %in% [-1,1]
			Xs <- scale(X, scale=scale)
				#Results are independent of scaling centre; GK do not
				#sweep out location
			Xpm <- cbind(Xs[,1]+Xs[,2], Xs[,1]-Xs[,2] )
			varXpm <- apply(Xpm, 2, cov.control$scalefn)^2
			rho <- -diff(varXpm)/sum(varXpm)
			rho * prod(scale)
		}


		#Get center. Choice of functions here so..
		#Does control$scalefn have a mu.too argument?
		scalefn.has.mu <- 'mu.too' %in% names(formals(cov.control$scalefn)) 
		center <- if(! is.null(locfn) ) {
			apply( X, 2, locfn) 
		} else {
			if(scalefn.has.mu) 
				apply( X, 2, function(x) cov.control$scalefn(x, mu.too=TRUE)[1]) 
			else 
				apply( X, 2, median ) 
		}
		covmat <- diag( scale^2 )
		covmat[1,2] <- covmat[2,1] <- Cov
		rv <- structure( list( method=cov.method, cov=covmat, cor=cov2cor(covmat),  
				center=center, scale=scale, n.obs=n.obs),
				class="cov.dellipse" )
	}


	if( cov.method=="mve" ) {
		cov.control$cor=TRUE
		m <- do.call( cov.mve, c( list( x=X ), cov.control) )
		rv <- structure( list( method=cov.method, cov=m$cov, cor=m$cor,  
				center=m$center, scale=sqrt(diag(m$cov)), n.obs=n.obs),
				class="cov.dellipse" )
	}

	if( cov.method=="mcd" ) {
		cov.control$cor=TRUE
		m <- do.call( cov.mcd, c( list( x=X ), cov.control) )
		rv <- structure( list( method=cov.method, cov=m$cov, cor=m$cor,  
				center=m$center, scale=sqrt(diag(m$cov)), n.obs=n.obs),
				class="cov.dellipse" )
	}

	#rv$center <- rv$center
		#Considering an alias to help folk familiar with cov.rob or cov.wt
		
	return( rv )

} ## End cov.dellipse

