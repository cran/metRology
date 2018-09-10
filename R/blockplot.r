#
# blockplot  -- Produces a labelled histogram
#

# v3 - adds group handling.

# v4 - supports plot=FALSE and includes level labels

# v5 separates creation and plotting

## Omitted for now: object name does not parse correctly in bxp()
##
## plot.blockplot <- function(x, ...) {
## 	bkp(x, ...)
## }

bplot <- function(x, ...) {
	UseMethod("blockplot")
}

blockplot <- function(x, ...) {
	UseMethod("blockplot")
}

nclass.23 <- function(x) ceiling(length(x)^(2/3))

blockplot.default <- function(x, breaks="23", labels=paste(1:length(x)), groups=NA, 
	xlim = NULL, ylim = NULL, main = NULL, xlab = NULL, ylab="Frequency", 
	grp.labs=FALSE,  include.lowest = TRUE, right = TRUE, nclass = NULL,  
	plot = TRUE,  add=FALSE, ...) {

	# Constructs the breaks and locations required to produce a block plot from 
	# numeric data x
	
	# groups allows multiple blockplots, vertically separated
	
	# uline specifies line under blocks. 
	#       TRUE: full width line
	#      FALSE: No line
	#    numeric: Number of box widths either side
	#             of data. Can be length 2 vector for
	#             left and right. Symmetric if length 1
	#
	
	# add. Add to an existing plot. TRUE means add. new plot if FALSE.
	
	#axes  Logical; specifies whether axes are drawn. Can be length >=2 to control 
	#               axis generation on each side. 
	
	# offset simply shifts everything up (or, if negative, down) by offset; useful for
	#	constructing with add.
	
	# grp.spacing is a minimum spacing between groups stacked vertically. 
	# The actual space allocated is the maximum height for any group plus grp.spacing
	
	# grp.at specifies the vertical location for group baselines; if present, 
	#	it overrides grp.spacing.
	
  
	if (!is.numeric(x)) 
	        stop("'x' must be numeric")
    	
    	xname <- paste(deparse(substitute(x)), collapse = "\n")
    	gname <- if(!is.na(groups[1])) {
    			if( !is.null(gn<-attr(groups, "gname")) ) gn
    			else deparse(substitute(groups)) 
    		} else NA
	
	if(is.null(main)) main <- paste("Block plot of", xname, if(!is.na(groups[1])) paste("by", deparse(substitute(groups))) )
	if(is.null(xlab)) xlab <- xname
	
	if(is.na(groups[1])) {
		g <- gl(1, length(x))
	} else {
		g <- factor(groups)
	}
	#Append the group name for later
	g <- structure(g, gname=gname)
	
	# Drop any missing/non-finite values from x, labels and groups
	keep <- is.finite(x)
	x0 <- x[keep]
	
	#match missing values in g to those in x
	g[!keep] <- NA 
	ng <- nlevels(g)

	
	#Get breaks from hist() code, modified for "23"
	n <- length(x0)
	n <- as.integer(n)
	if (is.na(n)) 
		stop("invalid length(x)")
	
	use.br <- !missing(breaks)
	if (use.br) {
		if (!missing(nclass)) 
		    warning("'nclass' not used when 'breaks' is specified")
	}
	else if (!is.null(nclass) && length(nclass) == 1L) 
	        breaks <- nclass
	use.br <- use.br && (nB <- length(breaks)) > 1L
	if (use.br) 
		breaks <- sort(breaks)
	else {
		if (!include.lowest) {
		    include.lowest <- TRUE
		    warning("'include.lowest' ignored as 'breaks' is not a vector")
		}
		if (is.character(breaks)) {
		    breaks <- match.arg(tolower(breaks), c("sturges", 
			"fd", "freedman-diaconis", "scott", "23"))
		    breaks <- switch(breaks, 
		    		sturges = nclass.Sturges(x0),
		    		`freedman-diaconis` = ,
		    		fd = nclass.FD(x0), 
		    		scott = nclass.scott(x0),
		    		"23" = nclass.23(x0),
				stop("unknown 'breaks' algorithm"))
		}
		else if (is.function(breaks)) {
		    breaks <- breaks(x0)
		}

		if (length(breaks) == 1) {
		    if (!is.numeric(breaks) || !is.finite(breaks) || 
			breaks < 1L) 
			stop("invalid number of 'breaks'")
		    breaks <- pretty(range(x0), n = breaks, min.n = 1)
		    nB <- length(breaks)
		    if (nB <= 1) 
			stop(gettextf("block plot error: breaks=%s", 
			  format(breaks)), domain = NA)
		} 
		else {
			#Length(breaks) >1
		    if (!is.numeric(breaks) || length(breaks) <= 1) 
			stop(gettextf("Invalid breakpoints produced by 'breaks(x)': %s", 
			  format(breaks)), domain = NA)
		    breaks <- sort(breaks)
		    nB <- length(breaks)
		    use.br <- TRUE
		}
	}
	nB <- as.integer(nB)
	if (is.na(nB)) stop("invalid length(breaks)")
	h <- diff(breaks)
	equidist <- all(abs(h-h[1]) <= .Machine$double.eps^(0.25))
	if(!equidist) stop("Breaks must be equal-spaced in a block plot")
	blockwidth <- h[1]
	mids <- 0.5 * (breaks[-1L] + breaks[-nB])
	
	#Check for missed x
	if(any(x0 > max(breaks)) || any(x0<min(breaks)) ) {
        	stop("Some values outside range of 'breaks'")
    	}
	
	cx <- cut(x, breaks=breaks, include.lowest=include.lowest, right=right)
    		#Note x, not x0. This retains length(x) for bxp consistency
    		
	x.left <- breaks[cx]
	x.mid <- mids[cx]
	x.height <- ave(x, cx, g, FUN=function(x) (1:length(x))[order(x, na.last=T)] -1 )
		#Note g as grouping factor here; sets h.height within each group, not globally
	
	rv <- list(x=x, groups=if(ng>1) g else NA, breaks=breaks, labels=labels,
		x.left=x.left, x.height=x.height, x.mid=x.mid)

	if(plot) bkp(rv, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab=ylab, 
		grp.labs=grp.labs,  add=add, ... )
	
	return(invisible( structure(rv, class=c("blockplot", "list") ) ))
#End blockplot
}

blockplot.formula <- function(x, data = NULL, ..., subset, main=NULL, xlab=NULL) {
    if (missing(x) || (length(x) != 3L)) 
        stop("Formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$main <- NULL
    m$xlab <- NULL
    m$formula <- m$x
    m$x <- NULL
    m$na.action <- na.pass
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    
    xname <- names(mf)[response]
    gname <- attr(attr(mf, "terms"), "term.labels")
    
    if(is.null(main)) main <- paste("Block plot of", xname, "by", paste(gname, collapse=", "))
    if(is.null(xlab)) xlab <- xname
    
    #Create a grouping factor with a 'gname' attribute for blockplot.default and bkp to find:
    g <- structure(do.call('interaction', c(mf[-response], drop=TRUE, sep=':')), gname=gname)
    	
    blockplot(mf[[response]], groups=g, main=main, xlab=xlab, ...)

}


bkp <- function(x, labels=x$labels, xlim = NULL, ylim = NULL, 
	main = NULL, xlab = NULL, ylab="Frequency", 
	square=FALSE, add=FALSE, offset=0, grp.spacing=2, grp.at=NA,
	fill = NA, border = NULL, density = NULL, angle = 45, lty=1, lwd=2, 
	label.col=1, cex=NA, adj=c(0.5, 0.4), 
	uline=2, uline.lwd=lwd, uline.lty=1, uline.col=if(!is.null(border)) border else par("fg"), 
	grp.labs=FALSE, grp.pos=1, glab.control=list(), axes = c(TRUE, FALSE), 
	asp=NA, frame.plot=any(axes),  drop.unused=TRUE, unused.label="[Missing]", ...) {

	# groups allows multiple blockplots, vertically separated
	
	# uline specifies line under blocks. 
	#       TRUE: full width line
	#      FALSE: No line
	#    numeric: Number of box widths either side
	#             of data. Can be length 2 vector for
	#             left and right. Symmetric if length 1
	#
	
	# add. Add to an existing plot. TRUE means add. new plot if FALSE.
	
	#axes  Logical; specifies whether axes are drawn. Can be length >=2 to control 
	#               axis generation on each side. 
	
	# offset simply shifts everything up (or, if negative, down) by offset; useful for
	#	constructing with add.
	
	# grp.spacing is a minimum spacing between groups stacked vertically. 
	# The actual space allocated is the maximum height for any group plus grp.spacing
	
	# grp.at specifies the vertical location for group baselines; if present, 
	#	it overrides grp.spacing.
	
	axis.local <- function (side, at = NULL, labels = TRUE, tick = TRUE, line = NA, 
		pos = NA, outer = FALSE, font = NA, lty = "solid", lwd = 1, 
		lwd.ticks = lwd, col = NULL, col.ticks = NULL, hadj = NA, 
		padj = NA, cex.axis=par("cex.axis"), las=par("las"), ...) {
		
		axis(side, at=at, labels=labels, tick=tick, line=line, 
		pos = pos, outer = outer, font = font, lty = lty, lwd = lwd, 
		lwd.ticks = lwd.ticks, col = col, col.ticks = col.ticks, hadj = hadj, 
		padj = padj, cex.axis=cex.axis, las=las)
	}

    	xname <- paste(deparse(substitute(x)), collapse = "\n")
	gname <- attr(x$groups, "gname")
	if(is.null(main)) main <- paste("Block plot of", xname, if(!is.na(x$groups[1])) paste("by", gname) )
	if(is.null(xlab)) xlab <- xname
	
	# Get any group labels if necessary
	if(is.logical(grp.labs[1])) {
		if(grp.labs[1]) {
			grp.labs <- levels(x$groups)
			label.groups <- TRUE
		} else {
			grp.labs <- NA
			label.groups <- FALSE
		}
	} else {
		#labels specified
		label.groups <- TRUE
		grp.labs<-as.character(grp.labs)
		if(length(grp.labs) != nlevels(x$groups)) stop("grp.labs wrong length for levels present in groups")
	} 

	if(label.groups) {
		grp.pos <- rep(grp.pos, length.out=nlevels(x$groups))
	}
	
	if(is.na(x$groups[1])) {
		g <- gl(1, length(x$x))
	} else {
		g <- x$groups
	}
	
	if(drop.unused) {
		g <- droplevels(g) #This also drops any gname attribute, but that is no longer used.
		if(label.groups) {
			#Remove unused labels and positions
			which.unused <- levels(x$groups) %in% levels(g) 
			grp.labs <- grp.labs[which.unused]
			grp.pos <- grp.pos[which.unused]
		}
	} else {
		if(label.groups) {
			#mark any 'missing' group labels
			missing.groups <- which( ! levels(x$groups) %in% levels(droplevels(g)))
			grp.labs[ missing.groups ] <- paste(grp.labs[ missing.groups ], unused.label)
			grp.pos[ missing.groups ] <- 1
		}
	}
	ng <- nlevels(g)
	
	## Drop any missing/non-finite values from x, labels and groups
	keep <- is.finite(x$x)
	x0 <- x$x[keep]
	g0 <- droplevels(g[keep]) #droplevels in case missing values dropped a complete group
	labels0 <- labels[keep]

	#... and lots of box and label colours to recycle explicitly 
	recycle.and.trim <- function(p) if(length(p)>0) rep(p, length.out=length(x$x.mid)) else p

	fill <- recycle.and.trim( fill )
	border  <- recycle.and.trim( border )
	density  <- recycle.and.trim( density )
	angle <- recycle.and.trim( angle )
	lty <- recycle.and.trim( lty )
	lwd  <- recycle.and.trim( lwd )
	label.col <- recycle.and.trim( label.col )

	if(!is.logical(uline)) uline <- rep(uline, length.out=2)
	
	#Get and check breaks from x
	h <- diff(x$breaks)
	equidist <- all(abs(h-h[1]) <= .Machine$double.eps^(0.25))
	if(!equidist) stop("Breaks must be equal-spaced in a block plot")
	blockwidth <- h[1]
	
    	if(is.na(grp.at[1])) {
		grp.at <- (0:max(0, ng-1)) * max(x$x.height + grp.spacing + 1, na.rm = TRUE)
					        # +1 because x.height is the bottom of a block
	}
	
	#Get heights for 'empty' plotting and ylim
	xh.g <- x$x.height+grp.at[g]
	
	if(is.null(ylim)) ylim<- c(0, max(xh.g, na.rm = TRUE)+offset+1.1) 
	if(is.null(xlim)) xlim<-range( x$breaks ) + c(-1,1)*pmax(0, uline)
	
	if(square) asp <- blockwidth
		#because block height == 1.0
	
	
	if(label.groups) {
		grp.pos <- rep(grp.pos, length.out=ng)
	}
	
	#Commence plot:
	if(!add) plot(x$x.mid, xh.g, type='n', ylim=ylim, xlim=xlim, 
		asp=asp, xlab=xlab, ylab=ylab, main=main,  axes=FALSE, ...) 

	#Set label cex to fit inside boxes
	if(is.na(cex)) {
		labelwidth <- max(strwidth(labels, cex=1))
		labelheight <- max(strheight(labels, cex=1))
		cex <- 0.8 * min(blockwidth/labelwidth, 1/labelheight)
	}

	if(length(axes) == 1) axes <- rep(axes, length.out=2)

	for(ig in 1:ng) {
		#For loop easier, for once, as we need the group index
		which.x <- which(g == levels(g)[ig])
		grp.empty <- length( which( !is.na( x[which.x] ) ) ) == 0
		if( !grp.empty ){
			#Plot the boxes
			xl <- x$x.left[which.x]
			xh <- x$x.height[which.x]
			y <- xh + grp.at[ig]+offset
			rect(xl, y, xl+blockwidth, y+1, 
				col=fill[which.x], border=border[which.x], 
				angle=angle[which.x], density=density[which.x], 
				lwd=lwd[which.x], lty=lty[which.x])
			text(x$x.mid[which.x], y+0.5, labels[which.x], col=label.col[which.x], cex=cex, adj)

			#Baseline
			if(is.logical(uline) && uline ) { #uline not tested if not logical
				abline(h=grp.at[ig]+offset, col=uline.col, lwd=uline.lwd, lty=uline.lty)
			} else if(is.numeric(uline)) {
				segments(min(xl, na.rm=TRUE)-uline[1]*blockwidth, grp.at[ig]+offset, 
					max(xl, na.rm=TRUE)+blockwidth*(1+uline[2]), grp.at[ig]+offset, 
					col=uline.col, lwd=uline.lwd, lty=uline.lty)
			} 

			#Left/right axes:
			bp.top <- max(xh, na.rm = TRUE)+1
			v.axis.at <- pretty(0:bp.top, n=min(5, bp.top))
			v.axis.at <- v.axis.at[ v.axis.at <= max(xh, na.rm = TRUE) + 1 ]
			if(axes[2]) {
				axis.local(2, at=grp.at[ig]+offset+v.axis.at, labels=paste(v.axis.at), ...)
			}
			if(length(axes)>3 && axes[4]) {
				axis.local(4, at=grp.at[ig]+offset+v.axis.at, labels=paste(v.axis.at), ...)
			}
		} 
		
		#Group labels
		if(label.groups) {
			if( !grp.empty ){
				grp.mid.x <- sum( range(x$x.mid[which.x], na.rm = TRUE) )/2
				grp.mid.y <- 0.5 + sum( range(y, na.rm = TRUE) )/2
				glab.x <- switch(grp.pos[ig],
					grp.mid.x,
					min(x$x.left[which.x]) - blockwidth/2,
					grp.mid.x,
					max(x$x.left[which.x]) + blockwidth*1.5			
				)
				glab.y <- switch(grp.pos[ig],
					min(y, na.rm = TRUE), grp.mid.y, max(y, na.rm = TRUE)+1, grp.mid.y)
				glab.args <- c(list(x=glab.x, y=glab.y, labels=grp.labs[ig]), glab.control,
					if(is.null(glab.control$pos) && is.null(glab.control$adj)) list(pos=grp.pos[ig]))
			} else {
				glab.x <- sum( range(x$x.mid, na.rm = TRUE) )/2
				glab.y <- offset + grp.at[ig] + 0.5 + sum( range(x$x.height, na.rm = TRUE) )/2
				glab.args <- c(list(x=glab.x, y=glab.y, labels=grp.labs[ig]), glab.control)
				#and disable pos and adj for this 'missing' label
				glab.args$pos <- NULL
				glab.args$adj <- NULL
			}
			
			do.call("text", c(glab.args))
		}


		if(!add) {
			if(axes[1]) axis.local(1, ...)
			if(length(axes)>2 && axes[3]) axis.local(3, ...)
		}
		if(frame.plot) box()
	}
	
	return(invisible(c(x, list(grp.at=grp.at, blockwidth=blockwidth))))
#End bbkp
}

