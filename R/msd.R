#Calculates median scaled difference for observations x given sd's s
#if s is a scalar function, it is applied to x to obtain an estimate of s

# This version (2018-02-03) is based on 
# MSD\R\msd_r_development\msd_v7.R
# [local to MSD paper]
#

#This version uses a beta formulation in .pmsd.xnorm to improve numerical stability
#at high N considerably

# It also provides full, if currently slow, support for exact odd N.

# Version 4 added asymptotic q, p, d and quantiles

# Version 4e added trap code 

# Version 5 modified internal odd N integral with pbeta()-integral(...) for more
# stability

# v6 uses a log form which proved stable for very high N and q up to 4
#
# NB: An exception at r==1 is needed as NaNs occur in the log form for r==1

# v7 is a substantial refactoring to 
# a) Include spline interpolation routines
# b) rename existing exact calculations to .xxx.exact
# c) provide new _msd parent functions with a 'method' switch, 
#    defaulting to 'fast'/'interpolate'
# d) Also removes optional scaling arguments; never used in practice 
#    for any distribution function

#Index-based msd
#Still O(n^2) for distance calculation, but a lot faster and lighter on RAM than outer()
msd<-function(x, s=mad , ...) {
        ss <- if(is.function(s)) {
                rep(s(x, ...), length(x))
        } else {
                if(length(s) == 1) {
                        rep(s, length(x))
                } else {
                        s
                }
        }
        
        ss <- ss^2
        
        N <- 1:length(x)
        
        structure( 
        	sapply(N, function(n) median( abs(x[n] - x[-n])/sqrt(ss[n]+ss[-n]) ) ),
        	names=names(x),
        	x=x,
        	s=sqrt(ss),
        	class=c("MSD", "numeric")
        )
}

#
# Basic probability, density and quantile functions for the homoscedastic case
# These are short wrappers allowing choice of method
#

# Note: interpolation functions use N rather than n 
# because the interpolation tables use N as column name

pmsd <- function(q, n, lower.tail=TRUE, method=c('fast', 'exact', 'even', 'asymp'), max.odd=199) {
	method=match.arg(method)
	if(missing(n) && method != "asymp") stop("argument \"n\" is missing, with no default")
	p <- switch(method,
		fast=.pmsd.interp(q=q, N=n),
		exact=.pmsd.exact(q=q, n=n, max.odd=max.odd),
		even=.pmsd.exact(q=q, n=n, max.odd=0),
		asymp=.pmsd.asymp(q=q)	
	)
	if(lower.tail) p else 1-p
}

dmsd <- function(q, n, method=c('fast', 'exact', 'even', 'asymp'), max.odd=199) {
	method=match.arg(method)
	if(missing(n) && method != "asymp") stop("argument \"n\" is missing, with no default")
	switch(method,
		fast=.dmsd.interp(q=q, N=n),
		exact=.dmsd.exact(q=q, n=n, max.odd=max.odd),
		even=.dmsd.exact(q=q, n=n, method="even"),
		asymp=.dmsd.asymp(q=q)	
	)
}

qmsd <- function(p, n, lower.tail=TRUE, method=c('fast', 'exact', 'even', 'asymp'), max.odd=199) {
	method=match.arg(method)
	if(missing(n) && method != "asymp") stop("argument \"n\" is missing, with no default")
	if(!lower.tail) p <- 1-p
	switch(method,
		fast=.qmsd.interp(p=p, N=n),
		exact=.qmsd.exact(p=p, n=n, max.odd=max.odd),
		even=.qmsd.exact(p=p, n=n, max.odd=0),
		asymp=.qmsd.asymp(p=p)	
	)
}


#
# Exact probability, density and quantile calculations using numerical integration
#

.pmsd.xnorm<-function(q, x, n, sd=1, scale=FALSE, method='exact') {
        #P for the median of abs(x-X) (optionally /sqrt(2))
        #based on Mood, Graybill and Boes (1974) pp252ff
        #q can be a vector, n can't
        #NOTE: n is the number of values, not the number of differences.
        
        #Arguments sd and scale retained as used (with default values only) 
        #in calling code
        
        if(scale) sd <- sd/sqrt(2)
        
        pxnorm<-function(q,x,sd=1) ifelse(q>0, pnorm(x+q, 0, sd)-pnorm(x-q, 0, sd), 0) 
        dxnorm<-function(q,x,sd=1) ifelse(q>0, dnorm(x+q, 0, sd)+dnorm(x-q, 0, sd), 0) 
        
        pxFd <- function(t, q, x, sd, r) {
        	FDd <- pxnorm(t, x, sd)
        	FD2d <- pxnorm(2*q-t,x, sd)
        	# Integral[ 2* FDd^(r-1) * ( (1-FDd)^r - (1-FD2d)^r  ) / beta(r,r) ]
        	#rewritten as:
        	#( pbeta(FDd, r, r+1)-  Integral[2* FDd^(r-1) * (1-FD2d)^r / beta(r,r)]dt
        	# with pbeta() taken outside the integral to leave
        	#(2 * ( FDd^(r-1) * (1-FD2d)^r ) / beta(r,r)) * dxnorm(t, x, sd)
        	#which is amenable to log form, so:
        	if(r>1) 2 * exp( log(FDd)*(r-1) + log(1-FD2d)*r - lbeta(r,r)) * dxnorm(t, x, sd)
        	else (2 * (1-FD2d)^r )  * dxnorm(t, x, sd) 
        		#log version fails for r==1 (NaN)
        		#expression for r=1 simplified by beta(1,1) = 1, log(FDd)*(r-1) = 0
        }		
        
        Fy<-rep(0, length(q))
        
        ph<-pxnorm(q,x,sd)
        
        if(method=="exact") {
        	if(n %% 2 == 0) {
			#even-n, odd-median case
			r <- n/2 
			Fy <- pbeta(ph, r, n-r)
				# pbeta(ph, r, n-r+1)*beta(r, n-r-1) is the incomplete beta function
				# as David & Nagaraja define it
				#This form is considerably faster for even quite modest n as it avoids 
				#summation; for n=50 about five times faster than an
				#explicit 'choose' summation for typical comparisons
        	} else {
 			#odd-n, even-median case
       			r <- ( n-1 ) / 2
       			Fy <- pbeta(ph, r, r+1) - 
       				integrate(pxFd, lower=0, upper=q, 
       				  rel.tol = .Machine$double.eps^0.75,  subdivisions = 500L,
       				    q=q, x=x, r=r, sd=sd)$value 
       				#density is zero below 0 so no integral in (-Inf, 0].
       				#Tolerance and subdivisions set low and high
       				#respectively to avoid problems in .qmsd.exact
        	} 
        } else {
        	stop(sprintf("method %s is not currently supported"), method) 
	}
	
        return(Fy)
        
}

.pmsd.exact<-function(q, n, max.odd=199) {
        #P for the median of abs(x-X)/sqrt(2) 
        #based on Mood, Graybill and Boes (1974) pp252ff
        #with additions (odd N) from H A David, H N magaraja, Order Statistics

        #NOTE:  In this implementation, n is the number of values, so msd 
        #       is the median of n-1 differences
        
        #WARNING: this uses a brute force double integration for 
        #         the odd n case. This is slow and can fail for large odd N 
        #         and high probabilities , though N=501 works on most platforms to q=3) 
        #
        #         To handle high N, odd n over max.odd is be replaced with the next even n case.
        #         At n=199, this results in probabilities in error by <1e-4 in the region
        #         q=0.4-0.6, and much less elsewhere. Similarly, quantiles in the same region
        #         differ by no more than 1e-5 on using n=200 instead of n=199
        
 	#Add warning for high max.odd
 	if(max.odd >= 200) warning("max.odd > 200 can result in integration instability for exact probabilities")
 	
        #Former arguments; no longer used but 
        #retained with former default values and with 
        #subsequent code unchanged in case code is reverted
        sd=1
        scale=TRUE 

        #The n max.odd replacement:
        n <- ifelse(n > max.odd & n %% 2 == 1, n+1, n)
        
        L <- max(length(q), length(n), length(sd))
        q<-rep(q, length.out=L)
        n<-rep(n, length.out=L)
        sd<-rep(sd, length.out=L) 
        
        if(scale) sd <- sd/sqrt(2)

        f<-function(x, q, n, sd=1, scale=FALSE, method) {
                L<-length(x)
                q<-rep(q, length.out=L)
                n<-rep(n, length.out=L)
                sd <- rep(sd, length.out=L)
                rv<-rep(0, length(x))
                for(i in 1:length(rv)) rv[i]<-.pmsd.xnorm(q[i], x[i],  n[i], sd[i], scale, method=method)
                return(rv*dnorm(x,0,sd))
        }
        
        p <- rep(0, length(q))
        #q.lower <- switch(method, asym=0.5, 0)
        q.lower=0.0
        for(i in 1:length(q)) {
        	p[i]<-tryCatch(
        		2* integrate(f, lower=q.lower, upper=Inf, 
        	        rel.tol = .Machine$double.eps^0.75, subdivisions = 500L,
                        q=q[i], n=n[i], sd=sd[i], scale=FALSE, method="exact")$value,
                        error=function(e) {
                        	warning(sprintf("Integration failed at q=%.4f, n=%d; NA returned", q[i], n[i]))
                        	NA
                        	}
                        )
                             #Note tolerance; .pmsd.exact and .qmsd.exact are quite inaccurate 
                             #(and .qmsd.exact even unstable) with default integrate() tolerance
        }
        return(p)
        
}

.qmsd.exact<-function(p, n, max.odd=199, ...) {
        #returns quantile(s) for msd using numerical root-finder
        #on .pmsd.exact.
        #p or n may be vectors
        #scale (passed to .pmsd.exact) allows suppression of scaling by sqrt(2)
        
        #Revised 2016-08-31 to use [0,1] scale for root finding and
        #a custom root finder (see .qmsroot)
        
        if(any(p<0) || any(p>1)) stop("p must be in [0,1]")

        #Former arguments; no longer used but 
        #retained with former default values and with 
        #subsequent code unchanged in case code is reverted
        sd=1
        scale=TRUE 

        n <- ifelse(n > max.odd & n %% 2 == 1, n+1, n)
        
        L <- max(length(p), length(n), length(sd))
        p<-rep(p, length.out=L)
        n<-rep(n, length.out=L)
        sd<-rep(sd, length.out=L)

        q<-rep(-1, length(p))
        
        q[p==0] <- 0
        q[ (1-p) < 1e-6 ] <- Inf 
	if(any( 1-p < 1e-6 )) warning("p > 0.999999 will return +Inf")
	
        qq1.upper <- 10/(1+10) #corresponds to  unrealistically high probability
        
        for(i in 1:length(p)) {

                if(q[i]<0) { #Note explicit p==0, p==1 values above
                        #Also note unusual tolerance; .qmsd.exact is  inaccurate with default tolerance
                        # uniroot(froot, interval=c(0,1), tol=.Machine$double.eps^0.75,
                        #        f.lower=-p[i], f.upper=1-p[i],
                        #        p=p[i], n=n[i], sd=sd[i], scale=scale)$root
                	#Custom rootfinder is appreciably more reliable than uniroot for this (!)
                	#Still in tryCatch because rootfinding can hit .pmsd.exact integration error
                	q[i] <- tryCatch( 
                		#.qmsdroot(p=p[i], n=n[i], sd=sd[i], scale=scale, max.odd=max.odd, ...),
                		.qmsdroot(p=p[i], n=n[i], max.odd=max.odd, ...),
                		error=function(e) {
                        		warning(sprintf("qmsd(..., method=\"exact\") failed for p=%e, n=%d; NA returned", p[i], n[i]))
                        		NA
                			}
                		)
                }
        }
        
        # .qmsdroot returns estimates on [0,Inf):
        return(q)
}

#Common problems:
#Error in integrate(pxFd, lower = 0, upper = q, rel.tol = .Machine$double.eps^0.75,  : 
#  the integral is probably divergent
#   Usually arises for odd, high, N and high probabilities. Related to numerical round-off
#   errors in double integration, coupled with very high values of beta(r,r) where r=(N-1)/2.
#
#Error in integrate(pxFd, lower = 0, upper = q, rel.tol = .Machine$double.eps^0.75,  : 
#  non-finite function value
#   Similarly arises in integration for high, odd, N, and for the same reasons.
# In both cases it is almost always sufficient to use the next higher even N; 
# probabilities and quantiles are likely to be essentially identical for most practical purposes 


.qmsdroot <- function(p, ..., tol=1e-6, rel.tol=0.0001, maxit=50, verbose=FALSE, max.odd=199) {
	#This simple bisection search is slower, but substantially more reliable,
	#than the faster search in uniroot, as uniroot seems to move into more
	#extreme regions which cause trouble for the odd-N case. 
	
	to.qq1 <- function(x) x/(1+x)
	from.qq1 <- function(x) x/(1-x)

	lo <- 0
	hi <- 1
	f.lo <- 0-p
	f.hi <- 1-p
	iter<-0
	while(iter <= maxit) {
		iter <- iter+1
		guess <- (lo + hi) / 2 
		f.guess <- .pmsd.exact(guess/(1-guess), ..., max.odd=max.odd) - p
		if(verbose) message(sprintf("iter %2d: estimate=%9.3e+-%6.0e, est.p.error=%7.1e", 
			iter, guess/(1-guess), 0.5 * (hi/(1-hi) - lo/(1-lo)), f.guess))
		if(f.guess < 0 ) {
			lo <- guess 
			f.lo <- f.guess
		} else {
			hi<-guess
			f.hi<-f.guess
		}
		#compound test for convergence
		if(  abs(f.guess) < tol && abs(f.guess) < rel.tol * min(p, (1-p)) ) break 
	}
	if(iter>maxit) warning(sprintf("qmsdroot did not converge in %d iterations; some quantiles may be inaccurate", maxit))

	#final linear interpolation improves results
	qnn.est <- lo - f.lo*(hi-lo)/(f.hi-f.lo)

	return(qnn.est/(1-qnn.est))
}

#
# Density functions
#
# From msd_density_v1.r
#

.dmsd.xnorm <-function(q, x, n, sd=1, scale=FALSE, method=c("exact", "omit", "even"), dx=FALSE) {
        #density given x for the median of abs(x-X) (optionally /sqrt(2))
        #for EVEN n. 
        #based on Mood, Graybill and Boes (1974) pp254

	#setting dx TRUE additionally multiplies by the normal density at x
	
        #scaling parameters sd and scale retained for this helper function,
        #though set values are used in calling code
        
        #NOTE: n is the number of values, not the number of differences.

	method <- match.arg(method)
	
        if(scale) sd <- sd/sqrt(2)
        
        # q must be zero or positive 
        q <- pmax(0,q)
       
        pxnorm<-function(q,x,sd=1) pnorm(x+q, 0, sd)-pnorm(x-q, 0, sd) 
        dxnorm<-function(q,x,sd=1) dnorm(x+q, 0, sd)+dnorm(x-q, 0, sd) 
        
	pxFd <- function(t, q, x, sd, r) {
		FDd <- pxnorm(t, x, sd)
		FD2d <- pxnorm(2*q-t,x, sd)
		2 * FDd^(r-1) * (1-FD2d)^(r-1) * dxnorm(2*q-t, x, sd) * dxnorm(t, x, sd) 
	}

        j <-floor((n+1)/2)   # j is the index for the median (formerly n.med)
                             # exact for odd samples, low for even
                             # Note that for n values there are n-1 differences,
                             # so an even-n case is an odd-median case
        
        #three n-adjustment methods used in development retained
        #for testing; 'omit' no longer an option in exported dmsd
        n <- switch(method,
        	omit=ifelse(n %% 2 ==1, NA, n),
        	even=ifelse(n %% 2 ==1, 2*j, n), #Substitutes the next higher n for the odd cases
        	exact=n 
        )
        
        ph<-pxnorm(q,x,sd)
        
        if(n %%2 == 0) {
	        # beta formulation uses j, n-j rather than j-1, n-j-1 as choose() does
        	fy <- dbeta(ph, j, n-j) * dxnorm(q,x,sd) 
        } else {
		r <- ( n-1 ) / 2
		fy <-  2 * r * integrate(pxFd, lower=0, upper=q, rel.tol = .Machine$double.eps^0.75, 
				q=q, x=x, r=r, sd=sd)$value / beta(r,r)
			# 2*r is required because beta(r,r) = (2r-1)!/(r-1)!(r-1)!
			# so needs multiplication by 2r to get the required (2r)!/(r-1)!(r-1)!
        }
        
        if( dx ) {  
        	 fy * dnorm(x, 0, sd)
        } else {
        	 fy
        }

}


.dmsd.exact<-function(q, n, method=c("exact", "omit", "even"), max.odd=199) {
        #density for the median of abs(x-X)/sqrt(2) 
        #based on Mood, Graybill and Boes (1974) pp252ff

        #NOTE:  In this implementation, n is the number of values, so msd 
        #       is the median of n-1 differences
        
	method <- match.arg(method)

        #Former arguments; no longer used but 
        #retained with former default values and with 
        #subsequent code unchanged in case code is reverted
        sd=1
        scale=TRUE 

        n <- switch(method,
        	omit=ifelse(n %% 2 ==1, NA, n),
        	even=ifelse(n %% 2 ==1, 2*floor((n+1)/2), n), #Substitutes the next higher n for the odd cases
        	exact=ifelse(n > max.odd & n %% 2 == 1, n+1, n)
        )

        if(scale) sd <- sd/sqrt(2)
        
        qns <- cbind(q, n, sd) #ensures same length
        
        f<-function(x, q, n, sd=1, scale=FALSE, method, dx=FALSE) {
                #vectorizes .dmsd.xnorm
                L<-length(x)
                q<-rep(q, length.out=L)
                n<-rep(n, length.out=L)
                sd <- rep(sd, length.out=L)
                rv<-rep(0, length(x))
                for(i in 1:length(rv)) rv[i]<-
                	.dmsd.xnorm(q[i], x[i],  n[i], sd[i], scale, method=method, dx=dx)
                return(rv*dnorm(x,0,sd))
        }

        int.x <- function( qns.i ) {
        	if(any(is.na(qns.i))) {
        		NA
        	} else {
                     2*integrate(f, lower=0, upper=Inf, rel.tol = .Machine$double.eps^0.75, 
                                 q=qns.i[1], n=qns.i[2], sd=qns.i[3], scale=FALSE, method=method)$value
                             #Note odd tolerance; .pmsd.exact and .qmsd.exact are quite inaccurate 
                             #(and .qmsd.exact even unstable) with default integrate() tolerance
        	}
        }
        apply(qns, 1, int.x)
}


#Asymptotic distribution for msd (very large N)

#See asymp_v1.r in Papers/MSD/R/msd_asymp_tests

#Need to find x(=x0) such that the median is q (ie that pxnorm=0.5).
#Then probability is simply 2 * pnorm( x0, 0, 1/sqrt(2) ) -1

#Separate x0 finder to help .dmsd.asymp
.findx<-function(q) {
	pxnorm<-function(q,x,sd=1) pnorm(x+q, 0, sd)-pnorm(x-q, 0, sd) 
	to.qq1 <- function(x) x/(x+1)
	from.qq1 <- function(x) x/(1-x)
	
	fn.findx <- function(x, q)  pxnorm(q, from.qq1(x), sd=1/sqrt(2))-0.5
	
	if(q == Inf) return(1.0) else
	if(q <= qnorm(0.75)/sqrt(2)) NA 
	else uniroot(fn.findx, c(0,1), q=q)$root
}

.pmsd.asymp <- function(q) {
	from.qq1 <- function(x) x/(1-x)
	x0 <- from.qq1( sapply(q, .findx) )
	ifelse(is.na(x0), 0 , 2 * pnorm( x0, sd=1/sqrt(2) ) - 1 )
}

.dmsd.asymp <- function(q, map=c("I", "xx1", "atan", "tanh")) {
	#A rough estimate by simple one-sided numerical differentiation
	#Performs a check for the critical value qnorm(0.75)/sqrt(2)
	#and returns Inf at that point.
	#map transforms the result and can be a function
	
	is.limiting <- function(x, limit) sapply(x, 
		function(xi, limit) isTRUE(all.equal(xi, limit, tolerance = .Machine$double.eps ^ 0.75)), 
		limit=qnorm(0.75)/sqrt(2))
		
	dxnorm<-function(q,x,sd=1) dnorm(x+q, 0, sd)+dnorm(x-q, 0, sd) 
	
	to.qq1 <- function(x) x/(x+1)
	from.qq1 <- function(x) x/(1-x)

	delta.q <- 1e-7
	x0 <- from.qq1( sapply(q, .findx) )
	x0.plus.delta <- from.qq1( sapply(q + delta.q, .findx) )
	dx0.by.dg <- (x0.plus.delta - x0)/delta.q
	dF.by.dx0 <- dxnorm(x0, 0, sd=1/sqrt(2))

	d <- dF.by.dx0 * dx0.by.dg 

	#Handle special cases:
	d[is.na(x0)] <- 0
	d[is.limiting(q)] <- Inf
	
	if(is.function(map) ) 
		map(d) 
	else
		switch(match.arg(map),
			I=d,
			xx1=ifelse(is.finite(d), d/(1+d), 1),
			atan=atan(d),
			tanh=tanh(d)
		)
}

.qmsd.asymp <- function(p) {
	to.qq1 <- function(x) x/(x+1)
	from.qq1 <- function(x) x/(1-x)
	fn.qq1 <- function(qq1, p) .pmsd.asymp( from.qq1(qq1) ) - p
	find.qq1 <- function(p) {
		if(isTRUE(all.equal(0, p))) {
			to.qq1(qnorm(0.75)/sqrt(2)) 
		} else {
			uniroot(fn.qq1, c(0,1), p=p, f.upper=1-p, f.lower=0-p)$root
		}
	}
	qq1 <- sapply(p, find.qq1)
	
	from.qq1(qq1)
}


#
# Fast interpolation functions (from msd_interp_v3.r)
#

# Developer Warning: interpolation functions use N 
# rather than n - mainly because the interpolation tables 
# use N as column name

.msd.ave <- function (x, by, FUN = mean, ...) {
    #Special version of ave allowing ... as args to FUN
    if (missing(by)) 
        x[] <- FUN(x, ...)
    else {
        g <- interaction(by)
        split(x, g) <- lapply(split(x, g), FUN, ...)
    }
    x
}

.quick.cubic <- function(x, y, newx, degree) {
   # In general, this is numerically poorer than using
   # an orthogonal polynomial basis. In this case, though,
   # the X matrix is never very poorly conditioned as 
   # NN1 is in [0,1] and the intervals on NN1 in ptable are 
   # reasonably chosen. And this is _much_ faster (roughly 60x) than 
   # multiple lm(y~poly(x,degree)) calls, or using poly() (8x slower)
	if(degree == 3) {
		m <- cbind(rep(1,4), x, x*x, x^3, deparse.level=0)
		coef <- solve(m,y)
		coef[1] + coef[2]*newx + coef[3]*newx^2 + coef[4]*newx^3
	} else if(degree == 2 ){
		m <- cbind(rep(1,3), x, x*x, deparse.level=0)
		coef <- solve(m,y)
		coef[1] + coef[2]*newx + coef[3]*newx^2
	} else stop("Degree must be 3 or 4")
}	

.pmsd.interp.spline <- function(N, verbose=FALSE) {
	#Generates a monotonic spline on q/(q+1) for
	#approximate, fast probability estimation for any N.

	#Requires scalar, integer N (checked by .pmsd.interp)
	to.NN1 <- function(x) x/(1+x)
	from.NN1 <- function(x) x/(1-x)
	
	#Locate the right table (odd or even)
	#NB: .msd.prob.tables is loaded on package loading;
	#    R's assignment means that no new copy is made
	ptable <- if( N %% 2 == 1 ) { 
		#N is odd
		.msd.prob.tables$odd	
	} else {
		#N is even
		.msd.prob.tables$even	
	}
	
	if ( N %in% attr(ptable, "N" ) ) {
		if(verbose) cat(sprintf("Found N=%d in table %s\n", N, if( N %% 2 == 1 ) "odd" else "even" ))
		#Probabilities for exactly N already exist in ptable:
		#Form spline on the N == N subset of ptable 
		with(ptable[ptable$N == N,], splinefun(qq1, p, method='h'))
	} else {
		if(verbose) cat(sprintf("Did not find N=%d. Interpolating on table %s\n", N, if( N %% 2 == 1 ) "odd" else "even" ))
		#Quantiles for N must be interpolated before
		#forming q-p spline
		ptable.N <- attr(ptable, "N" ) #For convenience
		
		#Find nearest 4 (or 3, if N >100000) values of N in ptable
		#NB: This can include Inf; but interpolation for p(q|N) will be on (N/1+N)
		nearest.N <- c( rev( ptable.N[ ptable.N < N ] )[2:1],
				ptable.N[ ptable.N > N] [1:2] ) 
		nearest.N <- nearest.N[ ! is.na(nearest.N) ] #Drops any trailing NA

		#Construct a table for required N, by aggregating over q using interpolation
		interp.N <- function(rows, ptable, N, degree) {
			if(any(ptable[rows,]$p < 2 * .Machine$double.xmin ) ) {
				#p effectively 0; don't interpolate in this region
				0.0
			} else {
			   #Interpolate on log(p) to avoid negative predictions 
			   #at very low p
			   
			   #This _can_ be done  using lm and poly: 
				#lm.qq1 <- lm(log(p) ~ poly(NN1, degree), ptable[rows,]) 
				#exp(predict(lm.qq1, newdata=list(NN1=to.NN1(N))))
			   #but that is slow due to multiple lm calls)

			   #Much quicker to use explicit cubic solution
			   new.NN1 <- to.NN1(N)
			   with(ptable[rows,], 
			   	exp( .quick.cubic(x=NN1, y=log(p), 
			   		newx=new.NN1, degree=degree) ) )
			}
		}
		
		#Need a subset of ptable for nearest.N
		ptable.nearest.N <- ptable[ ptable$N %in% nearest.N, ]
		
		#aggregate over row indices by qq1
		poly.degree <- length(nearest.N) - 1 #Simple way of controlling poly fit
		ptable.N <- aggregate(1:nrow(ptable.nearest.N), 
					by=list(qq1=ptable.nearest.N$qq1), 
					FUN=interp.N, 
					ptable=ptable.nearest.N, N=N, degree=poly.degree)
			#Generates a data frame with columns qq1 and x, 
			#where x is interpolated probability
		with(ptable.N, splinefun(qq1, x, method='h'))
		
	}
}

.msd.interp <- function(q, N, deriv=0, verbose=FALSE) {
	#Function to obtain interpolated probabilities or densities from tabulated 
	#MSD quantiles using spline interpolation
	#WARNING: The spline is constructed on q/(1+q), so
	#         densities are also wrt q/(1+q)
	
	#Setting deriv=0 returns probabilities; deriv=1 returns densities

	#N must be integers >= 2
	is.wholenumber <-
	         function(x, tol = .Machine$double.eps^0.75)  abs(x - round(x)) < tol
	
	if( any( ! is.wholenumber( N ) ) ) stop("N must be integers")
	if( any( N < 2 ) ) stop("N must be >= 2")

	to.qq1 <- function(x) x/(1+x)
	from.qq1 <- function(x) x/(1-x)

	pmsd.interpN <- function(iq, qN, deriv, verbose=FALSE) {
		#Apply spline interpretation to the 
		#data subset indicated by indices iq
		#N[iq] are guaranted equal by tapply()
		if(verbose) cat(sprintf("N[iq][1]=%f\n",N[iq][1]))
		with(qN, 
			.pmsd.interp.spline(N[iq][1], verbose=verbose)(to.qq1( q[iq]), 
				deriv=deriv) )
	}
	
	if(length(N1<-unique(N)) == 1) {
		#Single value of N can be dealt with slightly faster
		#without extending N or q
		.pmsd.interp.spline(N1, verbose=verbose)(to.qq1(q), deriv=deriv)
	} else {
		#Multiple values of N
		#Interpolation is for each N		
		#Match lengths for N and q and
		#form data frame for aggregate()
		L <- max(length(q), length(N))
		qN <- data.frame(q=rep(q, length.out=L), N=rep(N, length.out=L))
		#Modified ave applies pmsd.interpN over qN for all unique N
		#Uses indices rather than q itself as corresponding N must be passed
		.msd.ave(1:L, by=list(N=N), FUN=pmsd.interpN, qN=qN, 
			deriv=deriv, verbose=verbose)
	}
}

.pmsd.interp <- function(q, N) {
	.msd.interp(q, N, deriv=0)
}

.dmsd.interp <- function(q, N) {
	d.qq1 <- .msd.interp(q, N, deriv=1)
	#transform for q/(1+q) ...
	t <- 1/(1 + q) - q/(1 + q)^2
	t[is.na(t)] <- 0
	d.qq1 * t
}

.qmsd.interp <- function(p, N, verbose=FALSE) {
	#Function to obtain quantiles from tabulated 
	#MSD probabilities using spline interpolation and
	#root-finding
	#As above, the spline is constructed on q/(q+1), denoted qq1

	#N must be integers >= 2
	is.wholenumber <-
	         function(x, tol = .Machine$double.eps^0.75)  abs(x - round(x)) < tol
	
	if( any( ! is.wholenumber( N ) ) ) stop("N must be integers")
	if( any( N < 2 ) ) stop("N must be >= 2")

	to.qq1 <- function(x) x/(1+x)
	from.qq1 <- function(x) x/(1-x)

	#Objective function for root-finding
	#declared here to avoid repeated parsing in
	#qmsd.uniroot
	qmsd.objective <- function(qq1, p1, qsfun) p1-qsfun(qq1)

	#Single-p root-finding function
	qmsd.uniroot <- function(p1, qsfun) {
		#p1 is a scalar probability
		#qsfun is a spline function returning probability for qq1
		uniroot(qmsd.objective, interval=c(0,1), p1=p1, qsfun=qsfun,
			tol = .Machine$double.eps^0.75)$root
	}
	
	qmsd.interpN <- function(iq, pN, verbose=FALSE) {
		#Apply spline interpolation to the 
		#data subset indicated by indices iq
		#N[iq] are guaranted equal by .msd.ave()
		if(verbose) cat(sprintf("N[iq][1]=%f\n",N[iq][1]))
		#Get the spline function:
		qsfun <- with(pN, 
			.pmsd.interp.spline(N[iq][1], verbose=verbose) )
		#Find roots corresponding to pN[iq]
		qp <- sapply(pN[iq], qmsd.uniroot, qsfun=qsfun)
		from.qq1(qp)
	}
	
	if(length(N1<-unique(N)) == 1) {
		#Single value of N can be dealt with slightly faster
		#without extending N or q
		qsfun <- .pmsd.interp.spline(N1, verbose=verbose)
		qp <- sapply(p, qmsd.uniroot, qsfun=qsfun)
		from.qq1(qp)
	} else {
		#Multiple values of N
		#Interpolation is for each N		
		#Match lengths for N and q and
		#form data frame for aggregate()
		L <- max(length(q), length(N))
		pN <- data.frame(p=rep(p, length.out=L), N=rep(N, length.out=L))
		#Modified ave applies qmsd.interpN over pN for all unique N
		#Uses indices rather than p itself as corresponding N must be passed
		.msd.ave(1:L, by=list(N=N), FUN=qmsd.interpN, pN=pN, verbose=verbose)
	}
}
