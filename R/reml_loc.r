#Implements restricted maximum likelihood fit for
#a set of means and variances (or sd's)
#This implementation parallels mpaule in metRology
#in allowing a vector of raw data, with grouping factor

.neg.ll.reml <- function(mu, lambda, x, u, REML=TRUE, exact=FALSE, loglambda=FALSE) {
	#Returns the negative log-restricted likelihood or (if !REML) just the maxlik
	#Allows log-lambda specification for alternative optimisation methods.
	
	if(loglambda) lambda <- exp(lambda)
	
	vars <- lambda + u^2
	resids <- x - mu #Weighting occurs at sumsq stage (below)
	
	rv <- if(exact) {
		k <- length(x)
		p <- 1
		W <- diag(1/vars)
		X <- as.matrix(rep(1,k))
		#"exact" restricted likelihood (as used by metafor rma.uni)
		ll.reml <- if(REML) {
			-1/2 * (k - p) * log(2 * base::pi) + 
			(if(REML) 1/2 * determinant(crossprod(X), logarithm = TRUE)$modulus else 0) - 
			1/2 * sum(log(vars)) - 
			1/2 * determinant(crossprod(X, W) %*% X, logarithm = TRUE)$modulus - 
			1/2 * sum((resids^2)/vars) #checked exact against crossprod form
		} else {
			 -1/2 * (k) * log(2 * base::pi) - 1/2 * sum(log(vars)) 
			   - 1/2 * (sum((resids^2)/vars))
		}
	} else {
		#simplified - omits constant term 
		(-1/2) * ( sum((resids^2)/vars)  + sum(log(vars)) + 
			if(REML) log(sum(1/vars)) else 0 )
	}
	-rv  #return negated value for minimisation
}	

reml.loc<-function(x, ..., na.rm=FALSE) {
	UseMethod("reml.loc")
}

reml.loc.default<-function(x, s, n=NULL, groups=NULL, na.rm=FALSE, 
		tol=.Machine$double.eps^0.5, REML=TRUE, ...) {
	#x is a vector of (mean) observations
	#s is a vector of standard errors/standard uncertainties
	# or (if n is given) standard deviations
	#n is an optional mumber of observations per group
	#groups is a grouping factor, used to group x
	#tol is a convergence criterion
	
	#This uses a direct implementation using 1-d optimisation
	
	#Prepare to remove NA's 
	if(!is.null(n)) {
		n.in <- n
		n <- rep(n, length(x)) #Recycled if present
	}
	x.in <- x 
	s.in <- if(missing(s)) NA else s
	g.in <- if(is.null(groups)) NA else groups
	
	if(na.rm) {
		na.xs <- if(missing(s)) 
				is.na(x) 
			else 
				(is.na(x) | is.na(s))
		x <- x[!na.xs]
		if(!missing(s)) s <- s[!na.xs]
		if(!is.null(groups)) groups <- groups[!na.xs]
		if(!is.null(n)) n <- n[!na.xs]
	}
	
	if(!is.null(groups)) {
		if(!missing(s) || !is.null(n)) 
			warning("Using groups to calculate stderr: ignoring s and n")

		count<-function(x) length( x )
		groups <- factor(groups)
		x <- as.vector(tapply(xi<-x, groups, mean, na.rm=na.rm))
		s <- as.vector(tapply(xi, groups, sd, na.rm=na.rm))
		n <- as.vector(tapply(xi, groups, count))
		std.err <- s/sqrt(n)
	} else {
		if(missing(s)) stop("Either s or groups must be specified")
		std.err <- if(!is.null(n)) s/sqrt(n) else s
	}
	
	#Now have a vector of mean values x and stderr std.err
	reml1 <- function(lambda, x, u, REML=TRUE, exact=FALSE, log=FALSE) {
		#Calculates weighted mean and then calls .neg.ll.reml
		if(log) lambda <- exp(lambda) 
		vars <- lambda + u^2
		w <- 1/vars 
		mu <- sum(x*w)/sum(w)
		.neg.ll.reml(mu=mu, lambda=lambda, x=x, u=u, exact=exact, 
			REML=REML, loglambda=FALSE) #conversion above makes loglambda always FALSE
	}	
	lo <- 0.0
	hi <- max(c(var(x), std.err^2))
       	rv <- optimize(reml1, maximum=FALSE, tol=tol, lower=lo, upper=hi, x=x, u=std.err, REML=REML)
	lambda <- unname(rv$minimum)
        w <- 1 / (lambda + std.err^2)
        mu.reml <- unname(sum( x * w ) / sum (w))
	s.reml <- 1/sum(w)
	s.adj <- sqrt(std.err^2 + lambda)


	rv <- .construct.loc.est( x=mu.reml, u=s.reml, df=length(x)-1, 
		xi=x, ui=s, u.eff=s.adj, w=rep(1, length(x)), method=if(REML) "REML" else "ML", 
		method.details=list(mu=mu.reml, s=s.reml, tau=sqrt(lambda), REML=REML) )
	
	return(rv)
	
}

