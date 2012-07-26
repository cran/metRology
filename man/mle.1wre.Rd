\name{mle.1wre}
\Rdversion{1.1}
\alias{mle.1wre}
\title{
Vangel-Rukhin Maximum Likelihood Estimate
}
\description{
     Calculate a weighted mean, between-group standard deviation and
     standard error on the weighted mean using the Maximum likelihood 
     algorithm of Vangel-Rukhin.
}
\usage{
mle.1wre(x, s2, n, init.mu = mean(x), init.sigma2 = var(x), labels = c(1:length(x)), max.iter = 200, tol = .Machine$double.eps^0.5, trace = FALSE)
}

\arguments{
  \item{x}{
  	numeric vector of the sample mean values of each group
   
}
  \item{s2}{
	numeric vector of the sample variances of each group
}
  \item{n}{
	integer vector of sample size of each group
}
  \item{init.mu}{
	numeric initial value for the mean
}
  \item{init.sigma2}{
	numeric initial value for the between-group component of variance
}
  \item{labels}{
	vector of group names. Coerced to character on use.
}
  \item{max.iter}{
	numeric maximum number of iterations
}
  \item{tol}{
	numeric tolerance; iteration stops when the relative step size drops 
	below 'tol'
}
  \item{trace}{
	when TRUE shows the sequence of intermediate results
}
}
\details{
     The Vangel-Rukhin MLE algorithm finds the between-method variance by
     iteratively solving the equation relating the weighted mean to the
     weighting factor applied. The weighting factor is the inverse of
     the sum of the standard error in 'x' and the between-method
     variance, scaled by the between-method variance.

     For the default method, 's2'  is interpreted as a vector of
     sample variances. 'x' is interpreted as a vector of sample means
     and the algorithm is applied to the corresponding group means,
     variances, and sample sizes.

     The Vangel-Rukhin MLE algorithm shows an improvement in the number
     of iterations required to converge over the classical MLE based on
     the Score equations.

     The function mle.1wre implements the MLE for the one way random effects 
     based on the Fisher scoring equations and is provided for comparison 
     purpose only.
}
\value{
	\code{mle.1wre} returns an object of class "summary.mle.1wre" which contains
	the following fields:
	
	    \item{mu}{the estimated mean
	}
	    \item{var.mu}{the variance associated with the estimated mean
	}
	    \item{sigma2}{the estimated between variance component
	}
	    \item{llh}{the log likelihood of the estimates
	}
	    \item{tot.iter}{the total number of iterations ran
	}
	    \item{cur.rel.abs.error}{the current relative absolute error reached
	}
	    \item{sigmai2}{a vector with the estimates of the within variance components
	}
}
\references{
     Vangel, M. G. and Rukhin, A. L. (1999), Biometrics, Vol 55, No. 1 pp 129-136

     Searle, S. R., Cassella, G., and McCulloch, C. E. (1992). Variance Components. New York: Wiley.
}
\author{
     H. Gasca-Aragon
}

\seealso{
\code{\link{vr.mle}}, \code{\link{loc.est-class}}
}
\examples{

     ##===================================================================
     ## the dietary fiber in apples example in the Vangel and Rukhin paper
     ##===================================================================

     m1 <- c(12.46, 13.035, 12.44, 12.87, 13.42, 12.08, 13.18, 14.335, 12.23)
     s1 <- c(0.028, 0.233, 0.325, 0.071, 0.339, 0.325, 0.099, 0.064, 0.212)
     n1 <- c(2, 2, 2, 2, 2, 2, 2, 2, 2)

     mle.1wre(m1, s1^2, n1, tol=1e-6)


     # output:
     # 12.90585
     # 0.2234490
     # 0.4262122
     # 12.46790 13.34380
     # 6
 
}
\keyword{ univar }
