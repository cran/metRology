\name{algS}
\alias{algS}
\title{
\sQuote{Algorithm S} - robust estimate of pooled standard deviation
}
\description{
\sQuote{Algorithm S} calculates a robust estimate of pooled standard deviation 
from a set of standard deviations
}
\usage{
	algS(s, degfree, na.rm = FALSE, prob.eta = 0.9, 
		is.range = FALSE, tol = .Machine$double.eps^0.25, 
		maxiter = 25, verbose = FALSE)
}
%% $
\arguments{
  \item{s}{ A vector of standard deviations or, if \code{is.range} is \code{TRUE},
  		ranges.
}
  \item{degfree}{Scalar number of degrees of freedom associated with all 
  	values in \code{s}. If a vector is supplied,  \code{median(degfree)} will be used.
}
  \item{na.rm}{a logical value indicating whether 'NA' values should be
          stripped before the computation proceeds.
}
  \item{prob.eta}{prob.eta is set to specify the lower tail area of the chi-squared 
	distribution used as a cut-off.
}
  \item{is.range}{if is.range is TRUE, s is interpreted as a vector of positive differences of 
	duplcate observations and degfree is set to 1
}
  \item{tol}{Convergence tolerance Iteration continues until the relative 
  		change in estimated pooled sd drops below \code{tol}.
}
  \item{maxiter}{Maximum number of iterations permitted.
}
  \item{verbose}{Controls information displayed during iteration; see Details.
}
}
\details{
Algorithm S is suggested by ISO 5725-5:1998 as a robust estimator of 
pooled standard deviation \eqn{s_{pool}}{s.pool} from standard deviations of 
groups of size \eqn{\nu+1}{degfree}. 

The algorithm calculates a \sQuote{limit factor}, \eqn{\eta}{eta}, set to 
\code{qchisq(prob.eta, degfree)}. Following an initial estimate of 
\eqn{s_{pool}}{s.pool} as \code{median(s)}, the standard deviations \eqn{s_i}{s{i}}
are replaced with \eqn{w_i=min(\eta*s_{pool}, s_i)}{w[i]=min(eta*s.pool, s[i])} 
and an updated value for \eqn{s_{pool}}{s.pool} calculated as
	\deqn{\xi*\sqrt{\frac{\sum_{i=1}^{p} (w_i)^2}{p}}}{xi*sqrt(sum(w)^2)/p}
where \eqn{p} is the number of standard deviations and \eqn{\xi} is calculated as
	\deqn{\xi=\frac{1}{\sqrt{\chi_{p-1}^{2}\left(\nu\eta^{2}+\left(1-p_{\eta}\right)\eta^{2}\right)}}}{
		xi = 1/sqrt(pchisq(degfree*eta^2, degfree + 2) + (1-prob.eta)*eta^2)}

If the \eqn{s_i}{s[i]} are ranges of two values, ISO 5725 recommends 
carrying out the above iteration on the ranges and then dividing by 
\eqn{\sqrt{\nu+1}}{sqrt{degfree+1}}; in the implementation here, this
is done prior to returning the estimate.

If \code{verbose} is non-zero, the current iteration number
and estimate are printed; if \code{verbose>1}, the current set
of truncated values \eqn{w}{w} is also printed.
}
\value{
A scalar estimate of poooled standard deviation.
}
\references{
ISO 5725-5:1998 Accuracy (trueness and precision) of measurement methods and 
results - Part 5: Alternative methods for the determination of the precision 
of a standard measurement method
}
\author{
S L R Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link{algA}}
}
\examples{
#example from ISO 5725-5:1998 (cell ranges for percent creosote)

cdiff <- c(0.28, 0.49, 0.40, 0.00, 0.35, 1.98, 0.80, 0.32, 0.95)

algS(cdiff, is.range=TRUE)
	

#Compare with the sd of the two values (based on the range)
c.sd <- cdiff/sqrt(2)
algS(c.sd, degfree=1, verbose=TRUE)


}
\keyword{ robust }
