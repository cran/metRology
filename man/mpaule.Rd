\name{Mandel-Paule}
\alias{mpaule}
\alias{mpaule.default}
\alias{mandel.paule}
\title{Mandel-Paule calculation}
\description{
  Calculate a weighted mean, between-group standard deviation and standard error on the
  weighted mean using the Mandel-Paule algorithm. 
}
\usage{

mpaule(x, \dots, tol=.Machine$double.eps^0.25, maxiter=25) 

\method{mpaule}{default}(x, u=NULL, n=NULL, groups=NULL, 
                   tol=.Machine$double.eps^0.25, maxiter=25, \dots)

mandel.paule(x, \dots, tol=.Machine$double.eps^0.25, maxiter=25) 
}
\arguments{
  \item{x}{numeric vector of mean values for groups, or (if \code{groups} 
  	is given) of individual observations}
  \item{u}{numeric vector of standard deviations or standard uncertainties 
  	associated with the values \code{x}}
  \item{n}{integer vector of numbers in each group. If \code{NULL}, \code{u} are interpreted as 
      standard uncertainties or standard errors. \code{n} is recycled to \code{length(x)}.}
  \item{groups}{factor, or vetor which can be coerced to factor, of groups. If present, \code{x} is 
    interpreted as a vector of individual observations and \code{u} and \code{n} ignored.}
  \item{\dots}{Additional parameters passed to other functions.}
  \item{tol}{numeric tolerance; iteration stops when the variance step size drops below \code{tol*var(x)}}
  \item{maxiter}{numeric maximum number of iterations}
}
\details{
  The Mandel-Paule algorithm finds the between-method variance by iteratively solving the
  equation relating the weighted mean to the weighting factor applied. The weighting factor is
  the inverse of the sum of the standard error in \code{x} and the between-group variance.
  
  If the iterative procedure produces a negative estimate for the between-group variance,
  the between-group variance is set to zero.

  For the default method, if \code{u} is present and \code{n=NULL}, \code{u} is interpreted as 
  a vector of standard uncertainties or standard errors. If \code{n} is not \code{NULL}, \code{u} 
  is interpreted as a vector of standard deviations and standard errors are calculated as 
  \code{u/sqrt(n)}.
  If \code{groups} is not \code{NULL}, \code{x} is interpreted as a vector of individual 
  observations grouped by \code{groups}, and the algorithm is applied to the corresponding 
  group means and standard errors. 
  
  If \code{maxiter} is set less than 1, no iterations are performed and the consensus mean 
  is returned as \code{NA}.
  
  \code{mandel.paule} is an alias for \code{mpaule} retained for backward compatibility.
}
\value{
  A loc.est object; see loc.est for details. In the returned object, \code{df} is
  set to \eqn{n-1} where \eqn{n} is the number of non-\code{NA} observations or
  group means as appropriate, and \code{method.details} is returned as :
  \item{var.between}{the estimated between-group variance)}
  \item{iter}{the number of iterations taken}
  \item{converged}{\code{converged} indicates the convergence status. \code{0L} indicates
  failure to converge (\code{maxiter} reached before step size drops below tolerance); 
  \code{1L} indicates normal convergence; \code{2L} indicates that the final step size resulted 
  in a negative between-group variance, at which point the variance and step size are set to 0.0}
}
\references{
Paule, R. C.  and Mandel, J. (1982), J Res Nat Bur Stand, \bold{87}, (5) 377-385 
}

\author{
S. Cowen \email{simon.cowen@lgc.co.uk} with amendments by S. L. R. Ellison.
}
\examples{
## the second example in the paper cited above
x <- c(201.533, 216.55)
s <- c(0.154, 0.25)
n <- c(6, 2)

mpaule(x, s/sqrt(n))

## Cd heat of vapourisation example from the paper cited above
x2<-c(27.044, 26.022, 26.340, 26.787, 26.796)
v<-c(3, 76, 464, 3, 14)*1e-3
mpaule(x2, sqrt(v))

}
\keyword{univar}
