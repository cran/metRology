\name{update.uncert}
\Rdversion{1.1}
\alias{update.uncert}
\title{
Update and recalculate an uncertainty estimate
}
\description{
\sQuote{update} allows modification of components of \sQuote{uncert} or \sQuote{uncertMC} objects, 
including the uncertainty estimation method used, and will recalculate the estimate and return a new
 \sQuote{uncert} or \sQuote{uncertMC} object.  Individual elements of most components can be amended.
}
\usage{
\method{update}{uncert}(object, expr = NULL, method = NULL, x = NULL, u = NULL, c=NULL, 
		df = NULL, cov = NULL, cor = NULL, distrib = NULL, 
		distrib.pars = NULL, delta = NULL, B = NULL, keep.x = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \sQuote{uncert}}
  \item{expr}{An expression, formula or function.}
  \item{method}{
Uncertainty evaluation method. May be any of the methods listed for \code{\link{uncert}}. 
}
  \item{x, u, df}{Named list, vector or array of values to update elements of 
      object \code{uncert}. See Details for options.
}
  \item{c}{Update to \code{uncert$budget$c} No effect except for updates using \code{uncert.default}.
}
  \item{cov, cor}{A covariance or correlation matrix. Only one of \code{u} and \code{cov} 
  should be specified; if both are specified, \code{cov} will take precedence and a warning will be issued.}
  \item{distrib}{Named list or character vector of updated distribution names.}
  \item{distrib.pars}{Named list of updates for distribution names and parameters.}
  \item{delta}{Scalar value updating \code{delta} for numeric evaluation methods.}
  \item{B}{Updated number of Monte Carlo iterations for \sQuote{uncertMC} objects or specification of
  \code{B} when updating \sQuote{uncertMC} objects using \code{method="MC"}.
  }
  \item{keep.x}{Update to \code{keep.x} passed to \code{uncertMC} for Monte Carlo updates.}
  \item{\dots}{Other values passed to \code{uncert} or \code{uncertMC}}
}
\details{
Update will use the values provided to update the object given, call the original function 
with the revised parameters and return the result as an object of class \sQuote{uncert} 
or \sQuote{uncert} depending on the uncertainty evaluation method used. 

Note that updating with a different value of \code{method} may result in an object of 
different class. Updating an \sQuote{uncertMC} object with a method other than \code{"MC"} 
will return an object of class \sQuote{uncert}; similarly, updating an \sQuote{uncert} object 
using \code{method="MC"} will return an object of class \sQuote{uncertMC}.

Updates to vector or list elements of \code{uncert} such as \code{x}, \code{u}, \code{df} 
etc. can be specified as named lists, named vectors or arrays, with names corresponding to 
names of the input quantities in the uncertainty budget (that is, the names may correspond
to one or more of \code{row.names(uncert$budget)}). If names are present, only the 
corresponding individual members are updated.
If names are not present, the complete vector or list in \code{uncert} is replaced, and names added.

Matrix elements \code{cor} and \code{cov} must be specified completely; see
\code{\link{buildCor}}, \code{\link{updateCor}} and associated functions for
covariance matrices for compact update methods.
}
\value{
An object of class \sQuote{uncert} or, for \code{method="MC"} of class \sQuote{uncertMC}. 
See \code{\link{uncert-class}} and \code{\link{uncertMC-class}} for details.

}
\references{
None, yet.
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link{uncert-class}},  \code{\link{uncert-class}}, \code{\link{uncertMC}}, \code{\link{uncertMC-class}}
}

\examples{

  #From uncert:
  expr <- expression(a+b*2+c*3+d/2)
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.expr<-uncert(expr, x, u, method="NUM")
  u.expr

  update(u.expr, u=list(a=0.3))
  
  update(u.expr, method="MC")
  
}

\keyword{ univar }
