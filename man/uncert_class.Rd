\name{uncert-class}
\Rdversion{1.1}
\alias{print.uncert}
\alias{summary.uncert}
\alias{uncert-class}
\title{
The 'uncert' class
}
\description{
Object returned by \code{uncert} calls.
}
\usage{
   \method{print}{uncert}(x, digits=NULL, right=FALSE, \dots, simplify=TRUE)

   \method{summary}{uncert}(object, \dots, simplify=TRUE)

}
\arguments{
  \item{x, object}{An object of class uncert}
  \item{digits}{Number of digits to display in budget and (if present) distribution parameter lists.
     Passed to \code{format} for distribution parameter list and to \code{print.data.frame} 
     for output.}
  \item{right}{If TRUE, strings in uncertainty budget are right-justified. 
  	This differs from the default in \code{print.data.frame}.}
  \item{\dots}{Other parameters passed to \code{print.data.frame}}
  \item{simplify}{If \code{TRUE}, only the call, evaluation method, budget, value \code{y} 
    and combined uncertainty (u.y) are printed.}
}
\details{
\code{summary.uncert} simply calls \code{print.uncert}.

An object of class "uncert" contains:
\describe{
  \item{call}{The matched call}
  \item{y}{The calculated value (for function, expression or formula methods) or NA}
  \item{u.y}{The combined standard uncertainty}
  \item{method}{The uncertainty evaluation method used.}
  \item{budget}{A data frame consisting of:
     \describe{
     \item{x}{\code{x} (if supplied; otherwise a vector of NA's).}
     \item{u}{The standard uncertainties in input quantities (originally provided as \code{u})}
     \item{df}{The degrees of freedom asscociated with \code{u}}
     \item{c}{Sensitivity coefficients either provided as \code{c} or (for the formula, function 
        and expression methods) as calculated.}
     \item{u.c}{The product of \code{u} and \code{c}. These are the contributions to the 
        combined uncertainty for uncorrelated quantities.}
     }
  }
  \item{additional}{Any relevant parameters other than those in $budget$x (typically addditional 
    constants passed to function or expression methods)}
  \item{distrib}{If available, a named list of the distributions associated with \code{u}. 
     The list contains either root nams of distribution functions (e.g \code{"norm"} or
     function definitions.}
  \item{distrib.pars}{If available, a list of lists of parameters describing the 
     distributions associated with \code{u}.}
  \item{cov}{The covariance matrix used}
  \item{cor}{The correlation matrix used}
  \item{cov.xy}{A data frame of covariances between x and y. Row names correspond
  	to the correlation method used. For all uncertainty evaluation methods but MC,
  	the only correlation calculation is "theoretical"; for MC row names include
  	all methods supported by stats::cor at the time the object was created.}
  \item{cor.xy}{A data frame of correlations between x and y, of the same form as \code{cov.xy}}
  \item{deriv}{For the formula and expression methods, the result of a call to \code{deriv}; 
     an expression which evaluates to the value with attributes corresponding to the derivatives 
     (that is, an expression which can be evaluated to give the value and sensitivity coefficients)}
}     
}
\section{Methods}{
  \describe{
    \item{print}{The print method provides a formatted printout of the object. By default,
      \code{simplify=TRUE}; this displays a shortened listing. Columns in \code{$budget} are 
      	suppressed if all NA (typically df when not specified).}
    \item{summary}{summary is currently an alias for the print method.}
  }
}
\value{
	\code{print} and \code{summary} methods invisibly return the original object.
}
\references{
None.
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}
\seealso{
  \code{\link{uncert}}, especially for calculation methods; 
  \code{\link{plot.uncert}}, \code{\link{uncertMC-class}},
  \code{\link[base]{print.data.frame}}, \code{\link[base]{format}}.
}
\examples{

  expr <- expression(a+b*2+c*3+d/2)
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.expr<-uncert(expr, x, u, method="NUM")
  print(u.expr)

}

\keyword{ classes }

