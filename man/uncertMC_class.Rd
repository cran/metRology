\name{uncertMC-class}
\Rdversion{1.1}
\alias{uncertMC-class}
\alias{print.uncertMC}
\alias{summary.uncertMC}
\title{
The 'uncertMC' class
}
\description{
Object returned by \code{uncertMC} calls and by \code{uncertainty} with \code{method="MC"}.

\code{summary.uncertMC} is currently an alias for \code{print.uncertMC}.
}
\usage{
   \method{print}{uncertMC}(x, digits=NULL, right=FALSE, 
                      \dots, simplify=TRUE, minimise=FALSE)

   \method{summary}{uncertMC}(object, digits=NULL, right=FALSE, 
                      \dots, simplify=TRUE, minimise=FALSE)

}
\arguments{
  \item{x, object}{An object of class \code{"uncertMC"}}
  \item{digits}{Number of digits to display in budget and (if present) distribution parameter lists.
     Passed to \code{format} for distribution parameter list and to \code{print.data.frame} 
     for output.}
  \item{right}{If TRUE, strings in uncertainty budget are right-justified. 
  	This differs from the default in \code{print.data.frame}.}
  \item{\dots}{Other parameters passed to \code{print.data.frame}}
  \item{simplify}{If \code{TRUE}, only the call, evaluation method, budget, value \code{y} 
    and combined uncertainty (u.y) are printed.}
  \item{minimise}{If \code{TRUE}, the header, call, \code{expr} and evaluation method 
    are suppressed; this is the mode used when printing an \code{uncertMC} object as part of an 
    \code{uncert} object.} 
}
\details{
An object of class "uncertMC" inherits from class "uncert". In addition to the 
contents of the "uncert" object, it contains the results from the MC replication as 
a list MC. The complete description is:
\describe{
  \item{call}{The matched call}
  \item{expr}{The expression, formula or function supplied to \code{uncertMC}.}
  \item{method}{The uncertainty evaluation method used (always 'MC').}
  \item{B}{The number of Monte Carlo replicates used.}
  \item{budget}{A data frame consisting of:
     \describe{
     \item{x}{The starting values \code{x}.}
     \item{u}{The standard uncertainties in input quantities (originally provided as \code{u})}
     \item{df}{The degrees of freedom asscociated with \code{u}}
     \item{c}{Sensitivity coefficients estimated from the MC output (see \code{\link{uncertMC}} 
        for how this is done).}
     }
  }
  \item{distrib}{If available, a named list of the distributions associated with \code{u}. 
     The list contains either root nams of distribution functions (e.g \code{"norm"} or
     function definitions.}
  \item{distrib.pars}{If available, a list of lists of parameters describing the 
     distributions associated with \code{u}.}
  \item{additional}{If supplied, any relevant parameters other than those in \code{$budget$x} 
    (typically addditional constants passed to function or expression methods)}
  \item{cov}{The covariance matrix used}
  \item{cor}{The correlation matrix used}
  \item{cov.xy}{A data frame of covariances between x and y. The Row names correspond
  	to the correlation method used. For \code{uncertMC} objects only
  	\code{"pearson"} is currently supported (because \code{"kendall"} and \code{"spearman"}
  	take a very long time to compute)}
  \item{cor.xy}{A data frame of correlations between x and y, of the same form as \code{cov.xy}}
  \item{MC}{
     A list containing:
     \describe{
	  \item{seed}{The value of \code{.Random.seed} when \code{uncertMC} was called.}
	  \item{y}{The \code{B} Monte Carlo replicates of the standard uncertainty 
	     calculated as \code{sd(y)}.}
	  \item{x.mc}{If \code{uncertMC} is called with \code{keep.x=TRUE}, a data frame 
	     whose columns are the Monte Carlo replicates of the variables in \code{x}.}
	  }
  }
}     
}
\section{Methods}{
  \describe{
    \item{print}{The print method provides a formatted printout of the object. By default,
      \code{simplify=TRUE}; this displays a shortened listing. Columns in \code{$data} are 
      	suppressed if all NA.}
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
  \code{\link{uncert}}, \code{\link{uncertMC}}, \code{\link{uncert-class}}, 
  \code{\link[base]{print.data.frame}}, \code{\link[base]{format}}.
}
\examples{

  set.seed(13*17)
  expr <- expression(a+b*2+c*3+d/2)
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.expr<-uncertMC(expr, x, u, distrib=rep("norm", 4), method="MC")
  print(u.expr)

}

\keyword{ classes }

