\name{methods.ilab}
\Rdversion{1.1}
\alias{methods.ilab}
\alias{print.ilab}
\alias{plot.ilab}
\title{
Methods for the 'ilab' class. 
}
\description{
Functions for printing and plotting interlaboratory study objects objects of class \sQuote{ilab}.
}
\usage{
\method{print}{ilab}(x, \dots, digits = NULL, right = FALSE)

\method{plot}{ilab}(x, \dots)
}

\arguments{
  \item{x}{An object of class \sQuote{ilab}}
  \item{digits}{Number of digits to display in budget and (if present) distribution parameter lists.
     Passed to \code{format} for distribution parameter list and to \code{print.data.frame} 
     for output.}
  \item{right}{If TRUE, strings in uncertainty budget are right-justified. 
  	This differs from the default in \code{print.data.frame}.}
  \item{\dots}{Parameters passed to other functions}

}

\details{
The print method uses \code{print.data.frame} to display the data after formatting the 
\code{distrib} and \code{distrib.pars} elements. 

The plot method passes the object to \code{kplot}.

}

\value{
The print and plot methods are called for their side effects.

}
\references{
None, yet.}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{link{ilab-class}}, \code{link{subset.ilab}} \code{link{kplot}}.
}
\examples{
data(Pb)
il.pb<-construct.ilab(org=Pb$lab, x=Pb$value, measurand="Pb", item="none", 
                u=Pb$u, k=Pb$k, U=Pb$U, title=c("CCQM K30", "Lead in wine"), method=Pb$method)

print(il.pb)

plot(il.pb)

}

\keyword{ methods }
