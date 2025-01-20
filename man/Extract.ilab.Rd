\name{Extract.ilab}
\Rdversion{1.1}
\alias{Extract.ilab}
\alias{subset.ilab}
\alias{[.ilab}
\title{
The 'ilab' class. 
}
\description{
Functions for manipulating interlaboratory study objects objects of class \sQuote{ilab}.
}
\usage{
\method{subset}{ilab}(x, subset, drop=FALSE, \dots)

\method{[}{ilab}(x, i, j)
}

\arguments{
  \item{x}{An object of class \sQuote{ilab}}
  \item{subset}{logical expression indicating elements or rows to keep:
          missing values are taken as false.}
  \item{drop}{ passed on to '[' indexing operator.}
  \item{\dots}{Parameters passed to other functions}

  \item{i, j}{elements to extract. May be numeric or logical vectors.  }

}

\details{

For the subset method, \code{subset} is an expression evaluated in the frame
of \code{ilab$data} and in the parent environment if objects are not found in 
\code{ilab$data}. Note that since \code{ilab$distrib} and \code{ilab$distrib.pars} are
not in \code{ilab$data}, any operation on these must be specified in full.

The indexing method \code{'['} operates on both rows and columns of the object. However,
only the \code{$data} element can be addressed with the \code{j}; the \code{distrib} and 
\code{distrib.pars} elements are unaffected by \code{j} and will always be included 
in the returned object.
}

\value{
An object of class \sQuote{ilab} with fewer rows and (if \code{j} is present) 
fewer columns. 
}

\section{Warning}{
Removing the standard columns from \sQuote{ilab} objects using \code{'['} may have unforeseen 
consequences for other functions; only the print method is likely to operate successfully. 
}

\references{
None, yet.}
\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}

\seealso{
\code{\link{ilab-class}}. 
}
\examples{
data(Pb)
il.pb<-construct.ilab(org=Pb$lab, x=Pb$value, measurand="Pb", item="none", 
                u=Pb$u, k=Pb$k, U=Pb$U, title=c("CCQM K30", "Lead in wine"), method=Pb$method)

subset(il.pb, u < 0.03)

il.pb[1:6,]

}

\keyword{manip}
