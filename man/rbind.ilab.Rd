\name{rbind.ilab}
\Rdversion{1.1}
\alias{rbind}
\alias{cbind}
\alias{rbind.default}
\alias{cbind.default}
\alias{rbind.ilab}
\alias{cbind.ilab}
\alias{c.ilab}
\title{
Combine 'ilab' objects
}
\description{
Functions to combine ilab objects.
}
\usage{
rbind(..., deparse.level = 1)
\method{rbind}{default}(..., deparse.level = 1)
\method{rbind}{ilab}(..., deparse.level = 1)

\method{c}{ilab}(..., recursive=FALSE)

cbind(..., deparse.level = 1)
\method{cbind}{default}(..., deparse.level = 1)
\method{cbind}{ilab}(..., deparse.level = 1)
}

\arguments{
  \item{\dots}{For \code{rbind} and \code{c}, objects of class \dQuote{ilab}
     to be combined. For \code{cbind}, one \dQuote{ilab} object and vectors, scalars 
     or data frames to be appended to the \dQuote{ilab} object's $data element. See Details.
  }
  \item{deparse.level}{integer controlling the construction of labels. Passed 
  to \code{\link[base:cbind]{rbind.data.frame}} or \code{base::cbind}.
  }
  \item{recursive}{logical, controlling recursion in lists. Included solely for 
     consistency with \code{base::c}; no effect in \code{c.ilab}. 
  }
}
\details{
The generic and default cbind and rbind functions defined by metRology use the first 
object in \code{'...'} to decide which method to apply. This differs from the behaviour 
of these functions in the base package, which dispatch based on inspection of 
all objects in \code{'...'} (see \code{\link[base:cbind]{rbind}} in the base package for details). 
Control is, however, passed to the base package functions if an ilab object is not first 
in the list. See section \dQuote{Warning} below for a work-round if this causes difficulty.

The \code{rbind} method for class \dQuote{ilab} combines objects by applying
\code{rbind} to the \code{$data} elements in turn and then concatenating the 
\code{$distrib} and \code{$distrib.pars} elements using the default \code{c} 
method. All objects combined by \code{rbind.ilab} must be of class \dQuote{ilab}.

Combination of the \code{$data} elements follows the rules of \code{rbind.data.frame}; 
in particular, names must match, but need not be in the same order and the return value 
column classes will be coerced to match the first if necessary. 

\code{c.ilab} simply passes the objects to \code{rbind.ilab}, using the default value 
for \code{deparse.level}. \code{recursive} is ignored. An error is returned if any 
objects in \code{'...'} are not of class \dQuote{ilab}.

The \code{cbind} method for \dQuote{ilab} objects combines objects of class \dQuote{ilab} 
with atomic objects or data frames by applying \code{base::cbind} to \code{$data} in the 
(single) supplied ilab object and the items listed in \code{'...'}. \code{base::cbind} 
will extend scalars, vectors or columns in data frames to length \code{nrow(ilab$data)} if
their length is an integer fraction of \code{nrow(ilab$data)}. Unlike \code{base::cbind}, 
\code{cbind.ilab} does not permit vectors or data frames longer than \code{nrow(ilab$data)}
and will return an error in such cases. cbind.ilab will also return an error if any objects in 
\code{'...'} are not one of atomic, data frame or class \dQuote{ilab}, if more than one 
\sQuote{ilab} object is supplied or if none are. 
}

\value{
An object of class \sQuote{ilab}. The title for the returned object is the title for 
the first ilab object in the list.

For the \code{cbind} method, the returned object will have additional columns in the 
\code{$data} element, and the title will be unchanged.
}

\section{Warning}{
Because of the unusual method dispatch behaviour of \code{base::cbind} and 
\code{base::rbind}, which use the data frame method if \emph{any} objects in
\code{'...'} are data frames, metRology masks the base package functions 
in order to guarantee correct dispatch when data frames are included in \code{'...'}. 
No adverse effects are currently known, but please report any to the package maintainer.
Calling \code{base::cbind} directly will provide a work-around if necessary.
}

\references{
None, yet.}
\author{
S. L. R. Ellison \email{s.ellison@lgcgroup.com}
}


\seealso{
base package functions \code{\link[base]{cbind}}, \code{\link[base]{c}}.
  
}
\examples{
data(Pb)
il1<-construct.ilab(org=Pb$lab, x=Pb$value, measurand="Pb", item="none", 
                u=Pb$u, k=Pb$k, U=Pb$U, title=c("CCQM K30", "Lead in wine"),
                method=Pb$method)

rbind(il1, il1)

}

\keyword{manip}
