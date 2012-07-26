\name{ilab-class}
\Rdversion{1.1}
\alias{ilab-class}
\alias{construct.ilab}
\title{
The 'ilab' class. 
}
\description{
The \sQuote{ilab} class and its constructor function.
}
\usage{
construct.ilab(org, item, measurand, x, u, df, k, U, U.lower, U.upper,
           distrib=NULL, distrib.pars=NULL, study=NA, title=NA, p=0.95, ...) 

}

\arguments{
  \item{org}{Character vector or factor of organisation names.}
  \item{item}{vector or factor of identifiers for test items. Coerced to factor on storage.}
  \item{measurand}{vector or factor identifying the measurand(s) involved in the study.}
  \item{x}{numeric vector of reported values.}
  \item{u}{numeric vector of reported standard uncertainties or standard errors associated with x.}
  \item{df}{optional numeric vector of degrees of freedom associated with each reported uncertainty.}
  \item{k}{numeric vector of coverage factors. The coverage factor is the factor multiplying u to obtain U.}
  \item{U}{numeric or character vector of expanded uncertainties or confidence interval half-widths. 
      Coerced to numeric but may include a character representation of interval limits; see Details.}
  \item{U.lower, U.upper}{numeric vectors of lower and upper limits for the confidence 
      interval around x, allowing asymmetric intervals. Defaults to U or to the limits 
      specified by U. See Details.}
  \item{distrib}{A character vector of length \code{length(x)} or a named list 
        of names of distribution functions associated with \code{u}. If a character vector, 
        \code{distrib} is recycled to length \code{length(x)}.
        }
  \item{distrib.pars}{A named list of lists of parameters describing the distributions 
  	associated with \code{u} to be passed to the relevant distribution function.
  	If \code{distrib} is present but \code{distrib.pars} is not, an attempt is made 
  	to set defaults based on other parameters; see Details.}
  \item{study}{A character value or vector or a factor identifying different studies 
      or study populations within the data set. Typically used, for example, for identifying 
      participants in global and regional components of a combined study. Recycled to length
      \code{length(x)} if necessary.}
  \item{title}{An optional title for the study. May be a character vector, 
  	in which case each element is displayed on a separate line when printed.}
  \item{p}{Confidence level assumed to apply to \code{k}. Used only to set a default value 
        for \code{df} when \code{distrib} indicates a t-distribution and \code{df} is 
        unspecified.}
  \item{\dots}{Other \emph{named} factors or character vectors used to group observations.}

}

\details{
If \code{U} is a character vector, it may contain character representations of range.
Two forms are permitted:
  \describe{
    \item{"a-b"}{Interpreted as limits of a range from \code{a} to \code{b}. \code{U.lower}
       and \code{U.upper} are calculated from these limits and \code{x}}
    \item{"+a[/]-b" (or "-a[/]+b")}{\code{U.upper} is set to \code{a} in \code{"+a"}, 
       and \code{U.lower} is set to \code{b} in \code{"-b"}.}
  }

If \code{distrib.pars} is missing, an attempt is made to deduce appropriate 
distribution parameters from \code{x}, \code{u}, \code{df} and \code{distrib}. 
In doing so, the following assumptions and values apply for the respective distributions:
\describe{
   \item{norm}{\code{mean=x$name, sd=u$name}.}
   \item{unif}{\code{min=x-sqrt(3)*u, max=x+sqrt(3)*u}.}
   \item{tri}{\code{min=x-sqrt(6)*u, max=x+sqrt(6)*u, mode=x}.}
   \item{t, t.scaled}{\code{df=df, mean=x, sd=u}.}
}
In addition, if \code{distrib} contains \code{"t"} or \code{"t.scaled"}, and
\code{df} is \code{NA}, the corresponding degrees of freedom are chosen based on 
\code{k} and \code{p}.
}

\value{
An object of class \sQuote{ilab} consisting of:
\item{title}{A character value or vector describing the study}
\item{subset}{A character string describing any subset operation used to form the object.}
\item{data}{A data frame with columns:
   \tabular{ll}{
    org	\tab	Factor of organisations submitting results in the study  \cr
    item	\tab	Factor of test item identifiers.   \cr
    measurand	\tab	Factor of measurands determined for each item  \cr
    x	\tab	numeric vector of reported values.  \cr
    u	\tab	numeric vector of reported standard uncertainties or standard 
    		errors associated with x.  \cr
    df	\tab	numeric vector of degrees of freedom associated with each reported
    		uncertainty. Set to \code{NA} if not provided.  \cr
    k	\tab	numeric vector of coverage factors. The coverage factor is the 
    		factor multiplying u to obtain U.  \cr
    U	\tab	numeric or character vector of expanded uncertainties or confidence 
    		interval half-widths. \code{U} is coerced to numeric but may include
    		a character representation of interval limits; see Details.  \cr
    U.lower, U.upper	\tab	numeric vectors of lower and upper limits for the confidence 
      interval around x.  \cr
    study	\tab	Identifier for study groups (see Arguments above).  \cr
    \dots	\tab	Other grouping factors (supplied in \sQuote{\ldots} in \code{construct.ilab})
	which can be used for sub-categorisation.  \cr
  }
% \describe{
%    \item{org}{Factor of organisations submitting results in the study} 
%    \item{item}{Factor of test item identifiers. } 
%    \item{measurand}{Factor of measurands determined for each item} 
%    \item{x}{numeric vector of reported values.}
%    \item{u}{numeric vector of reported standard uncertainties or standard errors associated with x.}
%    \item{df}{numeric vector of degrees of freedom associated with each reported uncertainty. Set to 
%      \code{NA} if not provided.}
%    \item{k}{numeric vector of coverage factors. The coverage factor is the factor multiplying u to obtain U.}
%    \item{U}{numeric or character vector of expanded uncertainties or confidence interval half-widths. 
%      Coerced to numeric but may include a character representation of interval limits; see Details.}
%    \item{U.lower, U.upper}{numeric vectors of lower and upper limits for the confidence 
%      interval around x.}
%    \item{study}{Identifier for study groups (see Arguments above).}
%    \item{\dots}{Other grouping factors (supplied in \sQuote{\ldots} in \code{construct.ilab})
%	which can be used for sub-categorisation.}
%  }
}
\item{distrib}{An unnamed list of distribution names.}
\item{distrib.pars}{An unnamed list of lists of parameters describing the distributions 
  	associated with \code{u}.}
}
\references{
None, yet.}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link{print.ilab}},  \code{\link{subset.ilab}},  \code{\link{plot.ilab}}
}
\examples{
data(Pb)
construct.ilab(org=Pb$lab, x=Pb$value, measurand="Pb", item="none", 
                u=Pb$u, k=Pb$k, U=Pb$U, title=c("CCQM K30", "Lead in wine"), 
                method=Pb$method)

#Illustrate default for U and automatic distrib.pars
construct.ilab(org=Pb$lab, x=Pb$value, measurand="Pb", item="none", 
                u=Pb$u, k=Pb$k, distrib="norm")

construct.ilab(org=Pb$lab, x=Pb$value, measurand="Pb", item="none", 
                u=Pb$u, k=Pb$k, distrib="t.scaled")

}

\keyword{ methods }
