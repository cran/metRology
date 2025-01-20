\name{LCS}
\Rdversion{1.1}
\alias{LCS}

\title{LCS: Largest consistent subset}
\description{
Calculates a 'largest consistent subset' given values and associated uncertainty information.
}
\usage{
LCS(x, u, p = 0.05, method = "enum", simplify = FALSE, 
	verbose = FALSE)
}
\arguments{
  \item{x}{Vector of observations.}
  \item{u}{Vector of standard errors or standard uncertainties associated with \code{x}.}
  \item{p}{Significance level at which consistency is tested.}
  \item{method}{Subset identification method. Currently only 'enum' is supported.}
  \item{simplify}{If \code{simplify} is \code{TRUE}, only the lowest-uncertainty 
  	subset is returned even if several are of the same size.}
  \item{verbose}{Logical: Controls the level of reporting during the search.}
}
\details{
\code{LCS} obtains the largest subset(s) of \code{x} which pass a chi-squared 
test for consistency, taking the uncertainties \code{u} into account. 

\code{method} controls the search method used. Method "enum" uses complete enumeration 
of all subsets of size \code{n}, starting at \code{n==length(x)} and decreasing \code{n} 
until at least one consistent subset is found. No other method is currently supported; if 
a different method is specified, LCS provides a warning and continues with "enum".

There may be more than on consistent subset of size n. If so, LCS returns all such 
subsets unless \code{simplify} is \code{TRUE}, in which case LCS prints a short warning 
and returns the subset with smallest estimated uncertainty as estimated for the Graybill-Deal
weighted mean assuming large degrees of freedom in \code{u}. 

\code{verbose} controls the level of reporting. If \code{TRUE}, LCS prints the progress of 
the search.

The general idea of a Largest Consistent Subset as implemented here was suggested by Cox (2006),
though at least one other related method has been suggested by Heydorn (2006). It has,
however, been criticised as an estimator (Toman and Possolo (2009)) ; see Warning below.
}
\value{
If there is only one subset of maximum size, or if \code{simplify=TRUE}, a vector of indices 
for \code{x} representing the largest consistent subset.

If there is more than one subset of maximum size and \code{simplify=FALSE}, a matrix of indices 
in which the rows contain the indices of each subset.

}
\references{
Cox, M. G. (2007) The evaluation of key comparison data: determining the 
largest consistent subset. Metrologia \bold{44}, 187-200 (2007) 

Heydorn, K. (2006) The determination of an accepted reference value from 
proficiency data with stated uncertainties. Accred Qual Assur \bold{10}, 479-484 (2006)

Toman, B. and Possolo, A. (2009) Laboratory effects models for interlaboratory
comparisons. Accred. Qual. Assur. \bold{14}, 553-563 (2009)
}
\author{
S. Ellison \email{s.ellison@lgcgroup.com}
}

\section{Warning }{
LCS methods are essentially equivalent to unsupervised outlier rejection. In general,
this results in a possibly extreme low estimated variance for an arbitrarily small subset 
(in the limit of gross inconsistency, LCS will return subsets of size 1). The estimated uncertainty 
calculated for the Graybill-Deal weighted mean of the subset(s) does not generally take 
account of the subset selection process or the dispersion of the complete data set, so 
is \emph{not} an estimate of sampling variance. 

LCS is therefore not recommended for consensus value estimation. It is however, quite 
useful for identifying value/uncertainty outliers.

} 

\seealso{
None.
}
\examples{
data(Pb)
with(Pb, LCS(value, U/k))

}
\keyword{ misc }

