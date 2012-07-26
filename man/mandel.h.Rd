\name{mandel.h}
\alias{mandel.h}
\alias{mandel.h.default}
\alias{mandel.h.data.frame}
\alias{mandel.h.matrix}
\alias{mandel.h.array}
\alias{mandel.h.ilab}
\title{
Calculate Mandel's h statistics for replicate observations
}
\description{
\code{mandel.h} calculates Mandel's h statistics for replicate observations. 
Mandel's h is an indication of relative deviation from the mean value.
}
\usage{
	mandel.h(x, g = NULL, m = NULL, na.rm = T, rowname = NULL, 
			method=c("classical", "robust"), n = NA, ...)

	\method{mandel.h}{default}(x, g = NULL, m = NULL, na.rm = T, rowname = NULL, 
			method=c("classical", "robust"), n = NA, ...)

	\method{mandel.h}{ilab}(x, g = NULL, m = NULL, na.rm = T, rowname = NULL, 
			method=c("classical", "robust"), n = NA, ...)
}
\arguments{
  \item{x}{An R object (see Details below), which contains replicate observations or, 
  		if \code{g} is absent, means or standard deviations.}
  \item{g}{A primary grouping factor, usually corresponding to Laboratory in an
  		inter-laboratory study. If not present, \code{x} is taken as 
  		a set of means or standard deviations (depending on whether 
  		\code{type} is \code{"h"} or \code{"k"}.}
  \item{m}{A secondary grouping factor, usually corresponding to test item 
  		or measured quantity. \code{m} is ignored if \code{x} has
  		more than one column.}
  \item{na.rm}{A logical value indicating whether 'NA' values should be
          stripped before the computation proceeds. Passed to functions 
          such as \code{mean} and \code{sd}.}
  \item{rowname}{A single character label for the primary grouping factor 
  	(e.g. "Lab", "Organisation").}
  \item{method}{Character scalar giving the calculation method. \code{"classical"} gives the
  	traditional calculation; \code{"robust"} gives a robust variant (see Details).}
  \item{n}{scalar number of observations per group. Required only if \code{x} consists of 
  	calculated standard deviations.}
  \item{\dots}{Additional parameters passed to \code{hubers} when \code{method="robust"}.}
}
\details{
\code{mandel.h} is a convenience wrapper for mandel.kh(..., type="h"). It is generic, 
with methods for numeric vectors, arrays, data frames, matrices and objects of class 
\code{'ilab'}. All parameters are passed to \code{mandel.kh}.

Mandel's h is an indicators of relative deviation for grouped 
sets of observations. Given a set of observations \eqn{x_{ijl}}{x[i,j,l]} where \eqn{i, j, l}
denotes observation \eqn{l}, \eqn{l=1, 2, ... n} for measurand or test item \eqn{j} and group
(usually laboratory) \eqn{i}, \eqn{i=1, 2, ... p}, Mandel's \eqn{h} is given by:

\deqn{h=\frac{\bar{x_{ij}}-\bar{x_j}}{s_j}}

where
\eqn{ s_j=\sqrt{\sum_{i=1}^p{\frac{(\bar{x_{ij}}-\bar{x_j})}{p-1}}}}

If \code{x} is a vector, one-dimensional array or single-column matrix, values are aggregated 
by \code{g} and, if present, by \code{m}. If \code{x} is a data frame or matrix, each column 
is aggregated by \code{g} and \code{m} silently ignored if present. In all cases, if \code{g}
is \code{NULL} or missing, each row (or value, if a vector) in \code{x} 
is taken as a pre-calculated mean (for Mandel's \eqn{h}) or standard deviation (for Mandel's \eqn{k}).

If \code{x} is an object of class \code{'ilab'}, \code{g} defaults to \code{'$org'} and 
\code{m} to \code{$measurand}. 

The returned object includes a label (\code{'grouped.by'}) for the primary grouping factor. 
For the \code{'ilab'} method, this is \code{"Organisation"}. For other methods, If \code{rowname} is 
non-null, \code{rowname} is used. If \code{rowname} is NULL, the default is \code{deparse(substitute(g))};
if \code{g} is also NULL or missing, "Row" is used.

If \code{method="robust"}, Mandel's \eqn{h} is replaced by a robust z score calculated by 
replacing \eqn{\bar{x_j}} and \eqn{s_j} with the robust estimates of location and scale 
obtained using Huber's estimate with tuning constant \code{k} set to 1.5 (unless otherwise 
specified in \code{...}). 
}

\value{
mandel.h returns an object of class \code{"mandel.kh"}, which is  a data frame consisting
of the required Mandel's statistics and in which each row corresponds to a level of \code{g}
and each column to a level of \code{m} or (if \code{x} was a matrix or data frame) to the 
corresponding column in \code{x}. In addition to the class, the object has attributes:
\describe{
  \item{'mandel.type'}{\code{"h"} or \code{"k"}}
  \item{'grouped.by'}{Character scalar giving the label used for the grouping 
  		factor \code{g}; see Details above for the defaults.}
  \item{'n'}{Number of observations per group (\code{n} if specified}
 }
}
\references{
Accuracy (trueness and precision) of measurement methods and results -- Part 2: 
Basic method for the determination of repeatability and reproducibility of a 
standard measurement method. ISO, Geneva (1994).
}
\author{S Ellison \email{s.ellison@lgc.co.uk}
}
\seealso{
\code{\link{mandel.k}}, \code{\link{mandel.kh}};  
\code{\link{pmandelh}}, \code{\link{pmandelk}}  for probabilities, quantiles etc.;
\code{\link{plot.mandel.kh}}, \code{\link{barplot.mandel.kh}} for plotting methods.
}
\examples{
	data(RMstudy)

	#Data frame examples: note no secondary grouping factor
	h <- with(RMstudy, mandel.h(RMstudy[2:9], g=Lab))
	plot(h, las=2)

	#Vector variant
	RMstk <- stack(RMstudy[,2:9])
	names(RMstk) <- c("x", "meas")
		#names replace 'values' and 'ind'
	RMstk$Lab <- rep(RMstudy$Lab, 8)
	h2 <- with(RMstk, mandel.h(x, g=Lab, m=meas, rowname="Laboratory"))
		#Note use of rowname to override g
	plot(h2, las=2)
	
	#ilab method
	RM.ilab <- with(RMstk, construct.ilab(org=Lab, x=x, measurand=meas, 
		item=factor(rep("CRM", nrow(RMstk))) ) )

	plot(mandel.h(RM.ilab))
	
	#Robust variant
	hrob <- with(RMstudy, mandel.kh(RMstudy[2:9], g=Lab, type="h", method="robust"))
	plot(hrob, las=2)
	
}
\keyword{univar}
