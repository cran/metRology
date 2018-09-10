\name{metRology-package}
\alias{metRology-package}
\alias{metRology}
\docType{package}
\title{
\packageTitle{metRology}
}
\description{
\packageDescription{metRology}
}
\details{
The metRology package includes functions for:
  \itemize{ 
	\item Plotting for Key Comparisons (dot-and-bar, consistency)
	\item Uncertainty evaluation using algebraic or numeric differentiation, with support for correlation
	\item Monte Carlo evaluation of uncertainty (including correlation for normally distributed variables)
	\item Classes and functions for location estimates for metrology comparisons
	\item Mandel's h and k statistics and plots for interlaboratory studies
	\item Support functions for an excel interface
  } 

Changes in version 0.9-28-1 from version 0.9-28-0 include: 
  \itemize{	
	\item Fixed \code{blockplot} x-axis label, which was incorrect.
 
 }
Changes in version 0.9-28-0 from version 0.9-27-2 include: 
  \itemize{	
	\item A new plot, \code{blockplot}, added. A \dQuote{block plot} is a histogram variant identifiying 
	individual data points, which appear as \dQuote{blocks} in the plot. \code{blockplot} provides for grouped 
	data, which generates vertically separated subplots for each group. Fills and label colours can be specified 
	for each data point.
 
 }
Changes in version 0.9-27-2 from version 0.9-26-2 include: 
  \itemize{	
	\item \code{pmsd} and related functions will now use fast interpolation by default, and
	provide exact values for both odd- and even-\eqn{n} data sets up to \eqn{n=199}. 
	
  	\item \code{gplot} (called by \code{plot.mandel.kh}) now has a \code{spacing} parameter
	which allows finer control of vertical line spacing. 
 }

Changes in version 0.9-26-2 from version 0.9-26-1 include: 
  \itemize{
	\item Fix to a bug in \code{\link{reml.loc}} which failed to report the standard uncertainty \code{u} 
	correctly. 
	\item \code{cplot} now respects \code{cex.axis} as a plot parameter. 
	}

Changes in version 0.9-26-1 from version 0.9-26-0 include: 
  \itemize{
	\item Added plot and barplot methods for MSD class. 
	\item Minor correction to code in \code{msd} to prevent over-replication
	      of estimated \code{s} when \code{s} is a function and returns a vector.
	}

Changes in version 0.9-26 from version 0.9-25 include: 
  \itemize{
	\item \code{\link{msd}} now returns an object of class \code{"MSD"} which
	      includes the original data as attributes to permit subsequent bootstrapping.
	\item A new function, \code{\link{bootMSD}} that performs parametric
		bootstrapping for \code{MSD} objects to obtain critical values 
		and \emph{p}-values for the general case where standard uncertainties/standard 
		errors differ appreciably.
	}

Improvements in version 0.9-25 from version 0.9-23 include: 
  \itemize{
	\item \code{\link{plot.d.ellipse}} now takes default \code{xlab} and \code{ylab} from 
	      dimnames in the supplied \code{cov.dellipse}.
	}

Improvements in version 0.9-23 from version 0.9-22 include: 
  \itemize{
	\item A wholly new Youden plot (see \code{\link{yplot}}), with many options for confidence ellipses
	\item A REML location estimate, \code{\link{reml.loc}}, in addition to \code{\link{vr.mle}};
	      \code{reml.loc} can use means and standard uncertainties/standard errors instead of raw data 
	      and when doing so does not require degrees of freedom. 
	\item Incremental improvements in handling for the median scaled difference measure of anomalies.
	      \code{\link{msd}} is faster and less memory-intensive, and \code{\link{pmsd}} now uses a 
	      beta formulation to extend to very high \code{n} (at least 1e6 - if you feel \emph{very} 
	      patient).
	\item Support for \code{log} and \code{log.p} in \code{\link{dt.scaled}}.
	}
	      
Corrections and bugfixes include: 
  \itemize{
	\item amends \code{\link{plot.drop1.uncert}} to give a plot for each measure of change specified in \code{which}
	\item corrects a grep warning appearing in \code{drop1.uncert};
	\item corrects an unnecessary 'missing \code{u}' error message in version 0.9-22's \code{uncert()} 
              when \code{cov} was specified and \code{u} was not.
	}

}
\author{
\packageAuthor{metRology}

Maintainer: \packageMaintainer{metRology}
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
