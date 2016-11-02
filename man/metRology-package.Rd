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

Improvements in this version (0.9-23-2) from version 0.9-22 include: 
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
