\name{plot.uncertMC}
\Rdversion{1.1}
\alias{plot.uncertMC}
\title{
Plot method for 'uncertMC' objects
}
\description{
Plots for uncertainty evaluations produced by \code{uncertMC} or calls 
to \code{uncert} with \code{method=MC}.
}
\usage{
\method{plot}{uncertMC}(x, which = 1:2, 
    main=paste("Monte Carlo evaluation -",deparse(substitute(x))), 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    caption = list("Histogram", "Q-Q plot", "Density", 
    "Correlation x-y", "Covariance x-y"), 
    xlab = paste(deparse(substitute(x)), "$y", sep = ""), 
    ..., cex.caption = 1, cex.main = 1.25, lwd.y = 2, col.y = 2, lty.y,
    col.qqline = NULL, lty.qqline = NULL, lwd.qqline = NULL)
}

\arguments{
  \item{x}{An object of class \code{uncertMC} produced by \code{uncertMC()} or  
    \code{uncert()} with \code{method="MC"}.
}
  \item{which}{Integer in 1:5; the particular variant(s) of plot required. 
     A vector is permitted, in which case plots are produced in ascending order 
     of \code{which}.}
  \item{main}{Main title for the plot}
  \item{ask}{logical; if 'TRUE', the user is _ask_ed before each plot, see 'par("ask=")'}
  \item{caption}{A list of captions for \emph{all} the different plots. 
}
  \item{xlab}{x-axis label, currently passed only to the histogram plot.}
  \item{\dots}{Additional parameters passed to other functions. See details
     for which parameters are passed.}
  \item{cex.caption}{Expansion factor for individual plot captions; as \code{cex} in \code{par}.}
  \item{cex.main}{Expansion factor for main title; as \code{cex.main} in \code{par}.
}
  \item{lwd.y, col.y, lty.y}{Line width and colour for the location line in the histogram
    and density plots. Setting lwd.y=0 suppresses the location line.}
  \item{col.qqline, lty.qqline, lwd.qqline}{Graphical parameters for the Q-Q line in the 
     Q-Q plot.}
}
\details{
For uncert objects created with methods other than MC, the plot types are:
\describe{
  \item{\code{which=1}}{A histogram of the MC replicates in \code{x$MC$y}, 
     with optional line for \code{x$MC$y}. The histogram is produced using
     \code{hist.default}}
  \item{\code{which=2}}{A Q-Q plot of the MC replicates in \code{x$MC$y}, 
     with Q-Q line. The plot uses \code{qqnorm.default}. If \code{datax} is not
     present (in sQuote{\dots}), it is set to \code{TRUE}.}
  \item{\code{which=3}}{A density plot of the MC replicates in \code{x$MC$y}.
     The plot calls \code{density.default} to calculate the density and 
     \code{plot.density} to produce the plot.}
  \item{\code{which=4}}{A bar plot of eqn{cor(x_i,y)} if \code{x$y}
     is present. Any correlation method supported by \code{stats::cor} may 
     be included in \sQuote{\dots} (e.g as \code{method="pearson"}.}
  \item{\code{which=5}}{A bar plot of eqn{cov(x_i,y)} if \code{x$y}
     is present. Any correlation method supported by \code{stats::cov} may 
     be included in \sQuote{\dots} (e.g as \code{method="pearson"}.}
}
Values outside 1:5 are silently ignored.

Parameters in \sQuote{\dots} are passed to the various plot methods or calculations called. 
Only those parameters relevant to a given plot are passed to each calculation or
plotting function, so \sQuote{\dots} can include any parameter accepted by any of the 
functions called.

For the x-y correlation and x-y covariance plot, values in \code{x$cor.xy} are 
used if available. If not, \code{stats::cor} or \code{stats::cov} is called on values 
in \code{x$MC$y} and \code{x$MC$x} if the latter is available 
(i.e. \code{uncertMC} was called with \code{keep.x=TRUE}). If neither 
\code{x$cor.xy} nor \code{x$MC$x} is present, or if \code{method} is 
unknown, the plot is skipped with a warning.

}

\value{
\code{plot.uncertMC} invisibly returns \code{NULL}.
}
\references{
None.
}
\author{
S. L. R. Ellison, \email{s.ellison@lgc.co.uk}
}

\seealso{
\code{\link{uncertMC-class}}, \code{\link[graphics]{hist}}, 
\code{\link[stats]{qqnorm}}, \code{\link[stats]{qqline}}, 
\code{\link[stats]{density}}, \code{\link[stats]{plot.density}}

}
\examples{
  expr <- expression(a/(b-c))
  x <- list(a=1, b=3, c=2)
  u <- lapply(x, function(x) x/20)
  set.seed(403)
  u.invexpr<-uncertMC(expr, x, u, distrib=rep("norm", 3), B=999, keep.x=TRUE )
  par(mfrow=c(2,2))
  plot(u.invexpr, which=1:4, pch=20, method="k") 
                                # method="k" gives Kendall correlation 
}

\keyword{ hplot }
