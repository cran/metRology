% file RKV/man/kplot.Rd
% (c) 2009 S Ellison
\name{ kplot }
\Rdversion{1.1}
\alias{ kplot }
\alias{ kplot.default }
\alias{ kplot.ilab }
\alias{ kpoints }

\title{ Dot-and-bar plot for Key Comparisons }
\description{
  Produces a dot plot of typical metrology comparison data (value/uncertainty) with 
  error bars, assigned value and uncertainty and optional percentage deviation axis or 
  marginal density
}
\usage{
kplot(x, \dots)

\method{kplot}{default}(x,U=NULL, labels=names(x),  assigned=NULL, U.assigned=NULL, 
	U.lo=U, U.hi=U, k=2, strata=NULL,
	do.percent=!is.null(assigned) && !do.pdf, 
	ordered=TRUE, order.strata=levels(strata),
	xlim=c(0.5, length(x)+0.5), ylim,
	main=NULL, xlab=NULL, ylab=NULL,
	axis.main=2, axis.pct=4, at=1:length(x), at.main=NULL,
	cex.axis=0.8, las=2, las.pct=1, ylab.line=2.5, 
	ylab.line.pct=2.1, ci.width=0.03, col.ci=par("fg"), 
	lty.ci=par("lty"), lwd.ci=par("lwd"), pch=21, 
	col=par("fg"), bg="white", add.outliers=FALSE, 
	outlier.offset=0.2, mar=NULL, box=TRUE,
	do.pdf=FALSE, do.individual.pdf=do.pdf, 
	col.pdf=par("fg"), lwd.pdf=1, lty.pdf=1,
	do.total.pdf=TRUE, col.total.pdf=col.pdf[1], 
	lwd.total.pdf=2, lty.total.pdf=1, n.pdf=200, 
	pdf.layout=c(4,1), pdf.scale=0.7, pdf.offset=0.05, 
	xlim.pdf, pdf.axis=FALSE, las.pdf=0, 
	mgp.pdf=c(3,0.5,0), \dots)
	
\method{kplot}{ilab}(x, \dots)

kpoints(x,U=NULL, labels=names(x), U.lo=U, U.hi=U, k=2, 
	strata=NULL, ordered=TRUE, order.strata=levels(strata),
	at=1:length(x), ci.width=0.03, col.ci=par("fg"), 
	lty.ci=par("lty"), lwd.ci=par("lwd"), pch=21, 
	col=par("fg"), bg="white", add.outliers=FALSE, 
	outlier.offset=0.2, \dots)

}
\arguments{
  \item{x}{an R object. For the default method, a vector of reported values. For the 
  	ilab method, an object of class \sQuote{ilab}}
  \item{U}{Vector of length \code{length(x)} of expanded uncertainties}
  \item{labels}{Vector of of length \code{length(x)} labels for x-axis marks. }
  \item{assigned}{Assigned value for the comparison. Plotted as a horizontal line on the plot.}
  \item{U.assigned}{Expanded uncertainty for the assigned value}
  \item{U.lo, U.hi}{Vectors of of length \code{length(x)} of lower and upper 
                    limits for the uncertainty intervals around the reported 
                    values, to allow asymmetric intervals. Both default to U. }
  \item{k}{Coverage factor originally used in calculating U. Required only if 
           \code{do.pdf=TRUE}, as k is used to calculate standard uncertainties 
           from U. k can be a scalar (recycled to length  \code{length(x)} 
           if necessary) or a vector of length of length \code{length(x)}. }
  \item{strata}{A Factor identifying subsets of the data. Currently not implemented.}
  
  \item{do.percent}{Logical indicating whether percentage deviation should be 
           plotted as a secondary axis. Defaults to \code{TRUE} if an assigned 
           value is provided and \code{FALSE} if there is no assigned value 
           or if a marginal density is  required via \code{do.pdf=TRUE}.}
  \item{ordered}{If \code{TRUE}, values are plotted in ascending order.}
  \item{order.strata}{Character vector showing the order of plotting for strata. 
            Currently not implemented.}
  \item{xlim, ylim}{Plot limits as in \code{plot.default}.}
  \item{main, xlab, ylab}{Titles; see \code{?title} for details.}
  \item{axis.main, axis.pct}{Integers specifying on which side of the plot the relevant axis is to
          be drawn, passed to \code{axis} as \code{side}. The axis is placed as for follows: 
          1=below, 2=left, 3=above and 4=right. The main axis (\code{axis.main}) is provided 
          in the same units as \code{x}; the percentage axis (\code{axis.pct}) shows 
          corresponding percentage deviation.}
  \item{at}{Vector of x-axis locations for the data points, x-axis tick marks and labels. 
            Defaults to \code{1:length(x)}. }
  \item{at.main}{The points at which tick-marks are to be drawn on the main (y) axis. 
            Passed to \code{axis} via \code{at}. }
  \item{cex.axis}{The magnification to be used for axis annotation
          relative to the current setting of 'cex'. Passed to \code{axis}.}
  \item{las, las.pct}{Integers defining x- and y axis and percentage axis label
            orientation; see \code{par(las)}.}
  \item{ylab.line, ylab.line.pct}{Margin lines for main and percentage axis titles. 
  	Passed to \code{axis}.}
  \item{ci.width}{Width of error bar terminators, passed to \code{arrows}.}
  \item{col.ci, lty.ci, lwd.ci}{Graphical parameters for the error bars; passed 
            to \code{arrows}. }
  \item{pch, col, bg}{Graphical parameters for data points, passed to \code{points}. 
             \code{bg} specifies the fill colour for \code{pch} from 21 to 25. }
  \item{add.outliers}{If \code{TRUE}, points outside \code{ylim} are indicated as 
            an arrow indicating the direction in which the omitted points lie, with 
            a text label showing the reported value.}
  \item{outlier.offset}{X-offset (in x-axis units) specifying lateral location of outlier 
            tet labels relative to x-axis location of the outlier indicator.}
  \item{mar}{Plot margins as in \code{par(mar)}. The default varies depending on 
             \code{do.pct} and \code{do.pdf}.}
  \item{box}{If \code{TRUE}, a box is drawn round the plot region.}
  \item{do.pdf}{If \code{TRUE}, a marginal density and/or individual densities for the 
                individual reported values are plotted based on the reported 
                values \code{x} and standard uncertainties calculated as \code{U/k}.}
  \item{do.individual.pdf}{Logical controlling whether the individual densities 
                 are plotted as well as/instead of the combined density.}
  \item{col.pdf, lwd.pdf, lty.pdf}{Graphical parameters controlling the appearance of 
                 the marginal density plot(s). Vectors are permitted, allowing different 
                 styles for each individual pdf.}
  \item{do.total.pdf}{Logical controlling whether the sum of individual densities 
                 is plotted.}
  \item{col.total.pdf, lwd.total.pdf, lty.total.pdf}{Graphical parameters controlling the 
                 appearance of the marginal density plot for the combined density. }
  \item{n.pdf}{Number of points used to construct the marginal density.}
  \item{pdf.layout}{Vector of length 2 specifying the relative sizes of the main plot and 
                 marginal density plot. See \code{?layout} for details.}
  \item{pdf.scale, pdf.offset}{Offset and scaling factor used to control the location 
                 and height of the marginal density plot(s). }
  \item{xlim.pdf}{Controls the x-axis (i.e. the horizontal axis) for the marginal density 
                  plotting area. }
  \item{pdf.axis}{If \code{TRUE} and no other axis has been plotted to the right of the main plot, 
                 an axis is plotted with the marginal density.}
  \item{las.pdf, mgp.pdf}{Axis control parameters passed to \code{axis} to plot the axis for the 
                 marginal density.}
 \item{\dots}{Parameters passed to other functions; currently unused.}
}
\details{
	
	If \code{do.pdf=TRUE} a marginal density plot is added. This plot is constructed 
	from a set of (currently) normal densities centred at \code{x} with standard 
	deviation \code{U/k}.
	
	If a marginal density is plotted, \code{par("layout")} is changed to 
	\code{pdf.layout}; otherwise, \code{par("layout")} is set to \code{matrix(1)}. 
	Both override any previously set layout. \code{par("layout")} is preserved on exit.
	
	The \sQuote{ilab} method passes all parameters in \sQuote{\dots} to the default method
	with default values for \code{x}, upper and lower bounds \code{U.lo} and \code{U.hi}, 
	labels and title taken from the \code{ilab} object.
	
	\code{kpoints} is a convenience function for adding points with confidence 
	intervals to an existing plot. \code{kpoints} is not a generic function  
	and requires a vector \code{x}. Note that \code{kpoints} does not check for 
	the presence of a marginal density plot.
}
\value{
	Invisibly returns a list with components:
	\item{order}{The order for plotting the original data, as returned by \code{order}. }
	\item{at}{ x-axis locations used, in plotting order. }
	
}
\author{S Ellison \email{s.ellison@lgcgroup.com}}
\seealso{
  \code{\link{arrows}}.
}
\examples{

  data(Pb)
  kplot(Pb$value, Pb$U, assigned=2.99, U.assigned=0.06)
  kplot(Pb$value, Pb$U, assigned=2.99, U.assigned=0.06, do.pdf=TRUE)
	
  #Use of return value for annotation
  kp<-kplot(Pb$value, Pb$U, assigned=2.99, U.assigned=0.06)
  text(kp$at, Pb$value-Pb$U, Pb$lab, srt=90, pos=4, cex=0.7)

}

\keyword{ hplot }
