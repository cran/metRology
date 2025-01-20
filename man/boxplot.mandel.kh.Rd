\name{boxplot.mandel.kh}
\alias{boxplot.mandel.kh}
\title{
Box plot of Mandel's h or k statistics
}
\description{
Produces a box plot of Mandel's statistics, with optional outlier labels
and indicator lines for unusual values.
}
\usage{
	\method{boxplot}{mandel.kh}(x, probs=c(0.95, 0.99), main,  
		xlab=attr(x, "grouped.by"), ylab=attr(x, "mandel.type"),
		separators=FALSE, zero.line=TRUE, ylim,  p.adjust="none", 
		frame.plot = TRUE, horizontal=FALSE, at,
		... , 
		col.ind=1, lty.ind=c(2,1), lwd.ind=1, 
		col.sep="lightgrey", lwd.sep=1, lty.sep=1,
		lwd.zero=1, col.zero=1, lty.zero=1,
		outlier.labels=row.names(x), cex.lab=0.7, col.lab=1, 
		adj=NULL, pos=NULL, srt=0 )
}
\arguments{
  \item{x}{An object of class \code{"mandel.kh"}
}
  \item{probs}{
	Indicator lines are drawn for these probabilities. Note that
	\code{probs} is interpreted as specifying two-tailed probabilities
	for Mandel's h and one-sided (upper tail) probabilities for
	Mandel's k.
}
  \item{main}{
 	a main title for the plot. If missing, the default is
 	\code{ paste(deparse(substitute(x)), " - Mandel's", 
                     attr(x, "mandel.type"), 
                     if(attr(x, "mandel.method") == "robust") "(Robust variant)")}
  }
  \item{xlab}{
  	a label for the x axis; defaults to the \code{grouped.by} attribute for
  	\code{x}.
}
  \item{ylab}{
  	a label for the x axis; defaults to the \code{mandel.type} attribute for
  	\code{x}.
}
  \item{separators}{
	Logical; if \code{TRUE}, separator lines are drawn between groups of values.
}
  \item{zero.line}{
	logical; if \code{TRUE} a horizontal line is drawn at zero.
}
  \item{ylim}{
	the y limits of the plot. For Mandel's k, the default lower limit for y is zero.
}
  \item{p.adjust}{
	Correction method for probabilities. If not \code{"none"}, passed to \code{p.adjust}
	prior to calculating indicator lines. Usually, indicator lines are drawn without 
	correction (that is, with \code{p.adjust="none"}); specifying a p-value correction 
	effectively turns the Mandel's 	statistics into single outlier tests.
}
  \item{frame.plot}{
	Logical; If \code{TRUE} a box is drawn around the plot.
}
  \item{horizontal}{
	  if \code{TRUE} boxes are plotted horizontally and separators, indicators etc
	  adjusted accordingly.
}
  \item{at}{
	  numeric vector giving the locations where the boxplots should
          be drawn; defaults to \code{1:n} where \code{n} is the number of boxes.
}
  \item{\dots}{
	Other (usually graphical) parameters passed to \code{barplot}. Note that some 
	parameters appear after \dots to prevent spurious argument matching inside
	\code{barplot.default}.
}
  \item{col.ind, lty.ind, lwd.ind}{
	Graphical parameters used for the indicator lines, recyckled to \code{length(probs)}. 
	For \code{attr(x, "mandel.type")=="h"} the graphical parameters are applied to negative
	as well as positive indicator lines, applied outwards from zero.
}
  \item{col.sep, lwd.sep, lty.sep}{
	Graphical parameters used for the separator lines.
}
  \item{lwd.zero, col.zero, lty.zero}{
	Graphical parameters used for the zero line.
}
  \item{outlier.labels}{
	Either a logical indicating whether outliers should be labelled or a character vector
	of length nrow(x) giving labels. Defaults to row.names(x). 
}
  \item{cex.lab, col.lab}{
	Character size and colour for outlier labels, passed to \code{text} as 
	\code{col} and \code{cex} respectively.
}
  \item{adj, pos}{
	Position of outlier labels relative to outliers; passed to \code{text}.
}
  \item{srt}{
	Label rotation, in degrees, for outlier labels; passed to \code{text}.
}
}
\details{
This plot produces a box plot (using \code{boxplot.default}) of the variables in an object of 
class \code{"mandel.kh"}. 

If labels are specified for outliers (the default), outliers are first located based on the 
locations given by boxplot.default. WARNING: ties may be mislabelled, as the label 
allocated will be the _first_ point at that location.  

Indicator lines are, if requested, drawn as for \code{\link{plot.mandel.kh}}.

Vertical separators are drawn at midpoints of \code{at}. If 

}
\value{
\code{boxplot.mandel.kh} returns the box plot statistics returned by 
\code{boxplot}, invisibly.
}
\author{
S Ellison \email{s.ellison@lgcgroup.com}
}

\seealso{
\code{\link[graphics]{boxplot}} for box plot arguments, and \code{\link[graphics]{text}} for outlier label
location, colour and rotation.
\code{\link{mandel.h}}, \code{\link{mandel.k}}, \code{\link{mandel.kh}},
\code{\link{pmandelh}}, \code{\link{pmandelk}}  for probabilities, quantiles etc.

See \code{\link{plot.mandel.kh}} for the 'classic' Mandel plot.
}
\examples{
   data(RMstudy)

   h <- with(RMstudy, mandel.h(RMstudy[2:9], g=Lab))
   boxplot(h, las=2) 
   	#Recall that for normally distributed data mandel's h should 
   	#have the same dispersion (and IQR) for all groups. But outliers adversely 
   	#affect the estimate of dispersion, so the interquartile ranges differ.
   	#The same effect also accounts for the many boxplot outliers visible
   	#inside the classical Mandel indicator lines; the indicators also 
   	#assume normality.
   	
   #with separators:
   boxplot(h, las=2, separators=TRUE)
   
   #With different labels and label colours:
   boxplot(h, las=2, outlier.labels=paste(1:nrow(h)), col.lab=1:5) 
   
   #... and a horizontal variant (note use of pos to change label positions)
   par(omd=c(0.1,1,0,1))		#to make room for axis labels
   boxplot(h, las=1, separators=TRUE, horizontal=TRUE, pos=1)

}

\keyword{ hplot }
