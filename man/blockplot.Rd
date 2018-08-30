\name{blockplot}
\alias{blockplot.default}
\alias{blockplot.formula}
\alias{blockplot}
\alias{bplot}
\alias{nclass.23}

\title{
Generate a \dQuote{block plot} - a histogram variant identifiying individual data points.
}
\description{
A \dQuote{block plot} is a histogram variant identifiying individual data points. Histogram bars are replaced 
by a stack of rectangles (\dQuote{blocks}, each of which is labelled. \code{blockplot} provides for grouped 
data, which generates vertically separated subplots for each group. Fills and label colours can be specified 
for each data point.
}
\usage{
	blockplot(x, ...)

	bplot(x, ...)

	\method{blockplot}{default}(x, breaks = "23", labels = paste(1:length(x)), groups = NA, 
		xlim = NULL, ylim = NULL, 
		main = NULL, xlab = NULL, ylab = "Frequency", grp.labs = FALSE, 
		include.lowest = TRUE, right = TRUE, nclass = NULL, 
		plot = TRUE,  add=FALSE, ...)

	
	\method{blockplot}{formula}(x, data = NULL, ..., subset, main = NULL, xlab = NULL)
	
	
	nclass.23(x)

}
\arguments{
  \item{x}{ An R object. For the default method, a vector of values for which the blockplot is desired. 
            For the formula method, a valid formula, such as \code{y ~ grp} (see Details).
  }
  \item{data}{ For the formula method, a data frame or list from which the variables in \code{formula}
          should be taken.
  }
  \item{subset}{ For the formula method, an optional vector specifying a subset of observations to be
          used for plotting.
  }
  \item{breaks}{ Either a specification for choosing breakpoints for \dQuote{binning} 
  	the data, or a vector giving the breakpoints themselves. The 
  	specification can be a single number, a function, or a character 
  	string identifying a function. See \sQuote{Details} for detailed 
  	specification. 
  }
  \item{labels}{Labels for data points, passed to \code{text}; in principle of any type acceptable to \code{text}. 
        Labels are placed inside boxes so should be short for readability.
  }
  \item{groups}{An optional grouping variable, coerced to factor. If present, one subplot is produced for each non-empty group.
  }
  \item{xlim}{ Limits \code{(x1, x2)} for the horizontal (\eqn{x}) range of the plot. The default is the range of \code{breaks}, after computation if necessary.
  }
  \item{ylim}{ limits for the vertical range of the plot. Will be overridden if \code{square=TRUE} (see below).
  }
  \item{main}{Main title for the plot, passed to \code{\link{plot}}.
  }
  \item{xlab, ylab}{ x- and y-axis labels for the plot. As usual, either can be expressions (see \code{\link{plotmath}}).
  }
  \item{grp.labs}{ Logical, determining whether group labels are plotted, or 
  	a vector of labels. See \sQuote{Details}.
  }
  \item{include.lowest}{ Logical, indicating whether a value equal to the lowest 
        (or highest, for \code{right = FALSE}) \code{breaks} value should be included 
        in a given bin. Passed to \code{\link{cut}}.
  }
  \item{right}{Logical, indicating whether the bin intervals should be closed on the right (and open 
        on the left) or vice versa. Passed to \code{\link{cut}}.
  }
  \item{nclass}{Suggested number of classes for \code{breaks}; equivalent to a single 
  	numerical value for \code{breaks}. 
  }
  \item{plot}{If \code{FALSE}, no plot is produced. The return value is returned invisibly.
  }
  \item{add}{If \code{TRUE}, the plot is added to an existing figure. 
  }
  \item{\dots}{Further parameters passed to other functions, in particular, \code{\link{bkp}}, 
  	which creates the plot, and \code{\link{plot}}
  }
}
\details{
   \code{blockplot} produces a block plot - a histogram variant identifying individual
   data points. Histogram bars are replaced by a stack of rectangles (\dQuote{blocks}, 
   each of which can be (and by default, is) labelled.

   \code{bplot} is an alias for \code{blockplot}.

   For the formula method, \code{x} is a formula, such as \code{y ~ grp}, in which \code{y} is a numeric vector
   of data values to be split into groups according to the grouping variable \code{grp} (usually a factor). 
   More than one grouping variable can be specified, in which case subplots are produced for each level
   of the interaction between grouping factors.

   The specification for breakpoints, \code{breaks}, is modelled closely on that for 
   \code{\link{hist}}. \code{breaks} can be one of:

   \itemize{
       \item a vector giving the (equally spaced) breakpoints between bins;

       \item a function to compute the vector of breakpoints;

       \item a single number giving the suggested number of bins for the
            blockplot;

       \item a character string naming an algorithm to compute the
            number of cells. Values of \code{"23"} (the default),
            \code{"Sturges"}, \code{"Scott"}, \code{"FD"} and 
            \code{"Freedman-Diaconis"} are currently 
            supported; see below for their effect

       \item a function to compute the number of bins.
   }
   In the last three cases the number is a suggestion only, as the breakpoints 
   will be set to \dQuote{pretty} values. 
   
   The different character string specifications correspond to \dQuote{nclass} 
   functions, including those used by \code{\link{hist}}; see \code{\link{nclass.FD}} 
   for details of those. In addition, the default \code{"23"} corresponds to the 
   function \code{nclass.23}. This is just a wrapper for the one-line expression
   
   \code{ceiling(length(x)^(2/3))},
   
   which appears to provide good results for block plots.
   
   Considerable control of graphical elements is provided by the plotting function 
   \code{\link{bkp}}, which is called by \code{blockplot}. In particular, arguments 
   passed through \code{...} to \code{\link{bkp}} can control:
   \itemize{
   	\item The general shape of the plot, including the asbect ratio of the \dQuote{blocks};
   	\item whether a plot should be added to an existing figure (\code{add})
   	\item the fill colour and shading, the border width, type and colour, and the font size 
   	      and colour of individual blocks;
   	\item the vertical location of the plot in the figure region \code{offset};
   	\item the vertical spacing between multiple plots on the same figure when a grouping 
   	variable is provided (\code{grp.spacing} and \code{grp.at};
   	\item the presence, location and appearance of labels for individual subplots;
   	\item whether axes are plotted on any of the four sides of the plot;
   	\item the appearance or omission of empty groups.
   }
   See \code{\link{bkp}} for further details.  
}

\value{Blockplot currently returns an object of class \code{blockplot}, which is a list with elements:
	\item{x}{The original data}
	\item{groups}{If there is more than one group, a factor of groups for each data point,
		with additional attribute \code{"gname"} containing a default name for the grouping 
		variable(s). \code{groups} is set to \code{NA} if there is only one group.}
	\item{x.left}{Vector of x-coordinates for the left side of each block}
	\item{x.height}{Vector of y-coordinates for each box, relative to the group baseline}
	\item{x.mid}{Vector of x-coordinates for the middle of each block (the text location)}
	\item{x.mid}{Vector of x-coordinates for the middle of each block (the text location)}
}

\references{
ISO 5725-2:1994 \emph{Accuracy (trueness and precision) of measurement methods 
and results -- Part 2: Basic method for the determination of repeatability 
and reproducibility of a standard measurement method.} ISO, Geneva, 1994.
}
\note{
The name \dQuote{block plot} may not be in general use, but the package author 
has been unable to identify either an alternative designation or an original source
for this type of plot. An example - apparently hand drawn - was given in ISO 5725-2:1994 (referenced above).

}

\author{
S Ellison \email{s.ellison@lgcgroup.com}
}


\seealso{
For plotting and control of plot appearance:
\code{link{bkp}} 

For graphical elements:
\code{\link{text}}, \code{\link{rect}}

For specification of \code{breaks}:
\code{link{nclass.Sturges}}, \code{link{nclass.Scott}}, \code{link{nclass.FD}}

}
\examples{

#A simple blockplot
set.seed(55)
x<-rnorm(48, 15)
blockplot(x)

#Aspect ratio control for square blocks
blockplot(x, square=TRUE)

# Specifying groups 
grp <- gl(3, 16)
blockplot(x, groups=grp)

#Formula interface
blockplot(x~grp)

#Vectorised colour specification
blockplot(x~grp, square=TRUE, fill=ifelse(1:48 \%in\% c(15, 23, 24), "gold", "white") )

#Group labelling
blockplot(x~grp, square=TRUE, grp.labs=paste("Level", 1:3), grp.pos=2)

#A missing group
xm <- x
xm[ grp == "2" ] <- NA
blockplot(xm~grp, square=TRUE, grp.labs=paste("Level", 1:3), grp.pos=2)

blockplot(xm~grp, square=TRUE, grp.labs=paste("Level", 1:3), grp.pos=2, drop.unused=FALSE)

}

\keyword{ hplot }

