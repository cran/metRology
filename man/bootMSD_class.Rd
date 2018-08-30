% file man/bootMSD.Rd
% (c) 2017 S Ellison
\name{ bootMSD-class }
\Rdversion{1.1}
\alias{ bootMSD-class }
\alias{ print.bootMSD}
\alias{ summary.bootMSD}
\alias{ print.summary.bootMSD}
\alias{ barplot.bootMSD}
\alias{ plot.bootMSD}

\title{ Object returned by \code{bootMSD} and associated methods.}
\description{
  The object class returned by \code{\link{bootMSD}} and associated
  \code{print}, \code{summary}, and plotting classes.
}
\usage{
       \method{print}{bootMSD}(x, ...)

       \method{plot}{bootMSD}(x, ...)

       \method{barplot}{bootMSD}(height, ylab="MSD", names.arg=height$labels, 
	crit.vals=TRUE, lty.crit=c(2,1), col.crit=2, lwd.crit=c(1,2), ... )

       \method{summary}{bootMSD}(object, p.adjust="none", ...) 

       \method{print}{summary.bootMSD}(x, digits=3, ..., 
		signif.stars = getOption("show.signif.stars"), 
		signif.legend=signif.stars) 
}

\arguments{
  \item{x}{
            An R object. For \code{print.bootMSD} and \code{plot.bootMSD}, an object 
            of class \code{"bootMSD"}. For \code{print.summary.bootMSD}, an object 
            of class \code{"summary.bootMSD"}.
	}
  \item{height}{
            An object of class \code{"bootMSD"}. 
	}
  \item{object}{
            An object of class \code{"MSD"}. 
	}
  \item{p.adjust}{
	Multiple correction method for calculated \emph{p}-values, passed to
	\code{\link{p.adjust}}. 
	}
  \item{ylab}{
	Label for vertical axis, passed to \code{barplot}
	}
  \item{names.arg}{
	Labels for individual bars in bar plot, passed to \code{barplot}. If \code{names(height)} 
	is \code{NULL}, bars are numbered. 
	}
  \item{crit.vals}{
	If \code{TRUE}, individual critical values based on observation-specific 
	bootstrap quantiles are added to the plot. These are taken from \code{critical.values}
	in the supplied \code{bootMSD} object.
	}
  \item{lty.crit, col.crit, lwd.crit}{
	Vectors of line style parameters for plotted critical values, passed to 
	\code{\link{segments}}. Recycled to the length of \code{critical.values}
	in the supplied \code{bootMSD} object.}
   \item{digits}{integer; passed to \code{print}. The minimum number of
          significant digits to be printed in values. Change to \code{NULL} for default.}
   \item{signif.stars}{logical; if \code{TRUE}, P-values are additionally encoded
          visually as \sQuote{significance stars} in order to help scanning of
          long coefficient tables.  Defaults to the \code{show.signif.stars} 
          slot of \code{options}.}
   \item{signif.legend}{logical; if \code{TRUE}, a legend for the \sQuote{significance
          stars} is printed provided \code{signif.stars == TRUE}.}
   \item{\dots}{Parameters passed to other methods. }
}

\details{
	The default \code{plot} method is an alias for the \code{barplot} method. 
	For the plot methods, quantiles for each point are taken directly from the quantiles 
	calulated by \code{\link{bootMSD}} and retained in the returned object. 
	
	For the \code{summary} method, \emph{p}-values are initially calculated as the observed 
	proportion of simulated values exceeding the MSD value calculated by \code{msd}. The
	summary method additionally returns \emph{p}-values  after adjustment
	for multiple comparisons using the adjustment method specified. 
	
	The \code{print} method for the \code{summary.bootMSD} object prints the summary as a data 
	frame adjusted with columns for the calculated MSD values, data-specific upper quantiles 
	(one column for each probability supplied to \code{bootMSD} and the \emph{p}-values 
	after adjustment for multiple comparisons based on the proportion of simulated values 
	exceeding the observed MSD. Where that proportion is zero, the summary replaces the 
	raw zero proportion with \code{1/B}, corrects that proportion using the requested 
	adjustment method, andreports the \emph{p}-value as less than ("<") the resulting 
	adjusted value. 
	
}
\value{
	The \code{print} method returns the object, invisibly.
	
	The \code{plot} and \code{barplot} methods return the values at the midpoint of each bar.
	
	The \code{summary} method returns an object of class \code{"summary.bootMSD"} which 
	is a list with members:
	\itemize{
		\item{msd}{Calculated MSD values from \code{msd}}
		\item{labels}{character vector of labels for individual data points}
		\item{probs}{Probabilities used for quantiles}
		\item{critical.values}{matrix of quantiles. Each row corresponds to a probability 
			  in \code{probs} and each column to an individual data point.}
		\item{pvals}{\emph{p}-values estimated as the observed proportion of
			simulated values exceeding the MSD value calculated by \code{msd}.}
		\item{p.adjust}{Character value containing the name of the p-value adjustment method used.}
		\item{p.adj}{ \emph{p}-values adjusted using the given \emph{p}-value adjustment method 
			specified by \code{p.adjust}.}
		\item{B}{Number of bootstrap replicates used.}
		\item{method}{The sampling method used by the parametric bootstrap.}
	}
	
}

\seealso{
  \code{\link{msd}}, \code{\link{qmsd}}.
}
\author{
S. L. R. Ellison \email{s.ellison@lgc.co.uk}
}
\examples{

  \dontrun{
  data(Pb)
  msd.Pb<-msd(Pb$value, Pb$u)  # Uses individual standard uncertainties

  set.seed(1023)
  boot.Pb <- bootMSD(msd.Pb)
  summary(boot.Pb)

  # The default summary gives individual observation p-values. To 
  # avoid over-interpretation for the study as a whole, 
  # apply a sensible p-value adjustment:
  summary(boot.Pb, p.adjust="holm")

  plot(boot.Pb, crit=TRUE)

  }
}
\keyword{ univar }
