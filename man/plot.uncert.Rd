\name{plot.uncert}
\Rdversion{1.1}
\alias{plot.uncert}
\title{
Plot method for 'uncert' objects
}
\description{
Plots for uncertainty budgets produced by \code{uncert} calls.
}
\usage{
\method{plot}{uncert}(x, which = c(1,2,4,5), main = paste(deparse(substitute(x))), 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    caption = list("Variance and covariance contributions", 
    expression(sqrt(group("|", "Variance and covariance contributions", "|"))), 
    expression("Contribution " * u[i](y) == c[i] * u[i]), 
    "Combined contribution", "Correlation (x,y)", 
    "Covariances (x,y)"), cex.caption = 1, ...)
}
\arguments{
  \item{x}{An object of class \code{uncert} produced by a call to \code{uncert()}.
}
  \item{which}{Integer in 1:6; the particular variant(s) of plot required. A 
     vector is permitted, in which case plots are produced in ascending order.
}
  \item{main}{Main title for the plot
}
  \item{ask}{logical; if 'TRUE', the user is _ask_ed before each plot, see 'par("ask=")'}
  \item{caption}{A list of captions for \emph{all} the different plots. 
}
  \item{cex.caption}{Text size for captions. Note that if the number of figures per 
    page is over 2, captions are further scaled by 0.8 
}
  \item{\dots}{Further parameters passed to \code{plot()} (mostly barplot).
}
}
\details{
For uncert objects created with methods other than MC, the plot types are:
\describe{
  \item{\code{which=1}}{A barplot of all non-zero contributions to the combined 
    uncertainty. These are derived from the covariance matrix and the coefficients \eqn{c_i}. 
    For terms on the diagonal of the covariance matrix, these are \eqn{(c_i*u_i)^2}; for 
    off-diagonal terms (the correlation terms), \eqn{2*(c_i*u_i)*(c_j*u_j)*r_ij}. 
    The threshold for deciding an off-diagonal term is nonzero is that its magnitude 
    is greater than \code{2*.Machine$double.eps}. Note that off-diagonal contributions
    may be negative.}
  \item{\code{which=2}}{As for \code{which=1} except that the square root of the 
    absolute value is plotted. For the 'diagonal' terms, these are just \eqn{u_i(y)} 
    in the nomenclature used by the GUM.}
  \item{\code{which=3}}{A barplot of \eqn{u_i(y)^2}, without the correlation terms.}
  \item{\code{which=4}}{A barplot of the sum of all (co)variance contributions 
    associated with each \eqn{x_i}{x[i]}, that is, 
    \deqn{contrib(i)=(c_{i}u_{i})^2+\sum\nolimits _{j\neq{}i} 2(c_{i}u_{i})(c_{j}u_{j})r_{i,j}}{contrib[i]=(c[i]u[i])^2+sum( 2(c[i]u[i])(c[j]u[j])r[i,j], j != i )}. 
    }
  \item{\code{which=5}}{A barplot of the theoretical correlations 
  	\eqn{cov(x_i,y)/u(x_i)u_y}{cov(x[i],y)/(u[i]u.y)}. 
   }
  \item{\code{which=6}}{A barplot of the theoretical covariances \eqn{cov(x_i,y)}{cov(x[i],y)}. 
   }
}
Values of \code{which} outside this range are silently ignored.

For the X-Y correlation and covariance plots, the covariances are calculated from the 
covariance matrix \eqn{V} (supplied to \code{uncert()} as \code{cov} 
or calculated as \code{outer(u,u,"*")*cor}) and sensitivity coefficients \eqn{c_i}{c[i]} 
as \eqn{cov(x_i,y) = \sum\nolimits _{j} V_{j,i}c_j }{cov(x[i],y) = sum_j(V[j,i]c[j]) }. 
In fact the calculation used is simpler: \code{cov.xy <- V \%*\% ci}. The correlations 
are calculated in turn from these using \eqn{cov(x_i,y)/u(x_i)u_y}{cov(x[i],y)/(u[i]*u.y)}.

Perhaps the most informative plots are for \code{which=1}, \code{which=2},
\code{which=4} and \code{which=5}. The first of these includes all nonzero signed contributions, 
making the negative contributions visible; the second (\code{which=2}) makes direct
comparison of magnitudes easier. The combined contribution plot is the effect on
the total variance of removing all terms associated with a particular variable; it
shows how much \eqn{u_y^2}{u(y)^2} would \emph{reduce} if the uncertainty for \eqn{x_i}{x[i]} were 
reduced to zero. Note that in some cases with negative correlation the combined uncertainty can \emph{increase},
on dropping a variable, shown by a negative reduction in the plot. (\code{which=5}) is among the most 
direct indications of the relative importance of individual parameters.

Objects created with the MC method are passed to \code{\link{plot.uncertMC}}.
}
\value{Invisibly returns the default return value for the last plot produced.}
\references{
None.
}
\author{
S. L. R. Ellison, \email{s.ellison@lgcgroup.com}
}

\seealso{
  \code{\link{uncert}}, \code{\link[graphics]{barplot}}, \code{\link{plot.uncertMC}}.
}
\examples{
  
  #An example with negative correlation
  x <- list(a=1, b=3, c=2, d=11)
  u <- lapply(x, function(x) x/10)
  u.cor<-diag(1,4)
  u.cor[3,4]<-u.cor[4,3]<- -0.5
  u.form.c<-uncert(~a+b*2+c*3+d/2, x, u, method="NUM", cor=u.cor)
  
  par(mfrow=c(3,2))
  plot(u.form.c, which=1:6, las=1, horiz=TRUE)  #Note use of barplot parameters

}
\keyword{ hplot }
