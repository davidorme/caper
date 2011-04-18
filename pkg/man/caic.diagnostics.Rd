\name{caic.diagnostics}
\alias{caic.diagnostics}
\alias{print.caic.diagnostics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic plots of independent contrasts models}
\description{
Plots }
\usage{
caic.diagnostics(caicObj, which.terms=NULL, which.plots=c("NV","SD","AGE"), outlier.val=3, label=FALSE, ultrametric.tol=0.0001, ask, test.signif=TRUE, plot.signif=TRUE, alpha=0.05, ...)
\method{print}{caic.diagnostics}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{caicObj}{An object of class 'caic' produced by either \code{caic()} or \code{macrocaic()}}
  \item{which.terms}{A character vector of the terms in the linear model for which diagnostic plots should be generated. The default is to plot all terms.}
  \item{which.plots}{A character vector of the plots to be generated. The choices available plot the absolute contrast values against estimates of nodal values ("NV") and standard deviation ("SD") as well as, for ultrametric trees, the natural log of the age of the node ("AGE").}
  \item{outlier.val}{Contrasts with studentized residuals greater than this value are shown as open circles.}
  \item{label}{Logical value indicating whether the points should be labelled by node number.}
  \item{ultrametric.tol}{The tolerance used to determine whether the phylogeny is ultrametric.}
  \item{ask}{As in \code{par}, this is used to require user confirmation between successive plots.}
  \item{test.signif}{By default, the function tests for the significance of the relationships shown in the plots, but this behaviour can be suppressed by setting this to FALSE.}
  \item{plot.signif}{By default, significant relationships are indicated on the plots by showing the model line.}
  \item{alpha}{This value sets the significance at which model predictions are added to the diagnostic plots.}
  \item{x}{An object of class 'caic diagnostics' to be printed out}
  \item{\dots}{Arguments to be passed to \code{plot}, currently these arguments are not used by \code{print.caic.diagnostics}. }


}
\details{
This function generates plots of absolute values of calculated contrast against estimates of the nodal values, standard deviations and log age of each node. The user can select which of these plots are to be generated and also which terms in the contrasts model will be plotted. By default, the function also fits a linear model for each combination. These models are returned by the function and have a simple print method displaying the slopes. Significant models (using the specified significance level) are shown on the plotted graphs. In addition, significant outliers from the overall model, as judged using studentized residuals, are plotted as open circles.
}
\value{
If \code{test.signif} is TRUE, a list with class `caic.diagnostics' containing, for each of the model terms in \code{which.terms}, a list of \code{lm} objects testing the relationship between the values in each of the specified \code{which.plots}. Otherwise, the command only produces the plots and no value is returned.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ David Orme }
\seealso{\code{\link{crunch}}, \code{\link{brunch}}, \code{\link{macrocaic}}}
\examples{

par(mfrow=c(3,3))
modDiag <- caic.diagnostics(mod)
print(modDiag)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ models }% __ONLY ONE__ keyword per line
