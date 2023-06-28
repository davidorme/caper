\name{pgls-methods}
\alias{coef.pgls}
\alias{residuals.pgls}
\alias{fitted.pgls}
\alias{predict.pgls}
\alias{print.pgls}
\alias{summary.pgls}
\alias{print.summary.pgls}
\alias{nobs.pgls}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic model methods for 'pgls' models.
}
\description{
These are simple summary methods, accessor functions and summary and print methods for 'pgls' models.
}
\usage{
\method{coef}{pgls}(object, ...)
\method{residuals}{pgls}(object, phylo = FALSE, ...)
\method{fitted}{pgls}(object, ...)
\method{predict}{pgls}(object, newdata, ...)
\method{summary}{pgls}(object, ...)
\method{print}{pgls}(x, digits = max(3, getOption("digits") - 3), ...)
\method{print}{summary.pgls}(x, digits = max(3, getOption("digits") - 3), ...)
\method{nobs}{pgls}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class 'pgls'.}
  \item{x}{An object of class 'pgls'.}
  \item{phylo}{Return phylogenetically corrected residuals or ordinary residuals (see details).}
  \item{newdata}{Alternative data for predicting from 'pgls' models.}
  \item{digits}{Number of digits to show in summary methods.}
  \item{...}{Further arguments to methods.}
}
\details{
Phylogenetically corrected residuals from 'pgls' models [TODO].

Note that the r^2 values reported by \code{summary.pgls} have a specific interpretation. \code{pgls} fits the intercept-only model for the data using _exactly_ the same covariance matrix (phylogeny plugged through any branch length transformations) as the fitted model to get a null model. The r-squared and adjusted r-squared that are reported therefore hold the covariance matrix constant, so show percentage of variance explained between a null model and the actual model given that precise model of trait change.

The actual ML null model for the data (optimising the BL transformation independently) might be different from this - but then the r squared values confound change in explanatory power from changing the model parameters and from changing the trait model.
}
\value{
The 'summary' method returns an object of class 'summary.pgls' containing:
\item{call}{The original function call creating the model.}
\item{df}{A vector of the degrees of freedom used to estimate parameters and the residual degrees of freedom.}
\item{sigma}{The square root of the estimated variance of the random error.}
\item{residuals}{The phylogenetically corrected residuals.}
\item{coefficients}{A table of model coefficient, standard errors and t values.}
\item{param}{A vector of branch length parameters used in the model.}
\item{mlVals}{A vector showing which branch length parameters have been optimised.}
\item{param.CI}{A list of length three containing confidence intervals and p values on parameter bounds for each parameter.}
\item{fstatistic}{A vector of the F value, numerator and denominator degrees of freedom for the model.}
\item{r.squared}{The r^2 for the model.}
\item{adj.r.squared}{The adjusted r^2 for the model.}
}
\author{
Rob Freckleton, David Orme
}
\seealso{
\code{\link{pgls}}
}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird) 
print(mod1)

mod1.sum <- summary(mod1)
print(mod1.sum)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}
\keyword{stats}% __ONLY ONE__ keyword per line
