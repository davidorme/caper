\name{pglm-methods}
\alias{coef.pglm}
\alias{residuals.pglm}
\alias{fitted.pglm}
\alias{predict.pglm}
\alias{summary.pglm}
\alias{print.pglm}
\alias{plot.pglm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic model methods for 'pglm' models.
}
\description{
These are simple summary methods, accessor functions and summary and print methods for 'pglm' models.
}
\usage{
\method{coef}{pglm}(object, ...)
\method{residuals}{pglm}(object, phylo = FALSE, ...)
\method{fitted}{pglm}(object, ...)
\method{predict}{pglm}(object, pred.x, ...)
\method{summary}{pglm}(object, ...)
\method{print}{pglm}(x, ...)
\method{print}{pglm.summary}(x, digits = max(3, getOption("digits") - 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class 'pglm'.}
  \item{x}{An object of class 'pglm'.}
  \item{phylo}{Return phylogenetically corrected residuals or ordinary residuals (see details).}
  \item{pred.x}{Alternative data for predicting from 'pglm' models.}
  \item{digits}{Number of digits to show in summary methods.}
  \item{...}{Further arguments to methods.}
}
\details{
Phylogenetically corrected residuals from 'pglm' models [TODO].
}
\value{
The 'summary' method returns an object of class 'summary.pglm' containing:
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
\code{\link{pglm}}
}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
mod1 <- pglm(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird) 
print(mod1)

mod1.sum <- summary(mod1)
print(mod1.sum)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}
\keyword{stats}% __ONLY ONE__ keyword per line
