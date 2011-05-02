\name{pglm}
\alias{pglm}
\alias{pglm.likelihood}
\alias{pglm.blenTransform}


\alias{coef.pglm}
\alias{residuals.pglm}
\alias{fitted.pglm}
\alias{predict.pglm}
\alias{summary.pglm}
\alias{print.pglm}
\alias{plot.pglm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phylogenetic generalized linear models}
\description{
Fits a linear model, taking into account phylogenetic non-independence between data points. The strength of the phylogenetic signal in the data matrix can also be accounted for by adjusting the parameter lambda}
\usage{
pglm(formula, data, V = NULL, lambda = 1.0, kappa = 1.0,  delta= 1.0, param.CI = 0.95, control = list(fnscale=-1), bounds = list(lambda=c(1e-6,1), kappa=c(1e-6,3), delta=c(1e-6,3)))
pglm.likelihood(optimPar, fixedPar, y, x, V, optim.output=TRUE, names.optim=NULL)
pglm.blenTransform(V, par)

\method{coef}{pglm}(object, ...)
\method{residuals}{pglm}(object, phylo = FALSE, ...)
\method{fitted}{pglm}(object, ...)
\method{predict}{pglm}(object, pred.x, ...)
\method{logLik}{pglm}(object, k, ...)
\method{AICc}{pglm}(object)
\method{lambda}{pglm}(object)
\method{summary}{pglm}(object, ...)
\method{print}{pglm}(x, ...)
\method{plot}{pglm}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A model formula}
  \item{data}{A data frame containing the variables in the model}
  \item{phylomat}{A variance-covariance matrix representing the phylogenetic structure in the data.}
  \item{lambda}{A value between 0 and 1 determining the strength of the phylogenetic pattern in the data.}
  \item{plotit}{A logical indicating whether or not to plot a likelihood surface for the model with varying values of lambda.}
  \item{object}{An object of class 'pglm'.}
  \item{x}{An object of class 'pglm'.}
  \item{phylo}{Logical value indicating whether to extract standard or phylogenetic corrected residuals}
  \item{pred.x}{Values for which to extract predicted values from the model.}
  \item{k}{Value of k to use for calculating AIC - warning: this is currently fixed at 2 and the method needs integrating into the generic}
  \item{\dots}{Further arguments.}
}

\details{
This function fits a linear model controlling for the non-independence between cases resulting from phylogenetic structure in the data. The strength of the phylogenetic signal can be controlled by altering the parameter lambda, which scales the covariance between data points (the off-diagonals of \code{phylomat}, representing the amount of shared evolutionary history) as the product of this shared history and lambda. A lambda of 0 is equivalent to a standard linear model, with all shared phylogenetic history reduced to zero, and a lambda of 1 uses the original covariance matrix. Note that a \code{pglm} model with a lambda of 0 will not be the same as a \code{lm} model if the tree is not ultrametric: although the covariance terms will all be zero, the variances of the tips will differ. The function \code{pglmEstLambda} identifies the maximum likelihood value of lambda for a given model formula and also tests, using a log likelihood ratio test, whether this model is significantly more likely than a lambda of 0 or 1.

The implementation of the method is currently as described in Freckleton et al (2002). This section will be expanded to compare the results obtained in this implementation to that of Pagel's Continuous.

}
\section{Warning}{
The model is fitted using a data frame reduced to complete row cases to eliminate missing values. In order to ensure that the models fitted using different subsets of the data are comparable, the whole data frame \code{data} is reduced to complete cases. In the future, a scope argument may be provided to control this but at present the data frame should be reduced to only those variables used in the maximal model in order to avoid prevent redundant variables causing rows to be dropped unnecessarily.
}

\value{
Returns an object of class \code{pglm} containing the following:
 \item{model }{A summary of the fitted model containing:}
 \item{formula}{The model formula supplied.}
 \item{logLikY}{The log likelihood of the response variable given the model.}
 \item{RMS}{RMS} 
 \item{NMS}{NMS} 
 \item{NSSQ}{}, 
 \item{RSSQ}{}, 
 \item{aic}{The AIC score of the model} 
 \item{aicc}{The AICc score of the model, correcting for the number of cases and parameters estimated}
 \item{n}{The number of rows of data used in fitting the model}
 \item{k}{The number of parameter estimates} 
 \item{sterr}{The standard errors of the parameter estimates} 
 \item{vcv}{} 
 \item{fitted}{The predicted values} 
 \item{residuals}{The non-phylogenetic residuals}
 \item{phyres}{The phylogenetic residuals}
 \item{x}{}
 \item{data}{The data frame provided} 
 \item{varNames}{} 
 \item{y}{}
 \item{V}{} 
 \item{lambda}{The value of lambda used in fitting the model.} 
 \item{L0}{The log-likelihood of the model with lambda = 0.}
 \item{L1}{The log-likelihood of the model with lambda = 1.} 
 \item{LamOptimised}{A logical indicating whether the fitted lambda has been optimised using \code{pglmEstLambda()}}
 \item{namey}{The name of the y variable.}
}
\references{R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic analysis and comparative data: A test and 
review of evidence. American Naturalist, 160:712-726, 2002.}
\author{Rob Freckleton; David Orme}
%\note{}
\seealso{\code{\link{profile.lambda}} }
\examples{

data(shorebird)
shorebird.VCV <- vcv.phylo(shorebird.tree)

shorebird.modL1 <- pglm(Cl.size ~ M.Mass, data=shorebird.data, ph=shorebird.VCV, lambda=1)
shorebird.modL0 <- pglm(Cl.size ~ M.Mass, data=shorebird.data, ph=shorebird.VCV, lambda=0)
shorebird.mod <- pglmEstLambda(Cl.size ~ M.Mass, data=shorebird.data, ph=shorebird.VCV)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
