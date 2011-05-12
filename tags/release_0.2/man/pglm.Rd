\name{pglm}
\alias{pglm}
\alias{pglm.likelihood}
\alias{pglm.blenTransform}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phylogenetic generalized linear models}
\description{
Fits a linear model, taking into account phylogenetic non-independence between data points. The strength and type of the phylogenetic signal in the data matrix can also be accounted for by adjusting branch length transformations (lambda, delta and kappa). These transformations can also be optimised to find the maximum likelihood transformation given the data and the model.}
\usage{
pglm(formula, data, lambda = 1.0, kappa = 1.0,  delta= 1.0, param.CI = 0.95, control = list(fnscale=-1), bounds = list(lambda=c(1e-6,1), kappa=c(1e-6,3), delta=c(1e-6,3)))
pglm.likelihood(optimPar, fixedPar, y, x, V, optim.output=TRUE, names.optim=NULL)
pglm.blenTransform(V, fixedPar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{A model formula}
\item{data}{A 'comparative.data' object containing the covariance matrix and data to be used in the model.}
\item{lambda}{A value for the lambda transformation.}
\item{kappa}{A value for the kappa transformation.}
\item{delta}{A value for the delta transformation.}
\item{param.CI}{A p value used to calculate confidence intervals.}
\item{control}{A list of control parameters for the optim function.}
\item{bounds}{A list of length three giving upper and lower bounds for branch length parameters.}
\item{optimPar}{A named vector of branch length parameters to be optimised to find the maximum likelihood value.}
\item{fixedPar}{A named vector of fixed values for branch length parameters.}
\item{y}{A column matrix of the model response.}
\item{x}{The design matrix of the model.}
\item{V}{A phylogenetic covariance matrix.}
\item{optim.output}{A logical value. If true then 'pglm.likelihood' returns only the likelihood value for use in the 'optim' function.}
\item{names.optim}{The name of a single parameter being optimised. This is only required for estimating parameter confidence intervals, where the function 'uniroot' strips names from vectors.}
}
\details{
This function fits a linear model controlling for the non-independence between cases resulting from phylogenetic structure in the data. The stucture of the phylogenetic signal can be controlled by altering the parameters lambda, delta and kappa (see the 'caper' vignette for details). The implementation of the method is currently as described in Freckleton et al (2002). 

The 'pglm.likelihood' and 'pglm.blenTransform' methods are not primarily intended to be called by users. The 'pglm.likelihood' function provides a general method to calculate the likelihood of a model, given the covariance matrix, response, design matrix and branch length parameters. 
}
\section{Warning}{
The model is fitted using a data frame reduced to complete row cases to eliminate missing values. In order to ensure that the models fitted using different subsets of the data are comparable, the whole data frame \code{data} is reduced to complete cases. In the future, a scope argument may be provided to control this but at present the data frame should be reduced to only those variables used in the maximal model in order to avoid prevent redundant variables causing rows to be dropped unnecessarily.
}

\value{
The 'pglm' function returns an object of class \code{pglm} containing the following:

"na.action" "param.CI"
 \item{call}{The original call to the 'pglm' function}
 \item{model }{A summary of the fitted model containing:}
 \item{formula}{The model formula supplied.}
 \item{data}{The comparative data object provided.}
 \item{dname}{The name of the comparative data object.}
 \item{logLikY}{The log likelihood of the response variable given the model.}
 \item{RMS}{The residual mean square variance in the model.} 
 \item{RSSQ}{The residual sum of squares from the model.} 
 \item{NMS}{The null mean square variance for the model.} 
 \item{NSSQ}{The null sum of squares for the response.} 
 \item{aic}{The AIC score of the model} 
 \item{aicc}{The AICc score of the model, correcting for the number of cases and parameters estimated}
 \item{n}{The number of rows of data used in fitting the model}
 \item{k}{The number of parameter estimates} 
 \item{sterr}{The standard errors of the parameter estimates} 
 \item{Vt}{The phylogenetic covariance matrix used in the model, with branch length transformations applied.} 
 \item{fitted}{The predicted values} 
 \item{residuals}{The non-phylogenetic residuals}
 \item{phyres}{The phylogenetic residuals}
 \item{x}{The design matrix of the model }
 \item{varNames}{The variables include in the model.} 
 \item{y}{The response of the model.}
 \item{namey}{The name of the response variable.}
 \item{param}{A named numeric vector of length three giving the branch length transformations used in the model.}
 \item{mlVals}{A named logical vector of length three indicating which branch length values in 'param' are maximum likelihood estimates.}
 \item{bounds}{The bounds on branch length parameter estimates used in the model.}
 \item{param.CI}{A named list of length three giving confidence intervals and the p values at the parameter bounds for optimised branch length transformations. Fixed parameters will have a NULL entry in this list.}
 \item{na.action}{A named vector identifying any rows of missing data excluded from the model.}
}
\references{R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic analysis and comparative data: A test and 
review of evidence. American Naturalist, 160:712-726, 2002.}
\author{Rob Freckleton; David Orme}
%\note{}
\seealso{\code{\link{pglm.profile}}, \code{\link{anova.pglm}}, \code{\link{summary.pglm}}}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
mod1 <- pglm(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird, lambda='ML')
mod2 <- pglm(log(Egg.Mass) ~ log(M.Mass), data=shorebird, lambda='ML', delta='ML')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
