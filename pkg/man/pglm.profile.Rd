\name{pglm.profile}
\alias{pglm.profile}
\alias{plot.pglm.profile}
\alias{pglm.confint
}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood profiles and confidence intervals for 'pglm' models.
}
\description{
These functions create likelihood profiles for branch length transformations in phylogenetic generalised least squares models and fit confidence intervals to estimated branch length parameters.
}
\usage{
pglm.profile(pglm, which = c("lambda", "kappa", "delta"), N = 50, param.CI = NULL)
pglm.confint(pglm, which=c('lambda','kappa','delta'), param.CI=0.95)
\method{plot}{pglm.profile}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pglm}{A \code{pglm} object.}
  \item{which}{A choice of which branch length transformation ('lambda', 'kappa' or 'delta') to use.}
  \item{N}{The number of points used to profile the likelihood}
  \item{param.CI}{A p value used to add confidence intervals to a likelihood profile for a parameter.}
  \item{x}{A 'pglm.profile' object to plot.}
  \item{...}{Further arguments to plot functions.}
}
\details{
The 'pglm.profile' function calculates the likelihood of a 'pglm' model under different values of branch length transformations. A single parameter is chosen from 'lambda', 'kappa' or 'delta' to be profiled and the model likelihood is calculated at 'N' equally spaced  points between the parameter bounds used in the model. If the model contains a maximum likelihood estimate of the parameter (or if param.CI is not null) then the resulting 'pglm.profile' object will contain estimated confidence intervals. 

Only one parameter is profiled at a time and the other branch length parameters will be held at the fixed or ML estimates used to fit the model. The 'pglm.confint' function is used by either 'pglm' or 'pglm.profile' to find confidence intervals around a maximum likelihood estimate of a given branch length. The model must contain an ML estimate of the parameter for confidence intervals to be calculated.

The plot method simply draws an annotated profile plot, showing the location of the ML estimate and confidence intervals if present.
}
\value{
The 'pglm.profile' function returns a list containing:
\item{x}{Parameter values at which the likelihood has been calculated.}
\item{logLik}{The likelihood value at each value.}
\item{which}{The parameter being profiled.}
\item{pars}{The value of the other fixed parameters.}
\item{dname}{The name of the 'comparative.data' object used to fit the model.}
\item{formula}{The formula of the model being profiled}

If the model contains an ML estimate of the parameter being profiled, then the 'pglm.profile' object will also contain the output of 'pglm.confint':

\item{opt}{The maximum likelihood value of the parameter.}
\item{bounds.val}{The values of the bounds on the parameter.}
\item{bounds.p}{The p value of the likelihood at the bounds, given the ML value.}
\item{ci.val}{The values of the parameter at the confidence intervals.}
\item{ci}{The confidence interval value used.}
}
\author{
David Orme
}
\seealso{
\code{\link{pglm}}
}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
mod <- pglm(log(Egg.Mass) ~ log(M.Mass), shorebird, lambda='ML')
mod.l <- pglm.profile(mod, 'lambda')
plot(mod.l)
pglm.confint(mod, 'delta')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{util}
\keyword{stats}
