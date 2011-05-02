\name{anova.caic}
\alias{anova.caic}
\alias{anova.caiclist}
\alias{logLik.caic}
\alias{plot.caic}
\alias{predict.caic}
\alias{residuals.caic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Anova and model checking methods for independent contrast models.
}
\description{
These functions provide ANOVA tables and model comparison using ANOVA and AIC, along with standard model diagnostic plots and accessor functions for phylogenetic independent contrast objects.
}
\usage{
\method{anova}{caic}(object, ...)
\method{anova}{caiclist}(object, ..., scale=0, test='F')
\method{logLik}{caic}(object, ...)
\method{plot}{caic}(x, ...)
\method{predict}{caic}(object, ...)
\method{residuals}{caic}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class 'caic'.}
  \item{scale}{A character string specifying the test statistic to be used. Can be one of "F", "Chisq" or "Cp", with partial matching allowed, or NULL for no test.}
  \item{test}{numeric. An estimate of the noise variance sigma^2. If zero this will be estimated from the largest model considered.}
  \item{x}{An object of class 'caic'.}
  \item{\dots}{Further argument to be passed to methods.}
}
\details{
The 'anova' method provides access to single anova tables for a model and to comparison of lists of models. The 'logLik' method provides access to the log likelihood of the 'caic' model and hence to AIC comparison of models.

The 'plot' method uses the standard set of model diagnostic plots for linear models. It is also wise to check the evolutionary assumptions of independent contrast models using the 'caic' specific diagnostic plots. The 'predict' and 'residuals' functions provide access to these parts of the 'caic' object.
}
\author{
David Orme
}
\seealso{
\code{\link{crunch}}, \code{\link{brunch}},\code{\link{macrocaic}},\code{\link{caic.diagnostics}}
}
\examples{
data(shorebird)
shorebird.data$lgEgg.Mass <- log(shorebird.data$Egg.Mass)
shorebird.data$lgM.Mass <- log(shorebird.data$M.Mass)
shorebird.data$lgF.Mass <- log(shorebird.data$F.Mass)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)

cMod1 <- crunch(lgEgg.Mass ~ lgM.Mass * lgF.Mass, data=shorebird)
cMod2 <- crunch(lgEgg.Mass ~ lgM.Mass + lgF.Mass, data=shorebird)
cMod3 <- crunch(lgEgg.Mass ~ lgM.Mass , data=shorebird)

anova(cMod1, cMod2, cMod3)
AIC(cMod1, cMod2, cMod3)

plot(cMod3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}
\keyword{stats}% __ONLY ONE__ keyword per line
