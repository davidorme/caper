\name{plot.pgls}
\alias{plot.pgls}
\title{
Diagnostic plots for 'pgls' models.
}
\description{
The function generates four diagnostics plots for 'pgls' models.
}
\usage{
\method{plot}{pgls}(x, ...)
}

\arguments{
  \item{x}{An object of class 'pgls'.}
  \item{...}{Additional arguments to plot functions.}
}
\details{
The first two plots show the fit of the phylogenetic residuals from the model to a normal distribution: a density plot of the residuals and a normal Q-Q plot. The second two plots scatterplots show pattern in the distribution of the fitted values against the observed and residual values.
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
par(mfrow=c(2,2))
plot(mod1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}
\keyword{graphics}
