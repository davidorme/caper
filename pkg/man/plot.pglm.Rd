\name{plot.pglm}
\alias{plot.pglm}
\title{
Diagnostic plots for 'pglm' models.
}
\description{
The function generates four diagnostics plots for 'pglm' models.
}
\usage{
plot.pglm(x, ...)
}

\arguments{
  \item{x}{An object of class 'pglm'.}
}
\details{
The first two plots show the fit of the phylogenetic residuals from the model to a normal distribution: a density plot of the residuals and a normal Q-Q plot. The second two plots scatterplots show pattern in the distribution of the fitted values against the observed and residual values.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utils}
\keyword{graphics}
