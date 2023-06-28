\name{summary.caic}
\alias{summary.caic}
\alias{print.caic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize a crunch, brunch or macrocaic analysis}
\description{
The summary method simply returns the linear model summary from the 'caic' object. The print method prints some basic information about the analysis followed by the model summary.}
\usage{
\method{summary}{caic}(object, ...)
\method{print}{caic}(x,...)
}
\arguments{
  \item{object}{An object of class 'caic'.}
  \item{x}{An object of class 'caic'.}
  \item{\dots}{Arguments to be passed to 'summary.lm'.}
}
\value{
The summary method returns an object of class 'summary.lm'.
}

\author{David Orme}
\seealso{\code{\link{crunch}},\code{\link{brunch}}, \code{link{macrocaic}}}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)
crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data=shorebird)
print(crunchMod)
summary(crunchMod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
