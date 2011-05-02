\name{anova.pglm}
\alias{anova.pglm}
\alias{anova.pglmlist}
\alias{logLik.pglm}

\title{
Anova and AIC tables for 'pglm' models.
}
\description{
The 'anova' function creates ANOVA tables for a 'pglm' models using sequential sums of squares. 
}
\usage{
anova.pglm(object, ...)
anova.pglmlist(object, ..., scale = 0, test = "F")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A 'pglm' model object.}
  \item{\dots}{Additional 'pglm' models.}
  \item{scale}{A character string specifying the test statistic to be used. Can be one of "F", "Chisq" or "Cp", with partial matching allowed, or NULL for no test.}
  \item{test}{numeric. An estimate of the noise variance sigma^2. If zero this will be estimated from the largest model considered.}
}
}
\details{
The sequential sums of squares are calculated by refitting the model in the order of the terms of the formula and so can take a little time to calculate. Branch length transformations are held at the values of the initial object. The 'logLik.pglm' provides a simple accessor function that allows the use of AIC model comparisons. Note that the generic AIC methods do no checking to ensure that sensible models are being compared.
}
\value{
A table of class 'anova' and 'data.frame' that employs the generic plot methods for 'anova' tables.}
\author{
David Orme
}
\note{
The functions build heavily on the generic methods 'anova.lm' and 'anova.lmlist'.
}
\seealso{
\code{\link{pglm}}
}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)

mod1 <- pglm(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird) 
anova(mod1)

mod2 <- pglm(log(Egg.Mass) ~ log(M.Mass) + log(F.Mass), shorebird)  
mod3 <- pglm(log(Egg.Mass) ~ log(M.Mass) , shorebird)
mod4 <- pglm(log(Egg.Mass) ~ 1, shorebird)

anova(mod1, mod2, mod3, mod4)
AIC(mod1, mod2, mod3, mod4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
