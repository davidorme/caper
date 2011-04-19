\name{brunch}
\alias{brunch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparative analysis using the brunch algorithm.}
\description{
Calculate a linear model using the brunch algorithm.
}
\usage{
brunch(formula, data, phy, names.col, stand.contr = TRUE, ref.var=NULL, node.depth=NULL, equal.branch.length=FALSE)
}


\arguments{
\item{formula}{A model formula.}
\item{data}{An 'comparative.data' object. Alternatively, a data frame.}
\item{phy}{An object of class 'phylo', required when data is not a 'comparative.data' object.}
\item{names.col}{A name specifying the column in 'data' that matches rows to tips in 'phy', required when data is not a 'comparative.data' object.}
  \item{stand.contr}{A logical flag indicating whether or not to standardize contrasts}.
  \item{ref.var}{A reference variable present in the model that is used to specify the direction of calculation of contrasts. If null, this is assumed to be the first explanatory variable.}
  \item{node.depth}{A positive integer greater than 1 used to restrict the model to contrasts with a node depth less than or equal to the specified depth. Tips have a depth of 1.}
  \item{equal.branch.length}{If set to 'TRUE' then all branch lengths are set to 2.}
}

\details{
This function implements the 'brunch' algorithm for modelling the relationship between variables that are phylogenetically non-independent. This method was described and previously implemented in the Mac Classic computer programs CAIC, written by Purvis and Rambaut (1995) and updated by Nick Isaac and Paul-Michael Agapow. 

The 'brunch' algorithm calculates contrasts for models that include binary categorical variables. Contrasts are identified and calculated for all variables in the model for a set of nodes where each side can be unequivocally attributed to one or other of the categories. Unlike 'crunch', nested contrasts are not calculated and each row of data at the tips is used only once. This follows Burt (1989): contrasts whose paths do not meet or cross at any point will be phylogenetically independent.

Factors with more than two levels are supported but *must* be ordered to allow sensible contrasts to be drawn. In addition, there is no single best compromise set of contrasts with non-binary factors and implementations may differ in the set chosen.

The user provides a comparative dataset. The formula specifies the model to be fitted and contrasts are calculated in those variables. The specified reference variable is used to ensure that contrasts for multivariate models are calculated in a consistent direction at each node. The function \code{brunch} acts as a data preparation wrapper for the function \code{contrCalc}, which is not intended to be directly called by users. Missing data can be present in the explanatory variables: the algorithm makes use of the complete data available at each node as was the case with CAIC.

Polytomies - more detail here The Mac Classic program CAIC used 1 for both 'Brunch' and 'Crunch' analyses and this the default.


}
\value{
  A object of class 'caic'.
}
\references{
Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent contrasts (CAIC): an Apple Macintosh application for analysing comparative data.  Computer Appl. Biosciences 11, 247-251.

Burt, A. (1989). Comparative methods using phylogenetically independent contrasts. Oxford Surveys in Evolutionary Biology, 6:33–53.
}

\author{David Orme}

\seealso{\code{\link{caic-class}} for 'caic' object structure and methods.}

\examples{
data(perrisodactyla)
perriso <- comparative.data(perrisodactyla.tree, perrisodactyla.data, Binomial)
brunchMod <- brunch(log.female.wt ~ Territoriality, data=perriso)
summary(brunchMod)

# plot the contrasts
brunchTab <- caic.table(brunchMod)
plot(log.female.wt ~ Territoriality, brunchTab)

# for the actual model diagnostics
par(mfrow=c(3,1))
caic.diagnostics(brunchMod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
