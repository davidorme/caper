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
  \item{data}{A data frame containing the variables specified in the model formula. The data frame can contain missing data for explanatory variables.}
  \item{phy}{An object of class 'phylo' giving the phylogenetic structure in the data.}
  \item{names.col}{A name specifying the column in 'data' that matches rows to tips in 'phy'.}
  \item{stand.contr}{A logical flag indicating whether or not to standardize contrasts}.
  \item{ref.var}{A reference variable present in the model that is used to specify the direction of calculation of contrasts. If null, this is assumed to be the first explanatory variable.}
  \item{node.depth}{A positive integer greater than 1 used to restrict the model to contrasts with a node depth less than or equal to the specified depth. Tips have a depth of 1.}
  \item{equal.branch.length}{If set to 'TRUE' then all branch lengths are set to 2.}
}

\details{
This function implements the 'brunch' algorithm for modelling the relationship between variables that are phylogenetically non-independent. [The method was first described by Felsenstein (1985) and subsequently extended to permit the use of phylogenies with polytomies by Pagel (1992) - correct this...]. This method was previously implemented in the Mac Classic computer programs CAIC, written by Purvis and Rambaut (1995) and updated by Nick Isaac and Paul-Michael Agapow.

The user provides a data frame and a phylogeny. These are matched using \code{names.col} to identify a column in the data frame that contains tip labels on the phylogeny. The formula specifies the model to be fitted and contrasts are calculated in those variables. The specified reference variable is used to ensure that contrasts for multivariate models are calculated in a consistent direction at each node. The function \code{caic} acts as a data preparation wrapper, finding the union between the phylogeny and data and reducing the dataset to the correct model and response matrices. The prepared data is then sent to the \code{contrCalc} engine, which is not intended to be called by users.

The resulting table of contrasts is then used to fit the specified model - note that the intercept is automatically dropped from the model if present, following REF HERE. 

Missing data notes here...

Polytomies - more detail here The Mac Classic program CAIC used 1 for both 'Brunch' and 'Crunch' analyses and this the default.


}
\value{
  A list of class 'caic' containing the following:
	\item{contrast.data}{ A list of the following:
	\describe{
		\item{contr}{A list containing matrices of the contrasts in the response variable (contr\$response) and explanatory variables (contr\$explanatory).}
		\item{nodalVals}{A list containing matrices of the nodal values in the response variable (contr\$response) and explanatory variables (contr\$explanatory).}
		\item{contrVar}{A numeric vector of the calculated variances for each contrast.}
		\item{nChild}{A vector showing the number of nodes descending from each internal node}
		\item{nodeDepth}{A vector showing the maximum number of nodes between each internal node and the tips of the phylogeny (including both the node in question and the tip and hence always >=2) }
		\item{validNodes}{A logical vector showing which internal nodes on the tree have valid contrasts, given the available data and any user constraints.}
	}}
  \item{phy}{A 'phylo' object describing the phylogeny used to calculate contrasts. This is the phylogeny in the original call reduced to those tips that are also named in the provided data but possibly including tips with missing data.}
  \item{lm}{An 'lm' object containing a regression model through the origin for the calculated contrasts}
}
\references{
Felsenstein, J.  (1985).  Phylogenies and the comparative method.  Am. Nat.  125, 1-15

Pagel, M. D. (1992). A method for the analysis of comparative data.  J. theor. Biol. 156, 431-442.

Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent contrasts (CAIC): an Apple Macintosh application for analysing comparative data.  Computer Appl. Biosciences 11, 247-251.
}

\author{David Orme}

\seealso{\code{\link{caic.table}} for extracting a data table, \code{\link{summary.caic}} for a data and model summary and \code{\link{caic.diagnostics}} for model checking.}
\examples{

data(SmallTree)
brunchMod <- brunch(othervar ~ catX2, phy=diTree, data=SmallTreeDat, names.col=tip)

summary(brunchMod)

# plot the contrasts
brunchTab <- caic.table(brunchMod)
plot(othervar ~ catX2, brunchTab)

# for the actual model diagnostics
par(mfrow=c(3,1))
caic.diagnostics(brunchMod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
