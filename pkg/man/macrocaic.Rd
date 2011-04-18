\name{macrocaic}
\alias{macrocaic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparative analysis using independent contrasts on species richness data.}
\description{
Macroevolutionary hypotheses about correlates of species richness require testing in a phylogenetic framework in order to avoid phylogenetic autocorrelation. Independent contrasts as described by Felsenstein (1985) are appropriate for explanatory variables in such models but not for species richness as the response variable. This function implements two methods for calculating species richness constrasts described by Agapow and Isaac (2002) and originally implemented in the program MacroCAIC.
}
\usage{
macrocaic(formula, data, phy, names.col, macroMethod = "RRD", stand.contr = TRUE, ref.var = NULL, node.depth = NULL, macroMinSize = 3, equal.branch.length = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula describing a linear model predicting species richness.}
  \item{data}{A data frame containing the variables to be used in the model.}
  \item{phy}{An object of class 'phylo'.}
  \item{names.col}{A name identifying a column in \code{data} that contains the tip labels from \code{phy}.}
  \item{macroMethod}{One of either "RRD" for relative rate difference or "PDI" for proportion dominant index.}
  \item{stand.contr}{A logical flag indicating whether to standardize the contrasts.}
  \item{ref.var}{Identifies a predictor variable used for determining the direction of contrasts.}
  \item{node.depth}{A positive integer greater than 1 used to restrict the model to contrasts with a node depth less than or equal to the specified depth. Tips have a depth of 1.}
  \item{macroMinSize}{A positive integer giving the minimum species richness at a node for contrasts to be included in the model.}
  \item{equal.branch.length}{If set to 'TRUE' then all branch lengths are set to 2.}
}
\details{
  More details to come soon.
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

The list has attributes (DETAILS)
}
\references{
Felsenstein, J.  (1985).  Phylogenies and the comparative method.  Am. Nat.  125, 1-15

Agapow, P.-M. and Isaac, N. J. B. (2002) MacroCAIC: correlates of species richness. Diversity & Distributions, 8, 41-43
}
\author{David Orme}
\seealso{\code{\link{crunch}}, \code{vignette(caic)}}
\examples{

data(IsaacEtAl)
primatesBodySize <- macrocaic(species.rich ~ body.mass, data=IsaacEtAl.Primates)


chiropteraBodySize <- macrocaic(species.rich ~ body.mass, data=IsaacEtAl.Chiroptera, equal.branch.length=TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
