\name{crunch}
\alias{crunch}
\alias{contrCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparative analysis using the crunch algorithm.}
\description{
Calculate a linear model using the crunch algorithm.
}
\usage{
crunch(formula, data, phy, names.col, stand.contr = TRUE, ref.var=NULL, node.depth=NULL, polytomy.brlen=0, equal.branch.length=FALSE, factor.action="abort")
}


\arguments{
  \item{formula}{A model formula.}
  \item{data}{An 'comparative.data' object. Alternatively, a data frame.}
  \item{phy}{An object of class 'phylo', required when data is not a 'comparative.data' object.}
  \item{names.col}{A name specifying the column in 'data' that matches rows to tips in 'phy', required when data is not a 'comparative.data' object.}
  \item{stand.contr}{A logical flag indicating whether or not to standardize contrasts}.
  \item{ref.var}{A reference variable present in the model that is used to specify the direction of calculation of contrasts. If null, this is assumed to be the first explanatory variable.}
  \item{node.depth}{A positive integer greater than 1 used to restrict the model to contrasts with a node depth less than or equal to the specified depth. Tips have a depth of 1.}
  \item{polytomy.brlen}{The internal branch length used for calculating contrasts at a polytomy, following Pagel's (1992) method. }
  \item{equal.branch.length}{If set to 'TRUE' then all branch lengths are set to 2.}
 \item{factor.action}{One of "abort", "warn" or "allow", describing whether to stop if the formula contains a factor ("abort"), or continue after converting the factor to a numeric variable, either with ("warn") or without ("allow") a warning.}
}

\details{
This function implements the 'crunch' algorithm for modelling the relationship between variables that are phylogenetically non-independent. The method was first described by Felsenstein (1985) and subsequently extended to permit the use of phylogenies with polytomies by Pagel (1992). This method was previously implemented in the Mac Classic computer programs CAIC, written by Andy Purvis, Andy Rambaut (Purvis and Rambaut, 1995) and updated by Nick Isaac and Paul-Michael Agapow.

The user provides a comparative dataset. The formula specifies the model to be fitted and contrasts are calculated in those variables. The specified reference variable is used to ensure that contrasts for multivariate models are calculated in a consistent direction at each node. The function \code{crunch()} acts as a data preparation wrapper for the function \code{contrCalc()}, which is not intended to be directly called by users. Missing data can be present in the explanatory variables: the algorithm makes use of the complete data available at each node as was the case with CAIC.

The resulting table of contrasts is then used to fit the specified model - note that the intercept is automatically dropped from the model if present, following REF HERE. 

Contrasts at polytomies are calculated following Pagel (1992). The descendants from the node are split into two groups based on whether they are above or below the group mean in the reference variable. If there is no variation in the reference variable, then a 1:(N-1) split is used. Weighted means in the variables are then calculated for each subgroup and a contrast is calculated between these values using an arbitrary internal branch length.
}

\section{Warning}{
At a polytomy, subtracting the internal branch length from the real branch lengths can lead to negative branch lengths. CAIC used a hard-coded internal branch length of 1 for calculating crunch contrasts at polytomies. From version 2.6.9, CAIC issued a warning if this lead to negative branch lengths but allowed contrast calculation to continue. In contrast, the implementation in \code{crunch()} uses a default internal branch length (\code{polytomy.brlen}) of 0 and also treats a negative branch length in a polytomy calculation as an error. In either case, contrast calculation on negative branch lengths is not a desirable outcome. Duplication of CAIC results therefore requires \code{polytomy.brlen} to be set to 1 and an analyis \emph{cannot} be duplicated precisely if the phylogeny contains polytomies with descending branches shorter than 1. The method used by \code{pic.lm} to handle polytomies avoids such problems.
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

data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)
crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data=shorebird)
summary(crunchMod)

# plot the contrasts
crunchTab <- caic.table(crunchMod)
plot(Egg.Mass ~ F.Mass, crunchTab)

# for the actual model diagnostics
par(mfrow=c(3,2))
caic.diagnostics(crunchMod)
}

%% # old CAIC style arguments - don't publicize
%% crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, phy=shorebird.tree,  data=shorebird.data, names.col=Species)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
