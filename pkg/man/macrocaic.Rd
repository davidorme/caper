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
  \item{macroMethod}{One of either "RRD" or "PDI" (see Details).}
  \item{stand.contr}{A logical flag indicating whether to standardize the contrasts.}
  \item{ref.var}{Identifies a predictor variable used for determining the direction of contrasts.}
  \item{node.depth}{A positive integer greater than 1 used to restrict the model to contrasts with a node depth less than or equal to the specified depth. Tips have a depth of 1.}
  \item{macroMinSize}{A positive integer giving the minimum species richness at a node for contrasts to be included in the model.}
  \item{equal.branch.length}{If set to 'TRUE' then all branch lengths are set to 2.}
}
\details{
The 'macrocaic' function fits a regression to the formula provided using 'crunch' contrasts for continuous explanatory variables and species richness contrasts for the response. The species richness contrasts are either the relative rate difference (RRD) or proportion dominance index (PDI):

\deqn{RRD = \ln\left(\frac{N_1}{N_2}\right)}
     {RRD = ln(N_1/N_2)}

\deqn{PDI = \left(\frac{N_1}{N_1+N_2}\right)-0.5}
     {PDI = (N_1/(N_1 + N_2))-0.5}

The values \eqn{N_1} and \eqn{N_2} are the species richness of the two daughter nodes and \eqn{N_1} is the species richness of the clade with the larger value of the reference variable. Species richness contrasts are not calculated at polytomies. Nodal values for species richness are calculated as the sum of the richness of the daughter nodes.
}
\value{
  A object of class 'caic'.
}
\references{
Felsenstein, J.  (1985).  Phylogenies and the comparative method.  Am. Nat.  125, 1-15
Agapow, P.-M. and Isaac, N. J. B. (2002) MacroCAIC: correlates of species richness. Diversity & Distributions, 8, 41-43
Isaac, N., Agapow, P., Harvey, P., and Purvis, A. (2003). Phylogenetically nested com- parisons for testing correlates of species richness: A simulation study of continuous variables. Evolution, 57(1):18–26.
}
\author{David Orme}
\seealso{\code{\link{caic-class}} for 'caic' object structure and methods.}
\examples{
data(IsaacEtAl)
primates <- comparative.data(primates.tree, primates.data, binomial, na.omit=FALSE)
primatesBodySize <- macrocaic(species.rich ~ body.mass, data=primates)
summary(primatesBodySize)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
