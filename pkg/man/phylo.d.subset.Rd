\name{phylo.d.subset}
\alias{phylo.d.subset}
\alias{print.phylo.d.subset}
\alias{summary.phylo.d.subset}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates the phylogenetic D statistic across clades within a phylogeny}
\description{
Calculates the D value, a measure of phylogenetic signal in a binary trait, and tests the estimated D value for significant departure from both random association and the clumping expected under a Brownian evolution threshold model. Does this across clades within a phylogeny.
}
\usage{
phylo.d.subset(data, phy, names.col, binvar, permut = 1000, rnd.bias=NULL, 
	           min.tips=1, max.tips=length(data$phy$tip.label), min.nodes=1, 
			   max.nodes=data$phy$Nnode, verbose=FALSE)
\method{print}{phylo.d.subset}(x, ...)
\method{summary}{phylo.d.subset}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A 'comparative.data' or 'data.frame' object.}
  \item{phy}{An object of class 'phylo', required when data is not a 'comparative.data' object.}
  \item{names.col}{A name specifying the column in 'data' that matches rows to tips in 'phy', required when data is not a 'comparative.data' object.}
  \item{binvar}{The name of the variable in \code{data} holding the binary variable of interest.}
  \item{permut}{Number of permutations to be used in the randomisation test.}
  \item{rnd.bias}{An optional name of a variable in \code{data} holding probability weights to bias the generation of the random distribution. See 'destails'}
  \item{verbose}{Logical; do you want to know how many clades are being assessed, and see when each is being assessed?}
  \item{min.tips}{The minimum number of tips a clade should have for it to have a D value calculated. Defaults to 1 (i.e. no limit).}
  \item{max.tips}{The maximum number of species a clade should have for it to have a D value calculated. Defaults to the number of species in the whole phylogeny (i.e. no limit).}
  \item{min.nodes}{The minimum number of nodes a clade should have for it to have a D value calculated. Defaults to 1 (i.e. no limit).}
  \item{max.nodes}{The maximum number of nodes a clade should have for it to have a D value calculated. Defaults to the number of nodes in the whole phylogeny (i.e. no limit).}
  \item{x}{An object of class 'phylo.d.subset'}
  \item{object}{An object of class 'phylo.d.subset'}
  \item{\ldots}{Further arguments to print and summary methods}}

\details{
A wrapper function for \code{\link{phylo.d}}, calculating D values for clades within a given dataset. These clades can be filtered according to the number of species and nodes using the arguments above. See \code{\link{phylo.d}} for more details on the method itself.

Any clades for which there is no variation in the binary variable have \code{NA} values for all of the below slots.}
\value{
Returns an object of class 'phylo.d.subset', which is a list of the following:
\item{raw}{A list of the raw output from \code{\link{phylo.d}} for each clade}
\item{DEstimate}{A vector of the estimated D values}
\item{Pval1}{A vector of p values, giving the result of testing whether D is significantly different from one, for each clade}
\item{Pval0}{A vector of p values, giving the result of testing whether D is significantly different from zero, for each clade}
\item{phy.depth}{A numeric vector giving the age of the clade for which each value was calculated}}
\references{Fritz, S. A. and Purvis, A. (2010). Selectivity in mammalian extinction risk and threat types: a new measure of phylogenetic signal strength in binary traits. Conservation Biology, 24(4):1042-1051.}
\author{Susanne Fritz (SFritz@bio.ku.dk), Will Pearse and  David Orme}
\examples{
data(BritishBirds)
BritishBirds <- comparative.data(BritishBirds.tree, BritishBirds.data, binomial)
# Look at big clades only
\dontrun{
bigClades <- phylo.d.subset(BritishBirds, binvar=Red_list, verbose=TRUE, min.tips=10, min.nodes=5)
print(bigClades)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ htest }% __ONLY ONE__ keyword per line
