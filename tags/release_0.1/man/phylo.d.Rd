\name{phylo.d}
\alias{phylo.d}
\alias{print.phylo.d}
\alias{summary.phylo.d}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates the phylogenetic D statistic}
\description{
Calculates the D value, a measure of phylogenetic signal in a binary trait, and tests the estimated D value for significant departure from both random association and the clumping expected under a Brownian evolution threshold model.
}
\usage{
phylo.d(data, phy, names.col, binvar, permut = 1000)
\method{print}{phylo.d}(x, ...)
\method{summary}{phylo.d}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A 'comparative.data' or 'data.frame' object.}
  \item{phy}{An object of class 'phylo', required when data is not a 'comparative.data' object.}
  \item{names.col}{A name specifying the column in 'data' that matches rows to tips in 'phy', required when data is not a 'comparative.data' object.}
  \item{binvar}{The name of the variable in \code{data} holding the binary variable of interest.}
  \item{permut}{Number of permutations to be used in the randomisation test.}
  \item{x}{An object of class 'phylo.d'}
  \item{object}{An object of class 'phylo.d'}
  \item{\ldots}{Further arguments to print and summary methods}
}

\details{
The sum of changes in estimated nodal values of a binary trait along edges in a phylogeny (D) provides a measure of the phylogenetic signal in that trait (Fritz and Purvis, 2010). If a trait is highly conserved, with only a basal division between two clades expressing either trait value, then the only change will be along the two daughters at the root. This will give a summed value of 1: the two differences between the root nodal value of 0.5 and the ancestors of the 1 and 0 clades. In contrast, if the trait is labile, more differences will be observed and the sum will be higher.

This function calculates the observed D for a binary trait on a tree and compares this to the value of D found using an equal number of simulations under each of two models:
\describe{
\item{Phylogenetic randomness}{Trait values are randomly shuffled relative to the tips of the phylogeny and D is calculated.}
\item{Brownian threshold model}{A continuous trait is evolved along the phylogeny under a Brownian process and then converted to a binary trait using a threshold that reproduces the relative prevalence of the observed trait.}
}
The value of D depends on phylogeny size - more sister clades yield higher sums - and so the means of the  two sets of simulated data are used as calibrations to scale both observed and simulated values of D to set points of 0 (as phylogenetically conserved as expected under a Brownian threshold model) and 1 (random). The value of D can be both smaller than 0 (highly conserved) and greater than 1 (overdispersed) and the distributions of scaled D from the simulations are used to assess the significance of the observed scaled D.
}
\value{
Returns an object of class 'phylo.d', which is a list of the following:
\item{DEstimate}{The estimated D value}
\item{Pval1}{A p value, giving the result of testing whether D is significantly different from one}
\item{Pval0}{A p value, giving the result of testing whether D is significantly different from zero}
\item{Parameters}{A list of the Observed, MeanRandom and MeanBrownian sums of sister-clade differences}
\item{Permutations}{A list with elements random and brownian, containing the sums of sister-clade differences from random permutations and simulations of Brownian evolution under a threshold model}
\item{NodalVals}{A list with the elements observed, random and brownian, containing the nodal values estimated for the observed trait and permutations. The values are as matrices with rows labelled by the node names in the comparative data object.}
\item{binvar}{The binary variable used}
\item{phyName}{The name of the phylogeny object used}
\item{dsName}{The name of the dataframe used}
\item{nPermut}{The number of permutations used}
}
\references{Fritz, S. A. and Purvis, A. (2010). Selectivity in mammalian extinction risk and threat types: a new measure of phylogenetic signal strength in binary traits. Conservation Biology, 24(4):1042-1051.}
\author{Susanne Fritz (SFritz@bio.ku.dk) and David Orme}
\examples{
data(BritishBirds)
BritishBirds <- comparative.data(BritishBirds.tree, BritishBirds.data, binomial)
redPhyloD <- phylo.d(BritishBirds, binvar=Red_list)
print(redPhyloD)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ htest }% __ONLY ONE__ keyword per line

%% ## OLD example with plots - use real example
%% data(shorebird)
%% 
%% # add threatened status (threat data from IUCN Red List 2008)
%% # and a phylogenetically clumped dummy example to the data frame
%% shorebird.data$threat <- rep(c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0), times=c(22,2,2,2,1,1,4,1,8,1,3,1,1,2,20))
%% shorebird.data$clustered <- rep(c(1,0,1,0,1,0,1,0,1,0,1,0), times=c(4,13,2,12,1,2,1,7,2,2,2,23))
%% 
%% # comparative data object
%% shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)
%% 
%% # threat: D is significantly different a brownian threshold model
%% ex1 <- phylo.d(shorebird, binvar=threat, permut=100)
%% # clustered: D is significantly different from a random model
%% ex2 <- phylo.d(shorebird, binvar=clustered, permut=100)
%% 
%% # plot the data
%% par(mfrow=c(1,2))
%% tipOrd <- match(shorebird.data$Species, shorebird.tree$tip.label)
%% plot(shorebird.tree, show.tip.label=FALSE, main='Random', no.margin=TRUE) 
%% tiplabels(tip=tipOrd, pch=21, bg=ifelse(shorebird.data$threat==0, 'white','black'))
%% plot(shorebird.tree, show.tip.label=FALSE, main='Clustered', no.margin=TRUE) 
%% tiplabels(tip=tipOrd, pch=21, bg=ifelse(shorebird.data$clustered==0, 'white','black'))
%% 
