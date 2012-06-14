\name{fusco.test}
\alias{fusco.test}
\alias{summary.fusco}
\alias{print.fusco}
\alias{plot.fusco}


\title{Imbalance statistics using Fusco and Cronk's method.} 
\description{
Fusco and Cronk (1995) described a method for testing the imbalance of phylogenetic trees based on looking at the distribution of I. I is calculated using the number of tips descending from each side of a bifurcating node using the formula I = (B-m)/(M-m) and is bounded between 0 (a perfectly balanced node) and 1 (maximum imbalance). B is the larger number of tips descending from each branch, M is the maximum size of this larger group (i.e. a 1 : (S-1) split, where S is the total number of descendent tips), and m is the minimum size of the larger group (ceiling of S/2). The method can cope with small proportions of polytomies in the phylogeny and these are not used in calculating balance statistics. It can also incorporate information about species richness at the tips of the phylogeny and can therefore be used to distinguish between an unbalanced topology and the unbalanced distribution of diversity at the tips of a phylogeny. 

Purvis et al. (2002) demonstrated that I is not independent of the node size S, resulting in a bias to the expected median of 0.5. They proposed a modification (I') that corrects this to give a statistic with an expected median of 0.5 regardless of node size. The defaults in this function perform testing of imbalance using I', but it is also possible to use the original measure proposed by Fusco and Cronk (1995).
}
\usage{
fusco.test(phy, data , names.col, rich, tipsAsSpecies=FALSE, randomise.Iprime=TRUE,  reps=1000, conf.int=0.95)
\method{print}{fusco}(x, ...)
\method{summary}{fusco}(object, ...)
\method{plot}{fusco}(x, correction=TRUE, nBins=10, right=FALSE, I.prime=TRUE, plot=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{An object of class 'comparative.data' or of class 'phylo'.}
  \item{data}{A data frame containing species richness values. Not required if phy is a 'comparative.data' object.}
  \item{names.col}{A variable in \code{data} identifying tip labels. Not required if phy is a 'comparative.data' object.}
  \item{rich}{A variable  identifying species richness.}
  \item{tipsAsSpecies}{A logical value. If TRUE, a species richness column need not be specified and the tips will be treated as species.}
  \item{randomise.Iprime}{Use a randomization test on calculated I' to generate confidence intervals.}
  \item{reps}{Number of replicates to use in simulation or randomization.}
  \item{conf.int}{Width of confidence intervals required.}
  \item{x, object}{An object of class 'fusco'.}
  \item{correction}{Apply the correction described in Appendix A of Fusco and Cronk (1995) to the histogram of nodal imbalance.}
  \item{nBins}{The number of bins to be used in the histogram of nodal imbalance.}
  \item{right}{Use right or left open intervals in plotting the distribution and calculating the correction}
  \item{I.prime}{Plot distribution of I' or I.}
  \item{plot}{If changed to FALSE, then the plot method does not plot the frequency histogram. Because the method invisibly returns a table of histogram bins along with the observed and corrected frequencies, this isn't as dim an option as it sounds.}
  \item{...}{Further arguments to generic methods}
}
\details{
I is calculated only at bifurcating nodes giving rise to more than 3 tips (or more than 3 species at the tips): nodes with three or fewer descendants have no variation in I and are not informative in assessing imbalance. The expected distribution of the nodal imbalance values between 0 and 1 is theoretically uniform under a Markov null model. However, the range of possible I values at a node is constrained by the number of descendent species. For example, for a node with 8 species, only the values 0, 0.33, 0.66, 1 are possible, corresponding to 4:4, 5:3, 6:2 and 7:1 splits (Fusco and Cronk, 1995). As node size increases, this departure from a uniform distribution decreases. The plot method incorporates a correction, described by Fusco and Cronk (1995), that uses the distribution of all possible splits at each node to characterize and correct for the departure from uniformity.

The randomization option generates confidence intervals around the mean I'.
}
\value{
The function \code{fusco.test} produces an object of class 'fusco' containing:
	\item{observed}{A data frame of informative nodes showing nodal imbalance statistics. If the phylogeny has labelled nodes, then the node names are also returned in this data frame.}
	\item{median}{The median value of I.}
	\item{qd}{The quartile deviation of I.}
	\item{tipsAsSpecies}{A logical indicating whether the tips of the trees were treated as species or higher taxa.}
	\item{nInformative}{The number of informative nodes.}
	\item{nSpecies}{The number of species distributed across the tips.}
	\item{nTips}{The number of tips.}
	\item{reps}{The number of replicates used in randomization.}
	\item{conf.int}{The confidence levels used in randomization.}
If \code{randomise.Iprime} is TRUE, or the user calls \code{fusco.randomize} on a 'fusco' object, then the following are also present.
	\item{randomised}{A data frame of mean I' from the randomized observed values.}
	\item{rand.mean}{A vector of length 2 giving confidence intervals in mean I'.}
}
\references{
Fusco, G. & Cronk, Q.C.B. (1995) A New Method for Evaluating the Shape of Large Phylogenies. J. theor. Biol. 175, 235-243

Purvis A.,  Katzourakis A. & Agapow, P-M (2002) Evaluating Phylogenetic Tree Shape: Two Modifications to Fusco & Cronk's Method. J. theor. Biol. 214, 93-103.
}
\author{David Orme, Andy Purvis}
\examples{
data(syrphidae)
syrphidae <- comparative.data(phy=syrphidaeTree, dat=syrphidaeRich, names.col=genus)
summary(fusco.test(syrphidae, rich=nSpp))
summary(fusco.test(syrphidae, tipsAsSpecies=TRUE))
plot(fusco.test(syrphidae, rich=nSpp))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{htest}