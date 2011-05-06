\name{pd.calc}
\alias{pd.calc}
\alias{pd.bootstrap}
\alias{ed.calc}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate and bootstrap phylogenetic diversity measurements.}
\description{
These functions calculate various phylogenetic diversity measures for either a given set of nodes on a tree or for a randomly chosen set of nodes of a given size. The ed.calc function calculates a related species-level measurement of evolutionary distinctness.}
\usage{
pd.calc(cm, tip.subset = NULL, method = "TBL", root.edge=FALSE)
pd.bootstrap(cm, ntips, reps = 1000, method = "TBL", tip.weights = NULL)
edge.calc(cm, polytomy.cf=c("isaac","mooers","none"))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cm}{A object of class 'clade matrix'. Alternatively an object of class 'phylo', which will be converted to a clade.matrix.}
  \item{tip.subset}{An optional vector identifying the subset of tips to use for PD calculations. If no tip.subset is provided the method is applied to the complete phylogeny [Hmm.. this might be undesirable]. Can either be a character vector, in which case the elements are matched against tip labels, or a vector of positive integers in the range 1 to the number of tips, in which case the tips with those numbers are used.}
  \item{method}{ One of 'TBL', 'MST', 'UEH', 'SBL', defaulting to 'TBL'. See details.}
  \item{root.edge}{Logical indicating whether to include the root edge length in calculations, defaulting to FALSE.}
  \item{ntips}{A single integer giving the number of tips to be selected.}
  \item{reps}{The number of replicate values to calculate.}
  \item{tip.weights}{A numeric vector containing weights for all the tips in the phylogeny. Each element must be named in order to match weights to the tips.}
  \item{polytomy.cf}{Which correction factor to use for calculating ED at polytomies. One of 'isaac', 'mooers' or 'none'.}
}
\details{
There are five implemented PD measures:
\define{
  \item{Total Branch Length (TBL)}{The sum of all the edge lengths in the subtree given by the tip subset. This measure can be partitioned into the two next measures.}
  \item{Shared Branch Length (SBL)}{The sum of all edges in the subtree that are shared by more than one tip.}
  \item{Unique Evolutionary History (UEH)}{The sum of the edge lengths that give rise to only one tip in the subtree.}
  \item{Length of tip branch lengths (TIPS)}{Length of tip branch lengths (TIPS)}{Unlike UEH, this measure does not use the unique paths to each tips on the \strong{subtree} and instead gives the sum of the unique branches leading to the tips on the \strong{complete tree}.}
  \item{Minimum Spanning Tree (MST)}{The sum of the lengths of the edges for the smallest tree that links the subset tips, excluding any edges below the node of the most recent common ancestor.}}

These options are illustrated in the caper package vignette. The pd.calc function returns the PD value for a given set of tips, whereas the pd.bootstrap function returns a vector of PD values for randomly selected sets of tips of a given size.

The ed.calc function returns the evolutionary distinctness (ED) metric (Isaac et al, 2007) for the tips of a given phylogeny. The idea behind the ED measure is that the evolutionary history of each branch is shared equally between all tips descending from that branch. Each branch therefore has a per-tip values of the branch length divided by the number of descendants and the ED value for a tip is the sum of those per-tip contributions over the path to the root of the phylogeny. Polytomies inflate apparent ED since the branches of a properly resolved polytomy must be shorter than the branch lengths on the unresolved polytomy. The function provides two correction factors for this: 'isaac' uses a correction factor calibrated from simulations and 'mooers' uses empirical predictions from a pure birth model.
}
\value{
Both pd.calc and pd.bootstrap return a vector containing either a single value for the phylogenetic diversity of a given set of tips or a vector of length 'nrep' containing the pd values for a random set of tips of a given size. The method used is stored in the 'pd.method' attribute of the vector.

The ed.calc function returns a list containing:
\describe{
  \item{branch}{A data frame of the ED contributions arising from each branch.}
  \item{spp}{A data frame of the summed ED contributions for each species.}
}
}
\references{Faith, DP, 
Isaac, N. J. B., Turvey, S. T., Collen, B., Waterman, C., and Baillie, J. E. M. (2007). Mammals on the edge: Conservation priorities based on threat and phylogeny. Plos One, 2(3):e296}
\author{David Orme, Gavin Thomas, Nick Isaac}
\examples{
tre <- read.tree(text="((((A:1,B:1):1.5,C:2.5):0.5,(D:0.6,E:0.6):2.4):0.5,((F:1.9,G:1.9):0.8,(H:1.6,I:1.6):1.1):0.8):0.2;")
clmat <- clade.matrix(tre)
tips <- c("A","C","D","E","G","H")
pd.calc(clmat, tip.subset=tips)
pd.calc(clmat, tip.subset=c(1,3,4,5,7,8))
pd.calc(clmat, tip.subset=tips, root.edge=TRUE)

pd.bootstrap(clmat, 6, 'TBL')

data(IsaacEtAl)
primatesCM <- clade.matrix(primates.tree)
primatesED <- ed.calc(primatesCM)
}
\keyword{utilities}% at least one, from doc/KEYWORDS
