\name{vcv.array}
\alias{vcv.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a 2D or 3D variance-covariance matrix from a phylogeny}
\description{
The function turns a phylogeny into a variance-covariance matrix, as in the function \code{vcv.phylo} in the 'ape' package but can also return a 3D array retaining the individual branch lengths contributing to the shared branch lengths. This is useful for handling some branch length transformations, such as kappa, and has a lower overhead than repeatedly calling \code{vcv.phylo} on a phylogeny after transforming the vector of edge lengths.}
\usage{
vcv.array(phy, dim=2, compact=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{An object of class 'phylo'.}
  \item{dim}{Either 2, for a standard VCV matrix, or 3, for an array of branch lengths.}
  \item{compact}{A logical vector indicating the form to use for a 3D array.}
}
\details{
The compact form of the 3D array uses a shortened third dimension, which is only long enough to hold the maximum number of shared branches between root and tip for each pair of tips. Zeros are used to pad out this depth vector for tip pairs with shorter paths. The non-compact form returns 3D array showing, for each pair of tips and each node in the tree, either 0 if the node is not shared or the appropriate edge length if the node is shared. Note that, for maximally unbalanced trees, the size of the two forms will be identical.

The algorithm for the noncompact form is faster than for the compact form but it has very high memory overheads on big trees. The 2 dimensional algorithm is at least twice as fast as \code{vcv.phylo} on trees up to 2500 tips.

The \code{apply} function can be easily used to collapse the array down to a standard VCV matrix, as in the example. 
}
\value{
When dim = 2, a variance covariance matrix of class 'vcv.array' of dimension nTips by nTips with dimnames set from the tip labels of the phylogeny.

When dim = 3, a 3 dimensional array of class 'vcv.array' with dimensions of the number of taxa in the phylogeny for the rows and columns and either the maximum number of branches on the root to tip path or the number of internal nodes as the depth, depending on the setting of \code{compact}. The rows and columns are named using the tip labels from the phylogeny and the depth only named with node numbers if \code{compact} is TRUE. 
}
\author{David Orme}
\seealso{ \code{\link[ape]{vcv.phylo}},  \code{\link{pglm}}}
\examples{
tree <- rcoal(8)
tree.VCV <- vcv.phylo(tree)
tree.VCVA <- vcv.array(tree)

# reconstruct a simple VCV array
tree.VCVA.reduced <- apply(tree.VCVA, c(1,2), sum, na.rm=TRUE)

# minimal differences between the two
all((tree.VCVA.reduced - tree.VCV) < 1e-10)

# a kappa transformation of 0.5
apply(tree.VCVA ^ 0.5, c(1,2), sum, na.rm=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ utilities }% __ONLY ONE__ keyword per line
