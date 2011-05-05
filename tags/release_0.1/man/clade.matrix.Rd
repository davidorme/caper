\name{clade.matrix}
\alias{clade.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a clade matrix from a phylogeny}
\description{
Takes a phylogeny in the 'ape' package format and converts it into a binary matrix showing which tips (matrix columns) subtend from each node in the phylogeny. This is a useful format for quickly calculating branch length information for subsets of the phylogeny.}
\usage{
clade.matrix(phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{A object of class 'phylo'}
}
\details{
The clade matrix shows the tips from a phylogeny that subtend from each internal and external node. Each tip is represented as column showing the nodes of which it is a member and hence each row shows the tips that are members of a given node. Dropping columns gives a quick and easy way to find out which edges are retained in a particular subset of the tree and this structure is used for quickly calculating branch lengths calculations or clade statistics.
}
\value{
  A list of class 'clade.matrix' containing the following components:
  \item{clade.matrix}{A binary m x n matrix, where m is the total number of nodes in the phylogeny and n is the number of tips. An element is 1 if tip $n_i$ subtends from a node $m_j$.}
  \item{edge.length}{A numeric vector of length m showing the edge length leading to each node in the phylogeny and named with the node number.}
  \item{tip.label}{A character vector of length n giving the labels assigned to the tips of the phylogeny.}
}
\author{David Orme}

%\seealso{}%%\code{\link{pd.calc}}, \code{\link{pd.bootstrap}}, \code{\link{clade.stats}}}
\examples{
data(perrisodactyla)
clade.matrix(perrisodactyla.tree)
}

\keyword{manip}
\keyword{utilities}
