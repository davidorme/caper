\name{clade.members}
\alias{clade.members}
\alias{clade.members.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify tips descended from a node}
\description{
Obtains a vector of the tips subtending from either one node or all nodes in a phylogeny.}
\usage{
clade.members(x, phy, tip.labels = FALSE, include.nodes=FALSE)
clade.members.list(phy, tips = FALSE, tip.labels = FALSE, include.nodes=FALSE)
}

\arguments{
  \item{x}{A integer identifying the node for which a list of tips is required.}
  \item{phy}{An object of class 'phylo'. }
  \item{tips}{ A logical indicating whether to include external node membership in the list.}
  \item{tip.labels}{A logical flag indicating whether to return the node numbers of the tips or their tip labels.}
  \item{include.nodes}{A logical flag indicating whether to return the node number of descendent internal nodes}
}
\details{
The function \code{clade.members.list} runs \code{clade.members} over each node in the phylogeny, possibly including the external nodes as indicated by the tips argument, and returns a list of vectors showing the members of the clade defined by each node.
}

\value{
A numeric vector of external node (i.e. tip) numbers or a character vector of tip labels for a single internal node or, for \code{clade.members.list}, a list of such vector for all nodes in the phylogeny. If \code{include.nodes} is \code{TRUE} then \code{clade.members} returns a list of length two containing a vector of the descendent tips and a vector of the descendent internal nodes -  \code{clade.members.list} then contains a list of such lists.
}
\author{David Orme, Lynsey McInnes}


%\seealso{} %% \code{\link{clade.stats}} }
\examples{
data(perissodactyla)
# use comparative.data to add node labels
perisso <- comparative.data(perissodactyla.tree, perissodactyla.data, Binomial, na.omit=FALSE)
plot(perisso$phy, show.node.label=TRUE)
clade.members(22, perisso$phy, tip.labels=TRUE)
clade.members.list(perisso$phy, tip.labels=FALSE)

}
\keyword{ manip }% at least one, from doc/KEYWORDS
\keyword{ utilities}% __ONLY ONE__ keyword per line
