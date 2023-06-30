#' Create a clade matrix from a phylogeny
#'
#' Takes a phylogeny in the 'ape' package format and converts it into a binary
#' matrix showing which tips (matrix columns) subtend from each node in the
#' phylogeny (matrix rows). This is a useful format for quickly calculating
#' branch length information for subsets of the phylogeny.
#'
#' The clade matrix shows the tips from a phylogeny that subtend from each
#' internal and external node. Each tip is represented as column showing the
#' nodes of which it is a member and hence each row shows the tips that are
#' members of a given node. Dropping columns gives a quick and easy way to find
#' out which edges are retained in a particular subset of the tree and this
#' structure is used for quickly calculating branch lengths calculations or
#' clade statistics.
#'
#' @param phy A object of class 'phylo'
#' @return A list of class 'clade.matrix' containing the following components:
#' \item{clade.matrix}{A binary m x n matrix, where m is the total number of
#' nodes in the phylogeny and n is the number of tips. An element is 1 if tip
#' $n_i$ subtends from a node $m_j$.} \item{edge.length}{A numeric vector of
#' length m showing the edge length leading to each node in the phylogeny and
#' named with the node number.} \item{tip.label}{A character vector of length n
#' giving the labels assigned to the tips of the phylogeny.} \item{edge}{The
#' edge matrix from the original phylogeny.}
#' @author David Orme
#' @keywords manip utilities
#' @examples
#'
#' data(perissodactyla)
#' clade.matrix(perissodactyla.tree)
#' @export
clade.matrix <- function(phy) {
    # OLD2NEW: CONVERTED

    # returns the phylogeny as a table showing clade
    # membership below each node.

    # check the object is a phylogeny
    if (!inherits(phy, "phylo")) stop("Phylogeny required")

    # get the number of tips and nodes
    nb.tips <- max(phy$edge) - phy$Nnode
    nb.nodes <- phy$Nnode
    nb.all <- nb.tips + nb.nodes

    # set-up the clade matrix
    mat.names <- list(edges = c(1:nb.all), tips = 1:nb.tips)
    clade.mat <- matrix(0, nrow = nb.all, ncol = nb.tips, dimnames = mat.names)

    # the diagonal from [1,1] are the tips
    diag(clade.mat) <- 1

    # now deal with internals
    node.members <- clade.members.list(phy)

    node.id <- names(node.members)
    for (rows in seq(along = node.id)) {
        clade.mat[node.id[rows], node.members[rows][[1]]] <- 1
    }

    RET <- list(
        clade.matrix = clade.mat,
        tip.label = phy$tip.label, edge = phy$edge
    )
    class(RET) <- "clade.matrix"

    # if they exist, get edge lengths into correct order, inserting root edge
    if (!is.null(phy$edge.length)) {
        if (is.null(phy$root.edge)) {
            edge.len <- c(phy$edge.length, 0)
        } else {
            edge.len <- c(phy$edge.length, phy$root.edge)
        }

        names(edge.len) <- c(phy$edge[, 2], nb.tips + 1)
        edge.len <- edge.len[as.character(mat.names$edges)]
        RET$edge.length <- edge.len
    }

    return(RET)
}
