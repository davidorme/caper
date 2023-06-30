#' Identify tips descended from a node
#'
#' Obtains a vector of the tips subtending from either one node or all nodes in
#' a phylogeny.
#'
#' The function \code{clade.members.list} runs \code{clade.members} over each
#' node in the phylogeny, possibly including the external nodes as indicated by
#' the tips argument, and returns a list of vectors showing the members of the
#' clade defined by each node.
#'
#' @aliases clade.members clade.members.list
#' @param x A integer identifying the node for which a list of tips is
#' required.
#' @param phy An object of class 'phylo'.
#' @param tips A logical indicating whether to include external node membership
#' in the list.
#' @param tip.labels A logical flag indicating whether to return the node
#' numbers of the tips or their tip labels.
#' @param include.nodes A logical flag indicating whether to return the node
#' number of descendent internal nodes
#' @return A numeric vector of external node (i.e. tip) numbers or a character
#' vector of tip labels for a single internal node or, for
#' \code{clade.members.list}, a list of such vector for all nodes in the
#' phylogeny. If \code{include.nodes} is \code{TRUE} then \code{clade.members}
#' returns a list of length two containing a vector of the descendent tips and
#' a vector of the descendent internal nodes - \code{clade.members.list} then
#' contains a list of such lists.
#' @author David Orme, Lynsey McInnes
#' @keywords manip utilities
#' @examples
#'
#' data(perissodactyla)
#' # use comparative.data to add node labels
#' perisso <- comparative.data(perissodactyla.tree, perissodactyla.data, Binomial, na.omit = FALSE)
#' plot(perisso$phy, show.node.label = TRUE)
#' clade.members(22, perisso$phy, tip.labels = TRUE)
#' clade.members.list(perisso$phy, tip.labels = FALSE)
#'
#' @export
clade.members <- function(x, phy, tip.labels = FALSE,
                          include.nodes = FALSE) {
    # NEW2OLD: CONVERTED...

    # returns a vector of the tips that descend from an identified node
    if (!inherits(phy, "phylo")) stop("Phylogeny required")

    NallNodes <- max(phy$edge)
    Ntips <- max(phy$edge) - phy$Nnode

    if (!(x %in% 1:NallNodes)) stop("Node not in range for phylogeny")

    # find the children of the node, append them to the vector of nodes (x)
    # and remove the parent, until all the nodes in the vector are tips...
    # now updated to keep track of parents...

    intN <- x[x > Ntips]
    descNode <- numeric(length = 0)

    while (length(intN) > 0) {
        minIntN <- min(intN)
        childOfMinIntN <- with(phy, edge[, 2][which(edge[, 1] == minIntN)])

        descNode <- c(descNode, minIntN)
        x <- c(x[x != minIntN], childOfMinIntN)

        intN <- x[x > Ntips]
    }

    RET <- unique(x)

    if (tip.labels) {
        RET <- phy$tip.label[x]
    }

    if (include.nodes) {
        RET <- list(tips = RET, nodes = descNode)
    }

    return(RET)
}

clade.members.list <- function(phy, tips = FALSE, tip.labels = FALSE,
                               include.nodes = FALSE) {
    # OLD2NEW CONVERTED

    # returns a list of vectors showing the tips
    # subtending from each node in the tree
    if (!inherits(phy, "phylo")) stop("Phylogeny required")

    nodes <- 1:max(phy$edge)

    if (!tips) nodes <- nodes[nodes > length(nodes) - phy$Nnode]

    clade.list <- mapply(
        clade.members, nodes,
        MoreArgs = list(
            phy = phy, tip.labels = tip.labels,
            include.nodes = include.nodes
        ), SIMPLIFY = FALSE
    )
    names(clade.list) <- nodes

    return(clade.list)
}
