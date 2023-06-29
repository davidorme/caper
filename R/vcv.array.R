#' Create a 2D or 3D variance-covariance matrix from a phylogeny
#'
#' The function turns a phylogeny into a variance-covariance matrix, as in the
#' function \code{vcv.phylo} in the 'ape' package but can also return a 3D
#' array retaining the individual branch lengths contributing to the shared
#' branch lengths. This is useful for handling some branch length
#' transformations, such as kappa, and has a lower overhead than repeatedly
#' calling \code{vcv.phylo} on a phylogeny after transforming the vector of
#' edge lengths.
#'
#' The compact form of the 3D array uses a shortened third dimension, which is
#' only long enough to hold the maximum number of shared branches between root
#' and tip for each pair of tips. Zeros are used to pad out this depth vector
#' for tip pairs with shorter paths. The non-compact form returns 3D array
#' showing, for each pair of tips and each node in the tree, either 0 if the
#' node is not shared or the appropriate edge length if the node is shared.
#' Note that, for maximally unbalanced trees, the size of the two forms will be
#' identical.
#'
#' The algorithm for the noncompact form is faster than for the compact form
#' but it has very high memory overheads on big trees. The 2 dimensional
#' algorithm is at least twice as fast as \code{vcv.phylo} on trees up to 2500
#' tips.
#'
#' The \code{apply} function can be easily used to collapse the array down to a
#' standard VCV matrix, as in the example.
#'
#' @param phy An object of class 'phylo'.
#' @param dim Either 2, for a standard VCV matrix, or 3, for an array of branch
#' lengths.
#' @param compact A logical vector indicating the form to use for a 3D array.
#' @return When dim = 2, a variance covariance matrix of class 'VCV.array' of
#' dimension nTips by nTips with dimnames set from the tip labels of the
#' phylogeny.
#'
#' When dim = 3, a 3 dimensional array of class 'VCV.array' with dimensions of
#' the number of taxa in the phylogeny for the rows and columns and either the
#' maximum number of branches on the root to tip path or the number of internal
#' nodes as the depth, depending on the setting of \code{compact}. The rows and
#' columns are named using the tip labels from the phylogeny and the depth only
#' named with node numbers if \code{compact} is TRUE.
#' @author David Orme
#' @seealso \code{\link[ape]{vcv.phylo}}, \code{\link{pgls}}
#' @keywords manip utilities
#' @examples
#'
#' tree <- rcoal(8)
#' tree.VCV <- vcv.phylo(tree)
#' tree.VCVA <- VCV.array(tree)
#'
#' # reconstruct a simple VCV array
#' tree.VCVA.reduced <- apply(tree.VCVA, c(1, 2), sum, na.rm = TRUE)
#'
#' # minimal differences between the two
#' all((tree.VCVA.reduced - tree.VCV) < 1e-10)
#'
#' # a kappa transformation of 0.5
#' apply(tree.VCVA^0.5, c(1, 2), sum, na.rm = TRUE)
#'
VCV.array <- function(phy, dim = 2, compact = TRUE) {
    ## turns a phylogeny into a 3d array similar to a VCV matrix
    ## but keeping each beanch length separate. This is useful for
    ## handling branch length transformations in functions where
    ## VCVs are used to handle phylogenetic structure

    ## rewritten to use new ape, via the clade matrix structure

    if (!inherits(phy, "phylo")) {
        stop("object \"phy\" is not of class \"phylo\"")
    }

    if (is.null(phy$edge.length)) {
        stop(
            "Object \"phy\" is missing edge lengths. ",
            "Cannot provide VCV matrix."
        )
    }


    if (!dim %in% 2:3) {
        stop("dim must be 2 or 3, for a VCV matrix or array respectively. ")
    }


    cm <- clade.matrix(phy)
    cmM <- cm$clade.matrix
    cmE <- cm$edge.length
    cmEM <- cmM * cmE

    if (dim == 2) {
        V <- crossprod(cmM, cmEM)
        dimnames(V) <- list(phy$tip.label, phy$tip.label)
    } else {
        if (compact) {
            nTip <- dim(cmM)[2]
            max.node.depth <- max(colSums(cmM))

            V <- array(0,
                dim = c(nTip, nTip, max.node.depth),
                dimnames = list(phy$tip.label, phy$tip.label, NULL)
            )

            for (i in 1:nTip) {
                Vslice <- cmEM * cmM[, i]
                Vind <- which(Vslice > 0, arr.ind = TRUE)
                Vval <- Vslice[Vind]
                Vrle <- rle(Vind[, 2])
                Vind[, 1] <- unlist(mapply(seq, from = 1, Vrle$lengths))
                V[i, , ][Vind[, c(2, 1)]] <- Vval
            }
        } else {
            # returns a big 3d array showing, for each pair of tips, either 0
            # (not shared) or the appropriate edge length if the node is shared
            # - not good on big trees! but it is faster than the previous
            # version

            # multiply each column of the edge length matrix by the clade matrix
            V <- apply(cmEM, 2, function(X) X * cmM)
            dims <- dim(cmM)

            # gives a (Nnodes by Ntips) by Ntips matrix, which needs reshaping
            # into an array of Nnodes by Ntips by Ntips and then rotating to
            # Ntips by Ntips by Nnodes
            V <- array(V, rep(dims, c(1, 2)))
            V <- aperm(V, c(2, 3, 1))
        }
    }

    class(V) <- "VCV.array"
    return(V)
}
