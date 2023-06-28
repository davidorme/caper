VCV.array <- function(phy, dim = 2, compact = TRUE) {
    ## turns a phylogeny into a 3d array similar to a VCV matrix
    ## but keeping each beanch length separate. This is useful for
    ## handling branch length transformations in functions where
    ## VCVs are used to handle phylogenetic structure

    ## rewritten to use new ape, via the clade matrix structure

    if (inherits(phy, "phylo")) {
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

            # gives a (Nnodes * Ntips) by Ntips matrix, which needs reshaping
            # into an array of Nnodes by Ntips by Ntips and then rotating to
            # Ntips * Ntips * Nnodes
            V <- array(V, rep(dims, c(1, 2)))

            V <- aperm(V, c(2, 3, 1))
        }
    }

    class(V) <- "VCV.array"
    return(V)
}
