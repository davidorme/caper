#' Calculate and bootstrap phylogenetic diversity measurements.
#'
#' These functions calculate various phylogenetic diversity measures for either
#' a given set of nodes on a tree or for a randomly chosen set of nodes of a
#' given size. The ed.calc function calculates a related species-level
#' measurement of evolutionary distinctness.
#'
#' There are five implemented PD measures: \describe{ \item{Total Branch Length
#' (TBL)}{The sum of all the edge lengths in the subtree given by the tip
#' subset. This measure can be partitioned into the two next measures.}
#' \item{Shared Branch Length (SBL)}{The sum of all edges in the subtree that
#' are shared by more than one tip.} \item{Unique Evolutionary History
#' (UEH)}{The sum of the edge lengths that give rise to only one tip in the
#' subtree.} \item{Length of tip branch lengths (TIPS)}{Length of tip branch
#' lengths (TIPS)}Unlike UEH, this measure does not use the unique paths to
#' each tips on the \strong{subtree} and instead gives the sum of the unique
#' branches leading to the tips on the \strong{complete tree}.  \item{Minimum
#' Spanning Tree (MST)}{The sum of the lengths of the edges for the smallest
#' tree that links the subset tips, excluding any edges below the node of the
#' most recent common ancestor.}}
#'
#' These options are illustrated in the caper package vignette. The pd.calc
#' function returns the PD value for a given set of tips, whereas the
#' pd.bootstrap function returns a vector of PD values for randomly selected
#' sets of tips of a given size.
#'
#' The ed.calc function returns the evolutionary distinctness (ED) metric
#' (Isaac et al, 2007) for the tips of a given phylogeny. The idea behind the
#' ED measure is that the evolutionary history of each branch is shared equally
#' between all tips descending from that branch. Each branch therefore has a
#' per-tip values of the branch length divided by the number of descendants and
#' the ED value for a tip is the sum of those per-tip contributions over the
#' path to the root of the phylogeny. Polytomies inflate apparent ED since the
#' branches of a properly resolved polytomy must be shorter than the branch
#' lengths on the unresolved polytomy. The function provides two correction
#' factors for this: 'isaac' uses a correction factor calibrated from
#' simulations and 'mooers' uses empirical predictions from a pure birth model.
#'
#' @aliases pd.calc pd.bootstrap ed.calc
#' @param cm A object of class 'clade matrix'. Alternatively an object of class
#' 'phylo', which will be converted to a clade.matrix.
#' @param tip.subset An optional vector identifying the subset of tips to use
#' for PD calculations. If no tip.subset is provided the method is applied to
#' the complete phylogeny [Hmm.. this might be undesirable]. Can either be a
#' character vector, in which case the elements are matched against tip labels,
#' or a vector of positive integers in the range 1 to the number of tips, in
#' which case the tips with those numbers are used.
#' @param method One of 'TBL', 'MST', 'UEH', 'SBL', defaulting to 'TBL'. See
#' details.
#' @param root.edge Logical indicating whether to include the root edge length
#' in calculations, defaulting to FALSE.
#' @param ntips A single integer giving the number of tips to be selected.
#' @param reps The number of replicate values to calculate.
#' @param tip.weights A numeric vector containing weights for all the tips in
#' the phylogeny. Each element must be named in order to match weights to the
#' tips.
#' @param polytomy.cf Which correction factor to use for calculating ED at
#' polytomies. One of 'isaac', 'mooers' or 'none'.
#' @return Both pd.calc and pd.bootstrap return a vector containing either a
#' single value for the phylogenetic diversity of a given set of tips or a
#' vector of length 'nrep' containing the pd values for a random set of tips of
#' a given size. The method used is stored in the 'pd.method' attribute of the
#' vector.
#'
#' The ed.calc function returns a list containing: \describe{ \item{branch}{A
#' data frame of the ED contributions arising from each branch.} \item{spp}{A
#' data frame of the summed ED contributions for each species.} }
#' @author David Orme, Gavin Thomas, Nick Isaac
#' @references Faith, DP, Isaac, N. J. B., Turvey, S. T., Collen, B., Waterman,
#' C., and Baillie, J. E. M. (2007). Mammals on the edge: Conservation
#' priorities based on threat and phylogeny. Plos One, 2(3):e296
#' @keywords utilities
#' @examples
#'
#' treeString <- paste("((((A:1,B:1):1.5,C:2.5):0.5,(D:0.6,E:0.6):2.4):0.5,",
#'     "((F:1.9,G:1.9):0.8,(H:1.6,I:1.6):1.1):0.8):0.2;",
#'     sep = ""
#' )
#' tre <- read.tree(text = treeString)
#' clmat <- clade.matrix(tre)
#' tips <- c("A", "C", "D", "E", "G", "H")
#' pd.calc(clmat, tip.subset = tips)
#' pd.calc(clmat, tip.subset = c(1, 3, 4, 5, 7, 8))
#' pd.calc(clmat, tip.subset = tips, root.edge = TRUE)
#'
#' pd.bootstrap(clmat, ntips = 6, reps = 1000, method = "TBL")
#'
#' data(IsaacEtAl)
#' primatesCM <- clade.matrix(primates.tree)
#' primatesED <- ed.calc(primatesCM)
#' @export
pd.calc <- function(cm, tip.subset = NULL, method = "TBL", root.edge = FALSE) {
    # check we have a valid method
    method <- match.arg(method, c("TBL", "MST", "UEH", "SBL", "TIP"))

    # check we have a clade matrix and, if not, get one
    if (!inherits(cm, "clade.matrix")) {
        if (!inherits(cm, "phylo")) {
            warning("Converting phylo object to clade.matrix object")
            cm <- clade.matrix(cm)
        } else {
            stop("pd.calc requires a phylogeny")
        }
    }

    # if requested, drop the root edge
    nSpp <- dim(cm$clade.matrix)[2]
    if (!root.edge) {
        cm$edge.length[nSpp + 1] <- 0
    }

    # subset the tips if requested
    tip_ids <- seq_len(dim(cm$clade.matrix)[2])
    if (!is.null(tip.subset)) {
        # could be names - in which case they must match tip labels
        # could be numbers - in which case they must be in range
        switch(mode(tip.subset),
            "character" = {
                tip.subset <- match(tip.subset, cm$tip.label)
                if (any(is.na(tip.subset))) {
                    stop("Unmatched names in tip.subset")
                }
            },
            "numeric" = {
                if (any(tip.subset %in% tip_ids == FALSE)) {
                    stop(
                        "numeric tip.subset contains outside the ",
                        "range 1 to number of tips"
                    )
                }
            },
            stop("tip.subset must be either a vector of names or numbers")
        )
    } else {
        tip.subset <- tip_ids
    }

    # choose method
    switch(method,
        "TBL" = {
            edge.in.matrix <- cm$clade.matrix[, tip.subset]
            if (is.array(edge.in.matrix)) {
                edge.in.matrix <- rowSums(edge.in.matrix)
            }
            edge.in.matrix <- edge.in.matrix > 0
        },
        "MST" = {
            edge.in.matrix <- cm$clade.matrix[, tip.subset]
            if (is.array(edge.in.matrix)) {
                edge.in.matrix <- rowSums(edge.in.matrix)
            }
            edge.in.matrix <- (
                edge.in.matrix > 0 &
                    edge.in.matrix < length(tip.subset)
            )
        },
        "TIP" = {
            edge.in.matrix <- cm$clade.matrix[, tip.subset]
            if (is.array(edge.in.matrix)) {
                edge.in.matrix <- rowSums(edge.in.matrix)
            }
            edge.in.matrix <- edge.in.matrix > 0
            edge.in.matrix[
                (dim(cm$clade.matrix)[2] + 1):length(edge.in.matrix)
            ] <- FALSE
        },
        "UEH" = {
            edge.in.matrix <- cm$clade.matrix[, tip.subset]
            if (is.array(edge.in.matrix)) {
                edge.in.matrix <- rowSums(edge.in.matrix)
            }
            edge.in.matrix <- edge.in.matrix == 1
        },
        "SBL" = {
            edge.in.matrix <- cm$clade.matrix[, tip.subset]
            if (is.array(edge.in.matrix)) {
                edge.in.matrix <- rowSums(edge.in.matrix)
            }
            edge.in.matrix <- edge.in.matrix > 1
        }
    )

    pd <- sum(cm$edge.len[edge.in.matrix])
    RET <- structure(.Data = pd, pd.method = method)

    return(RET)
}


pd.bootstrap <- function(cm, ntips, reps = 1000, method = "TBL",
                         tip.weights = NULL) {
    # check we have a valid method
    method <- match.arg(method, c("TBL", "MST", "UEH", "SBL", "TIP"))

    # check we have a clade matrix and, if not, get one
    if (!inherits(cm, "clade.matrix")) {
        if (inherits(cm, "phylo")) {
            warning("Converting phylo object to clade.matrix object")
            cm <- clade.matrix(cm)
        } else {
            stop("pd.calc requires a phylogeny")
        }
    }

    # check for sensible sample
    total.nb.tips <- dim(cm$clade.matrix)[2]
    if (!(ntips %in% 1:(total.nb.tips - 1))) {
        stop(
            "'sample' must be a positive integer ",
            "lower than the number of tips"
        )
    }

    # set up the store
    pd.store <- numeric(reps)
    tips <- 1:total.nb.tips

    # if there are weights make sure they go in the right place...
    if (!is.null(tip.weights)) {
        # if the vector is named then match the order to tips
        if (!is.null(names(tip.weights))) {
            wght.match <- match(cm$tip.label, names(tip.weights))

            # this is not elegant but can't work out how to stop and return
            if (any(is.na(wght.match))) {
                warning(
                    "The returned tip labels have no matching ",
                    "named element in tip.weights"
                )
                return(cm$tip.label[is.na(wght.match)])
            }

            tip.weights <- tip.weights[wght.match]
        } else {
            stop(
                "'weights' must be a vector of weights, named ",
                "to match the tip labels"
            )
        }
    }

    # get the pd values
    for (rep in seq(along = pd.store)) {
        which.tips <- sample(tips, ntips, prob = tip.weights)
        pd.store[rep] <- pd.calc(cm, tip.subset = which.tips, method = method)
    }

    return(structure(.Data = pd.store, pd.method = method))
}

ed.calc <- function(cm, polytomy.cf = c("isaac", "mooers", "none")) {
    # Nick Isaac, March 2009 + David Orme 2011 takes the phylogeny and returns a
    # list containing ED scores of a) species and b) branches the polytomy.cf
    # argument specifies which set of polytomies should be applied. There are
    # three options: "isaac" is as the EDGE paper of 2007, based on logarithmic
    # decay with node size (too harsh on large nodes) "mooers" is empirical,
    # based on a pure-birth process "none" : no correction

    # CHANGES TO MAKE:
    # 1) Optional data frame containing IUCN categories
    # 2) Optional data frame containing richness weights - to account for
    #    missing species

    # check we have a clade matrix and, if not, get one

    if (inherits(cm, "phylo")) {
        warning("Converting phylo object to clade.matrix object")
        cm <- clade.matrix(cm)
    } else if (!inherits(cm, "clade.matrix")) {
        stop("pd.calc requires a phylogeny")
    }

    polytomy.cf <- match.arg(polytomy.cf)

    # get raw edge scores (allocates equal proportions of each branch length
    # between all descendent species)
    branch <- data.frame(len = cm$edge.length, nSp = rowSums(cm$clade.matrix))
    branch$ED <- with(branch, len / nSp)

    ## polytomy corrections (the branch lengths at soft polytomies
    ## overestimate the amount of evolution going on)

    # get the group size at the parent nodes (i.e. the number of siblings at a
    # node)
    branch$parent <- cm$edge[, 1][match(rownames(branch), cm$edge[, 2])]
    node.size <- as.data.frame(table(cm$edge[, 1]))
    branch$node.size <- node.size$Freq[match(branch$parent, node.size$Var1)]

    branch$ED.cor <- switch(polytomy.cf,
        "isaac" = {
            with(
                branch,
                ifelse(node.size > 57, 0, ED * (1.081 - 0.267 * log(node.size)))
            )
        },
        "mooers" = {
            with(
                branch,
                ED / node.size *
                    (node.size - 1) / sapply(node.size, function(n) {
                        if (is.na(n)) NA else sum(1 / 2:n)
                    })
            )
        },
        "none" = branch$ED
    )
    branch$ED.cor <- with(branch, ifelse(node.size > 2, ED.cor, ED))

    # get species edge sums
    edge.matrix <- branch$ED.cor * cm$clade.matrix
    spp.ED <- data.frame(
        species = cm$tip.label,
        ED = colSums(edge.matrix, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    return(list(spp = spp.ED, branch = branch))
}
