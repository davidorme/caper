#' Calculates the phylogenetic D statistic across clades within a phylogeny
#'
#' Calculates the D value, a measure of phylogenetic signal in a binary trait,
#' and tests the estimated D value for significant departure from both random
#' association and the clumping expected under a Brownian evolution threshold
#' model. Does this across clades within a phylogeny.
#'
#' A wrapper function for \code{\link{phylo.d}}, calculating D values for
#' clades within a given dataset. These clades can be filtered according to the
#' number of species and nodes using the arguments above. See
#' \code{\link{phylo.d}} for more details on the method itself.
#'
#' Any clades for which there is no variation in the binary variable have
#' \code{NA} values for all of the below slots.
#'
#' @aliases phylo.d.subset print.phylo.d.subset summary.phylo.d.subset
#' @param data A 'comparative.data' or 'data.frame' object.
#' @param phy An object of class 'phylo', required when data is not a
#' 'comparative.data' object.
#' @param names.col A name specifying the column in 'data' that matches rows to
#' tips in 'phy', required when data is not a 'comparative.data' object.
#' @param binvar The name of the variable in \code{data} holding the binary
#' variable of interest.
#' @param permut Number of permutations to be used in the randomisation test.
#' @param rnd.bias An optional name of a variable in \code{data} holding
#' probability weights to bias the generation of the random distribution. See
#' 'destails'
#' @param verbose Logical; do you want to know how many clades are being
#' assessed, and see when each is being assessed?
#' @param min.tips The minimum number of tips a clade should have for it to
#' have a D value calculated. Defaults to 1 (i.e. no limit).
#' @param max.tips The maximum number of species a clade should have for it to
#' have a D value calculated. Defaults to the number of species in the whole
#' phylogeny (i.e. no limit).
#' @param min.nodes The minimum number of nodes a clade should have for it to
#' have a D value calculated. Defaults to 1 (i.e. no limit).
#' @param max.nodes The maximum number of nodes a clade should have for it to
#' have a D value calculated. Defaults to the number of nodes in the whole
#' phylogeny (i.e. no limit).
#' @param object An object of class 'phylo.d.subset'
#' @param ... Further arguments to print and summary methods
#' @return Returns an object of class 'phylo.d.subset', which is a list of the
#' following:
#'
#' \item{raw}{
#'      A list of the raw output from \code{\link{phylo.d}} for each clade
#' }
#' \item{DEstimate}{A vector of the estimated D values}
#' \item{Pval1}{
#'      A vector of p values, giving the result of testing whether D is
#'      significantly different from one, for each clade
#' }
#' \item{Pval0}{
#'      A vector of p values, giving the result of testing whether D is
#'      significantly different from zero, for each clade
#' }
#' \item{phy.depth}{
#'      A numeric vector giving the age of the clade for which each value was
#'      calculated
#' }
#' @author Susanne Fritz (SFritz@@bio.ku.dk), Will Pearse and David Orme
#' @references Fritz, S. A. and Purvis, A. (2010). Selectivity in mammalian
#' extinction risk and threat types: a new measure of phylogenetic signal
#' strength in binary traits. Conservation Biology, 24(4):1042-1051.
#' @keywords utilities htest
#' @examples
#'
#' data(BritishBirds)
#' BritishBirds <- comparative.data(
#'     BritishBirds.tree, BritishBirds.data, binomial
#' )
#' # Look at big clades only
#' \dontrun{
#' bigClades <- phylo.d.subset(
#'     BritishBirds,
#'     binvar = Red_list, verbose = TRUE, min.tips = 10, min.nodes = 5
#' )
#' print(bigClades)
#' }
#' @export
phylo.d.subset <- function(data, phy, names.col, binvar, permut = 1000,
                           rnd.bias = NULL, min.tips = 1,
                           max.tips = length(data$phy$tip.label),
                           min.nodes = 1, max.nodes = data$phy$Nnode,
                           verbose = FALSE) {
    # - test to see if there is a comparative data object and if not then
    #   retrofit the remaining arguments into a comparative data object.
    if (!missing(data)) {
        if (!inherits(data, "comparative.data")) {
            if (missing(names.col)) stop("names column is missing")
            names.col <- deparse(substitute(names.col))
            data <- caicStyleArgs(data = data, phy = phy, names.col = names.col)
        }
    }

    # look for binary variable
    binvar <- deparse(substitute(binvar))
    bininds <- match(binvar, names(data$data))
    if (is.na(bininds)) (stop("'", binvar, "' is not a variable in data."))

    # get the variable out and do a general test for binarity
    ds <- data$data[, bininds]
    if (length(unique(ds)) != 2) {
        stop("'", binvar, "' doesn't contain two states.")
    }
    if (any(is.na(ds))) {
        stop("'", binvar, "' contains missing values.")
    }

    # check for a number
    if (!is.numeric(permut)) {
        stop("'", permut, "' is not numeric.")
    }

    # look for probaility weights argument and get its value if found
    if (missing(rnd.bias)) {
        rnd.bias <- NULL
    } else {
        rnd.bias <- deparse(substitute(rnd.bias))
        rnd.ind <- match(rnd.bias, names(data$data))
        if (is.na(rnd.ind)) {
            stop("'", rnd.bias, "' is not a variable in data.")
        }
        rnd.bias <- data$data[, rnd.bias]
    }

    # make the subtrees
    phy.subtrees <- ape::subtrees(data$phy)
    if (verbose) phy.max.subtrees <- length(phy.subtrees)

    # filter them
    phy.subtrees <- phy.subtrees[
        sapply(
            phy.subtrees,
            function(x) {
                max.tips >= length(x$tip.label) &
                    length(x$tip.label) >= min.tips &
                    x$Nnode <= max.nodes &
                    min.nodes <= x$Nnode
            }
        )
    ]

    # prepare output
    output.raw <- vector(mode = "list", length(phy.subtrees))
    output.D <-
        output.depth <-
        output.P0 <-
        output.P1 <-
        output.tips <-
        output.nodes <-
        output.bin.freq <- numeric(length(phy.subtrees))

    ## TO-DO:
    #  - make a 'clean' version where the raw output isn't saved
    #    (when used on a big phylogeny this could be big)
    #  - allow filtering according to presence:absence
    #  - record node labels/tip labels within each subset in a convenient format

    # run phylo.d on each subtree
    if (verbose) {
        cat(
            "\nCalculating D values for ", length(phy.subtrees),
            " out of a possible ", phy.max.subtrees, "clades"
        )
    }
    for (i in seq(along = phy.subtrees)) {
        if (verbose) cat(".")
        # make a temporary comparative.data.frame and check it's worth passing
        # through phylo.d - if not, fill it with NAs
        # - assuming most users would want to know about subsets that fail this
        #   test---I would!
        t.data <- data[rownames(data$data) %in% phy.subtrees[[i]]$tip.label, ]
        if (length(unique(t.data$data[, bininds])) != 2) {
            output.raw[[i]] <-
                output.depth[i] <-
                output.P0[i] <-
                output.P1[i] <- NA

            next
        }
        # get the raw output
        output.raw[[i]] <- eval(substitute(
            phylo.d(t.data, binvar = XXX),
            list(XXX = as.name(names(data$data)[bininds]))
        ))
        # get the summaries
        output.D[i] <- output.raw[[i]]$DEstimate
        output.P0[i] <- output.raw[[i]]$Pval0
        output.P1[i] <- output.raw[[i]]$Pval1
        clade.mat <- clade.matrix(phy.subtrees[[i]])
        output.depth[i] <- sum(clade.mat$edge.length[
            as.logical(clade.mat$clade.matrix[, 1])
        ])
        output.tips[i] <- length(phy.subtrees[[i]]$tip.label)
        output.nodes[i] <- phy.subtrees[[i]]$Nnode
        output.bin.freq[i] <- sum(t.data$data[, bininds])
    }

    output <- list(
        raw = output.raw,
        DEstimate = output.D,
        Pval1 = output.P1,
        Pval0 = output.P0,
        phy.depth = output.depth,
        tips = output.tips,
        nodes = output.nodes,
        bin.freq = output.bin.freq
    )
    class(output) <- "phylo.d.subset"
    return(output)
}

print.phylo.d.subset <- function(x, ...) {
    summary(x)
}

##############################################
## TO-DO:
## MAKE THIS MORE PLEASANT TO INTERPRET
## ADD IN PLOTS OF D vs. PHYLOGENETIC DEPTH
##############################################

#' @describeIn phylo.d.subset Print a summary of phylo.d.subset object
#' @export
summary.phylo.d.subset <- function(object, ...) {
    cat(
        "\nCalculation of D statistic for the ",
        "phylogenetic structure of a binary variable"
    )
    cat("\n...across multiple clades within a phylogeny\n")
    cat("\n  Data : ", object$raw[[1]]$data$data.name)
    cat("\n  Binary variable : ", object$raw[[1]]$binvar)
    cat("\n  Phylogeny : ", object$raw[[1]]$data$phy.name)
    cat("\n  Number of permutations : ", object$raw[[1]]$nPermut)

    cat("\n\nEstimated D values: \n")
    t <- format(round(object$DEstimate, digits = 2), trim = TRUE)
    cat(format(object$DEstimate, width = nchar(t)))
    if (is.null(object$raw[[1]]$rnd.bias)) {
        cat(
            "\nProbabilities of E(D) resulting from no (random) ",
            "phylogenetic structure : \n"
        )
    } else {
        cat(
            "\nProbability of E(D) resulting from no (biased random) ",
            "phylogenetic structure : \n"
        )
    }
    cat(format(object$Pval0, width = nchar(t)))
    cat(
        "\nProbabilities of E(D) resulting from Brownian ",
        "phylogenetic structure    : \n"
    )
    cat(format(object$Pval1, width = nchar(t)))
    cat(
        "\nAges of clades                                  ",
        "                      : \n"
    )
    cat(format(object$phy.depth, width = nchar(t)))
    cat("\n\n")
}
