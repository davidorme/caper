#' Calculates the phylogenetic D statistic
#'
#' Calculates the D value, a measure of phylogenetic signal in a binary trait,
#' and tests the estimated D value for significant departure from both random
#' association and the clumping expected under a Brownian evolution threshold
#' model.
#'
#' The sum of changes in estimated nodal values of a binary trait along edges
#' in a phylogeny (D) provides a measure of the phylogenetic signal in that
#' trait (Fritz and Purvis, 2010). If a trait is highly conserved, with only a
#' basal division between two clades expressing either trait value, then the
#' only change will be along the two daughters at the root. This will give a
#' summed value of 1: the two differences between the root nodal value of 0.5
#' and the ancestors of the 1 and 0 clades. In contrast, if the trait is
#' labile, more differences will be observed and the sum will be higher.
#'
#' This function calculates the observed D for a binary trait on a tree and
#' compares this to the value of D found using an equal number of simulations
#' under each of two models: \describe{ \item{Phylogenetic randomness}{Trait
#' values are randomly shuffled relative to the tips of the phylogeny and D is
#' calculated.} \item{Brownian threshold model}{A continuous trait is evolved
#' along the phylogeny under a Brownian process and then converted to a binary
#' trait using a threshold that reproduces the relative prevalence of the
#' observed trait.} } The value of D depends on phylogeny size - more sister
#' clades yield higher sums - and so the means of the two sets of simulated
#' data are used as calibrations to scale both observed and simulated values of
#' D to set points of 0 (as phylogenetically conserved as expected under a
#' Brownian threshold model) and 1 (random). The value of D can be both smaller
#' than 0 (highly conserved) and greater than 1 (overdispersed) and the
#' distributions of scaled D from the simulations are used to assess the
#' significance of the observed scaled D. The \code{plot} method generates
#' density plots of the distributions of the two simulations relative to the
#' observed D value.
#'
#' \code{rnd.bias} is passed to \code{\link{sample}} as the \code{prob}
#' argument to weight the random shuffles of the observed trait. The weights
#' are not checked for validity.
#'
#' @aliases phylo.d print.phylo.d summary.phylo.d plot.phylo.d
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
#' @param x An object of class 'phylo.d'
#' @param object An object of class 'phylo.d'
#' @param bw The bandwidth to be used for the density plots
#' @param list() Further arguments to print and summary methods
#' @return Returns an object of class 'phylo.d', which is a list of the
#' following:
#'
#' \item{DEstimate}{The estimated D value}
#' \item{Pval1}{
#'      A p value, giving the result of testing whether D is significantly
#'      different from one
#' }
#' \item{Pval0}{
#'      A p value, giving the result of testing whether D is significantly
#'      different from zero
#' }
#' \item{Parameters}{
#'      A list of the Observed, MeanRandom and MeanBrownian sums of sister-clade
#'      differences
#' }
#' \item{Permutations}{
#'      A list with elements random and brownian, containing the sums of
#'      sister-clade differences from random permutations and simulations of
#'      Brownian evolution under a threshold model
#' }
#' \item{NodalVals}{
#'      A list with the elements observed, random and brownian, containing the
#'      nodal values estimated for the observed trait and permutations. The
#'      values are as matrices with rows labelled by the node names in the
#'      comparative data object.
#' }
#' \item{binvar}{The binary variable used}
#' \item{phyName}{The name of the phylogeny object used}
#' \item{dsName}{The name of the dataframe used}
#' \item{nPermut}{The number of permutations used}
#' \item{rnd.bias}{
#'      If a bias was introduced to the calculation of the random distribution,
#'      the bias used, else \code{NULL}
#' }
#'
#' @author Susanne Fritz <Susanne.Fritz@senckenberg.de> and David Orme
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
#' redPhyloD <- phylo.d(BritishBirds, binvar = Red_list)
#' print(redPhyloD)
#' plot(redPhyloD)
#' @export
phylo.d <- function(data, phy, names.col, binvar,
                    permut = 1000, rnd.bias = NULL) {
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

    # get the variable out
    ds <- data$data[, bininds]
    if (any(is.na(ds))) stop("'", binvar, "' contains missing values.")

    # sort out character variables
    if (is.character(ds)) ds <- as.factor(ds)

    # test for binary states
    if (length(unique(ds)) > 2) {
        stop("'", binvar, "' contains more than two states.")
    }
    if (length(unique(ds)) < 2) {
        stop("'", binvar, "' only contains a single state.")
    }


    # get proportion and table of classes
    propStates <- unclass(table(ds))
    propState1 <- propStates[1] / sum(propStates)
    names(dimnames(propStates)) <- binvar

    # convert factors to numeric for calculation
    if (is.factor(ds)) ds <- as.numeric(ds)

    # check for a number
    if (!is.numeric(permut)) (stop("'", permut, "' is not numeric."))

    # look for probaility weights argument and get its value if found
    if (!is.null(rnd.bias)) {
        rnd.bias <- deparse(substitute(rnd.bias))
        rnd.ind <- match(rnd.bias, names(data$data))
        if (is.na(rnd.ind)) {
            stop("'", rnd.bias, "' is not a variable in data.")
        }
        rnd.bias <- data$data[, rnd.bias]
    }

    # check tree branch lengths
    el <- data$phy$edge.length
    elTip <- data$phy$edge[, 2] <= length(data$phy$tip.label)

    if (any(el[elTip] == 0)) {
        stop(
            "Phylogeny contains pairs of tips on zero branch lengths, ",
            "cannot currently simulate"
        )
    }
    if (any(el[!elTip] == 0)) {
        stop("Phylogeny contains zero length internal branches. Use di2multi.")
    }

    ## This is rewritten away from the original version with internal functions
    ##  - structure was slowing and the functions aren't externalised ever

    ## Random Association model random data
    ##  - with weighted shuffling if weights are given
    ds.ran <- replicate(permut, sample(ds, prob = rnd.bias))

    ## Brownian Threshold model random data

    ## there was a call to lambdaTree(phy,1) - why???
    ## - get the variance covariance for the tree
    if (is.null(data$vcv)) {
        vcv <- VCV.array(data$phy)
    } else {
        vcv <- data$vcv
    }

    # Simulate traits up the tree, no class to avoid throwing method dispatch
    ds.phy <- mvtnorm::rmvnorm(permut, sigma = unclass(vcv))
    ds.phy <- as.data.frame(t(ds.phy))

    ## - find the threshold in each variable.
    ## - quantile interpolates between values
    ds.phy.thresh <- apply(ds.phy, 2, stats::quantile, propState1)

    ## sweep out the thresholds
    ds.phy <- sweep(ds.phy, 2, ds.phy.thresh, "<")
    ds.phy <- as.numeric(ds.phy) ## bah! kills dims so reinstate
    dim(ds.phy) <- dim(ds.ran)

    ## Get change along edges

    ## insert observed and set dimnames for contrCalc
    ds.ran <- cbind(Obs = ds, ds.ran)
    ds.phy <- cbind(Obs = ds, ds.phy)
    dimnames(ds.ran) <- dimnames(ds.phy) <- list(
        data$phy$tip.label, c("Obs", paste("V", 1:permut, sep = ""))
    )

    ## being careful with the edge order - pre-reorder the phylogeny
    ## because the method won't reorder an already matching order.
    ## Plus we need the pruningwise order later.
    phy <- reorder(data$phy, "pruningwise")

    # now run that through the contrast engine
    # - in fact, the change calculation requires a tree traversal to compare
    #   change along the edges from the nodal values of the daughters to the
    #   parent and this traversal is what contrCalc does. So create a new
    #   contrCalc method.
    ds.ran.cc <- contrCalc(
        vals = ds.ran, phy = phy, ref.var = "V1",
        picMethod = "phylo.d", crunch.brlen = 0
    )
    ds.phy.cc <- contrCalc(
        vals = ds.phy, phy = phy, ref.var = "V1",
        picMethod = "phylo.d", crunch.brlen = 0
    )

    ## get sums of change and distributions

    ransocc <- colSums(ds.ran.cc$contrMat)
    physocc <- colSums(ds.phy.cc$contrMat)
    # double check the observed, but only to six decimal places or you can get
    # floating point errors
    if (round(ransocc[1], digits = 6) != round(physocc[1], digits = 6)) {
        stop("Problem with character change calculation in phylo.d")
    }
    obssocc <- ransocc[1]
    ransocc <- ransocc[-1]
    physocc <- physocc[-1]

    soccratio <- (obssocc - mean(physocc)) / (mean(ransocc) - mean(physocc))
    soccpval1 <- sum(ransocc < obssocc) / permut
    soccpval0 <- sum(physocc > obssocc) / permut


    dvals <- list(
        DEstimate = soccratio, Pval1 = soccpval1, Pval0 = soccpval0,
        Parameters = list(
            Observed = obssocc,
            MeanRandom = mean(ransocc), MeanBrownian = mean(physocc)
        ),
        StatesTable = propStates,
        Permutations = list(random = ransocc, brownian = physocc),
        NodalVals = list(
            observed = ds.ran.cc$nodVal[, 1, drop = FALSE],
            random = ds.ran.cc$nodVal[, -1, drop = FALSE],
            brownian = ds.phy.cc$nodVal[, -1, drop = FALSE]
        ),
        binvar = binvar, data = data, nPermut = permut, rnd.bias = rnd.bias
    )

    class(dvals) <- "phylo.d"
    return(dvals)
}

#' @describeIn phylo.d Print a summary of a phylo.d object
#' @export
print.phylo.d <- function(x, ...) {
    summary(x)
}

#' @describeIn phylo.d Print a summary of a phylo.d object
#' @export
summary.phylo.d <- function(object, ...) {
    cat(
        "\nCalculation of D statistic for the phylogenetic ",
        "structure of a binary variable\n"
    )
    cat("\n  Data : ", object$data$data.name)
    cat("\n  Binary variable : ", object$binvar)
    stCounts <- paste(
        names(object$StatesTable), " = ", object$StatesTable,
        sep = ""
    )
    cat("\n  Counts of states: ", stCounts[1])
    cat("\n                    ", stCounts[2])
    cat("\n  Phylogeny : ", object$data$phy.name)
    cat("\n  Number of permutations : ", object$nPermut)

    cat("\n\nEstimated D : ", object$DEstimate)
    if (is.null(object$rnd.bias)) {
        cat(
            "\nProbability of E(D) resulting from no (random) ",
            "phylogenetic structure : ", object$Pval1
        )
    } else {
        cat(
            "\nProbability of E(D) resulting from no (biased random) ",
            "phylogenetic structure : ", object$Pval1
        )
    }
    cat(
        "\nProbability of E(D) resulting from Brownian ",
        "phylogenetic structure    : ", object$Pval0
    )
    cat("\n\n")
}

#' @describeIn phylo.d Plot a phylo.d object
#' @export
plot.phylo.d <- function(x, bw = 0.02, ...) {
    brownian <- x$Permutations$brownian
    random <- x$Permutations$random

    centre <- x$Parameters$MeanBrownian
    scale <- x$Parameters$MeanRandom - centre

    brownian <- (brownian - centre) / scale
    random <- (random - centre) / scale
    obs <- (x$Parameters$Observed - centre) / scale


    xlim <- range(obs, brownian, random)

    bdens <- as.data.frame(stats::density(brownian, bw = bw)[c("y", "x")])
    rdens <- as.data.frame(stats::density(random, bw = bw)[c("y", "x")])

    ylim <- range(bdens$y, rdens$y)
    plot(y ~ x,
        data = bdens, xlim = xlim, ylim = ylim, type = "l", col = "blue",
        xlab = "D value", ylab = "Density", ...
    )
    graphics::lines(y ~ x, data = rdens, col = "red", ...)

    graphics::abline(v = c(0, 1, obs), col = c("blue", "red", "black"))
}
