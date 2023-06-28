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
    ds.phy.thresh <- apply(ds.phy, 2, quantile, propState1)

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

print.phylo.d <- function(x, ...) {
    summary(x)
}

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

plot.phylo.d <- function(x, bw = 0.02, ...) {
    brownian <- x$Permutations$brownian
    random <- x$Permutations$random

    centre <- x$Parameters$MeanBrownian
    scale <- x$Parameters$MeanRandom - centre

    brownian <- (brownian - centre) / scale
    random <- (random - centre) / scale
    obs <- (x$Parameters$Observed - centre) / scale


    xlim <- range(obs, brownian, random)

    bdens <- as.data.frame(density(brownian, bw = bw)[c("y", "x")])
    rdens <- as.data.frame(density(random, bw = bw)[c("y", "x")])

    ylim <- range(bdens$y, rdens$y)
    plot(y ~ x,
        data = bdens, xlim = xlim, ylim = ylim, type = "l", col = "blue",
        xlab = "D value", ylab = "Density", ...
    )
    lines(y ~ x, data = rdens, col = "red", ...)

    abline(v = c(0, 1, obs), col = c("blue", "red", "black"))
}
