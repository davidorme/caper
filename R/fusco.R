#' Imbalance statistics using Fusco and Cronk's method.
#' 
#' Fusco and Cronk (1995) described a method for testing the imbalance of
#' phylogenetic trees based on looking at the distribution of I. I is
#' calculated using the number of tips descending from each side of a
#' bifurcating node using the formula I = (B-m)/(M-m) and is bounded between 0
#' (a perfectly balanced node) and 1 (maximum imbalance). B is the larger
#' number of tips descending from each branch, M is the maximum size of this
#' larger group (i.e. a 1 : (S-1) split, where S is the total number of
#' descendent tips), and m is the minimum size of the larger group (ceiling of
#' S/2). The method can cope with small proportions of polytomies in the
#' phylogeny and these are not used in calculating balance statistics. It can
#' also incorporate information about species richness at the tips of the
#' phylogeny and can therefore be used to distinguish between an unbalanced
#' topology and the unbalanced distribution of diversity at the tips of a
#' phylogeny.
#' 
#' Purvis et al. (2002) demonstrated that I is not independent of the node size
#' S, resulting in a bias to the expected median of 0.5. They proposed a
#' modification (I') that corrects this to give a statistic with an expected
#' median of 0.5 regardless of node size. The defaults in this function perform
#' testing of imbalance using I', but it is also possible to use the original
#' measure proposed by Fusco and Cronk (1995).
#' 
#' I is calculated only at bifurcating nodes giving rise to more than 3 tips
#' (or more than 3 species at the tips): nodes with three or fewer descendants
#' have no variation in I and are not informative in assessing imbalance. The
#' expected distribution of the nodal imbalance values between 0 and 1 is
#' theoretically uniform under a Markov null model. However, the range of
#' possible I values at a node is constrained by the number of descendent
#' species. For example, for a node with 8 species, only the values 0, 0.33,
#' 0.66, 1 are possible, corresponding to 4:4, 5:3, 6:2 and 7:1 splits (Fusco
#' and Cronk, 1995). As node size increases, this departure from a uniform
#' distribution decreases. The plot method incorporates a correction, described
#' by Fusco and Cronk (1995), that uses the distribution of all possible splits
#' at each node to characterize and correct for the departure from uniformity.
#' 
#' The randomization option generates confidence intervals around the mean I'.
#' 
#' @aliases fusco.test summary.fusco print.fusco plot.fusco
#' @param phy An object of class 'comparative.data' or of class 'phylo'.
#' @param data A data frame containing species richness values. Not required if
#' phy is a 'comparative.data' object.
#' @param names.col A variable in \code{data} identifying tip labels. Not
#' required if phy is a 'comparative.data' object.
#' @param rich A variable identifying species richness.
#' @param tipsAsSpecies A logical value. If TRUE, a species richness column
#' need not be specified and the tips will be treated as species.
#' @param randomise.Iprime Use a randomization test on calculated I' to
#' generate confidence intervals.
#' @param reps Number of replicates to use in simulation or randomization.
#' @param conf.int Width of confidence intervals required.
#' @param x,object An object of class 'fusco'.
#' @param correction Apply the correction described in Appendix A of Fusco and
#' Cronk (1995) to the histogram of nodal imbalance.
#' @param nBins The number of bins to be used in the histogram of nodal
#' imbalance.
#' @param right Use right or left open intervals in plotting the distribution
#' and calculating the correction
#' @param I.prime Plot distribution of I' or I.
#' @param plot If changed to FALSE, then the plot method does not plot the
#' frequency histogram. Because the method invisibly returns a table of
#' histogram bins along with the observed and corrected frequencies, this isn't
#' as dim an option as it sounds.
#' @param ... Further arguments to generic methods
#' @return The function \code{fusco.test} produces an object of class 'fusco'
#' containing: \item{observed}{A data frame of informative nodes showing nodal
#' imbalance statistics. If the phylogeny has labelled nodes, then the node
#' names are also returned in this data frame.} \item{median}{The median value
#' of I.} \item{qd}{The quartile deviation of I.} \item{tipsAsSpecies}{A
#' logical indicating whether the tips of the trees were treated as species or
#' higher taxa.} \item{nInformative}{The number of informative nodes.}
#' \item{nSpecies}{The number of species distributed across the tips.}
#' \item{nTips}{The number of tips.} \item{reps}{The number of replicates used
#' in randomization.} \item{conf.int}{The confidence levels used in
#' randomization.} If \code{randomise.Iprime} is TRUE, or the user calls
#' \code{fusco.randomize} on a 'fusco' object, then the following are also
#' present.  \item{randomised}{A data frame of mean I' from the randomized
#' observed values.} \item{rand.mean}{A vector of length 2 giving confidence
#' intervals in mean I'.}
#' @author David Orme, Andy Purvis
#' @references Fusco, G. & Cronk, Q.C.B. (1995) A New Method for Evaluating the
#' Shape of Large Phylogenies. J. theor. Biol. 175, 235-243
#' 
#' Purvis A., Katzourakis A. & Agapow, P-M (2002) Evaluating Phylogenetic Tree
#' Shape: Two Modifications to Fusco & Cronk's Method. J. theor. Biol. 214,
#' 93-103.
#' @keywords utilities htest
#' @examples
#' 
#' data(syrphidae)
#' syrphidae <- comparative.data(phy=syrphidaeTree, dat=syrphidaeRich, names.col=genus)
#' summary(fusco.test(syrphidae, rich=nSpp))
#' summary(fusco.test(syrphidae, tipsAsSpecies=TRUE))
#' plot(fusco.test(syrphidae, rich=nSpp))
#' 
fusco.test <- function(phy, data, names.col, rich, tipsAsSpecies = FALSE,
                       randomise.Iprime = TRUE, reps = 1000, conf.int = 0.95) {
    # want to be able to just run it on a phylogeny for the topology
    if (inherits(phy, "phylo") && missing(data)) {
        tipsAsSpecies <- TRUE
        data <- data.frame(
            nSpp = rep(1, length(phy$tip.label)),
            tips = phy$tip.label
        )
        phy <- comparative.data(phy = phy, data = data, names.col = "tips")
    } else if (inherits(phy, "phylo") && !missing(data)) {
        if (missing(names.col)) stop("Names column not specified")
        names.col <- deparse(substitute(names.col))
        # old style use with data and phylogeny separated
        # let comparative.data() match and handle class checking
        phy <- eval(substitute(
            comparative.data(
                phy = phy, data = data, names.col = XXX
            ), list(XXX = names.col)
        ))
    } else if (!inherits(phy, "comparative.data")) {
        stop("phy must be a phylo or comparative.data object.")
    }

    # check for the richness column
    if (!tipsAsSpecies) {
        if (missing(rich)) {
            stop("The name of a column of richness values must be provided")
        } else {
            rich <- deparse(substitute(rich))
            if (!rich %in% names(phy$data)) {
                stop(
                    "The column '", rich, "' was not found in the data from '",
                    phy$data.name, "'."
                )
            }
        }
    }

    # get into pruningwise order
    # method dispatch works on either 'phylo' or 'comparative.data'
    phy <- reorder(phy, "pruningwise")
    if (tipsAsSpecies) {
        rich <- rep(1, length(phy$phy$tip.label))
    } else {
        rich <- phy$data[, rich, drop = TRUE] # vector of values
    }

    nSpecies <- sum(rich)
    nTips <- length(rich)

    # calculate fusco
    intNodes <- unique(phy$phy$edge[, 1])
    nTip <- length(phy$phy$tip.label)
    nNode <- phy$phy$Nnode
    rich <- c(rich, rep(NA, nNode))

    # Data store
    observed <- data.frame(
        polytomy = logical(nNode), N1 = numeric(nNode),
        N2 = numeric(nNode), row.names = intNodes
    )

    # loop tips
    for (ind in seq(along = intNodes)) {
        # grab the daughters of the node
        daughters <- phy$phy$edge[, 2][phy$phy$edge[, 1] == intNodes[ind]]
        richD <- rich[daughters]

        if (length(daughters) > 2) {
            observed$polytomy[ind] <- TRUE
        } else {
            observed[ind, 2:3] <- richD
        }

        rich[intNodes[ind]] <- sum(richD)
    }

    observed$S <- with(observed, N1 + N2)
    observed <- observed[
        !observed$polytomy & observed$S > 3, names(observed) != "polytomy"
    ]

    observed$B <- with(observed, pmax(N1, N2))
    observed$M <- observed$S - 1
    observed$m <- ceiling(observed$S / 2)

    observed$I <- with(observed, (B - m) / (M - m))
    observed$S.odd <- (observed$S %% 2) == 1
    observed$w <- with(
        observed,
        ifelse(S.odd, 1, ifelse(I > 0, M / S, 2 * M / S))
    )
    observed$I.w <- with(observed, (I * w) / mean(w))
    observed$I.prime <- with(observed, ifelse(S.odd, I, I * M / S))

    # add node labels if present [Arne Mooers]
    if (!is.null(phy$phy$node.label)) {
        observed$node <- phy$phy$node.label[
            as.numeric(rownames(observed)) - nTip
        ]
    }

    ret <- list(
        observed = observed,
        median.I = median(observed$I),
        mean.Iprime = mean(observed$I.prime),
        qd = IQR(observed$I) / 2,
        tipsAsSpecies = tipsAsSpecies,
        nInformative = dim(observed)[1],
        nSpecies = nSpecies,
        nTips = nTips
    )

    class(ret) <- "fusco"

    if (randomise.Iprime) {
        expFun <- function(x) {
            y <- runif(length(x))
            rand.I.prime <- ifelse(y > 0.5, x, 1 - x)
            ret <- mean(rand.I.prime)
            return(ret)
        }

        randomised <- with(ret, replicate(reps, expFun(observed$I.prime)))
        randomised <- as.data.frame(randomised)
        names(randomised) <- "mean"
        rand.twotail <- c((1 - conf.int) / 2, 1 - (1 - conf.int) / 2)
        rand.mean <- quantile(randomised$mean, rand.twotail)

        ret <- c(
            ret,
            list(
                randomised = randomised, rand.mean = rand.mean,
                reps = reps, conf.int = conf.int
            )
        )
    }

    class(ret) <- "fusco"
    return(ret)
}


print.fusco <- function(x, ...) {
    print(x$observed)
}

summary.fusco <- function(object, ...) {
    cat("Fusco test for phylogenetic imbalance\n\n")

    cat(
        "  Tree with", object$nInformative, "informative nodes and",
        object$nTips, "tips.\n"
    )
    if (object$tipsAsSpecies) {
        cat("  Tips are treated as species.\n\n")
    } else {
        cat("  Tips are higher taxa containing", object$nSpecies, "species.\n")
    }

    if (!is.null(object$randomised) || !is.null(object$simulated)) {
        cat(
            " ", sprintf("%2.1f%%", object$conf.int * 100),
            "confidence intervals around 0.5 randomised using",
            object$reps, "replicates.\n"
        )
    }
    cat("\n")

    cat("  Mean I prime:", round(object$mean.Iprime, 3))
    if (!is.null(object$randomised)) {
        cat(sprintf(" [%1.3f,%1.3f]", object$rand.mean[1], object$rand.mean[2]))
    }
    cat("\n")

    cat("  Median I:", round(object$median.I, 3))
    if (!is.null(object$simulated)) {
        cat(sprintf(
            " [%1.3f,%1.3f]", object$sim.median[1],
            object$sim.median[2]
        ))
    }
    cat("\n")
    cat("  Quartile deviation in I:", round(object$qd, 3))
    if (!is.null(object$simulated)) {
        cat(sprintf(" [%1.3f,%1.3f]", object$sim.qd[1], object$sim.qd[2]))
    }
    cat("\n")

    print(wilcox.test(object$observed$I.prime, mu = 0.5))
}

plot.fusco <- function(x, correction = TRUE, nBins = 10, right = FALSE,
                       I.prime = TRUE, plot = TRUE, ...) {
    breaks <- seq(0, 1, length = nBins + 1)
    if (I.prime) {
        fuscoDist <- hist(x$observed$I.prime,
            breaks = breaks,
            plot = FALSE, right = right
        )
        xLab <- "Nodal imbalance score (I')"
    } else {
        fuscoDist <- hist(x$observed$I,
            breaks = breaks,
            plot = FALSE, right = right
        )
        xLab <- "Nodal imbalance score (I)"
    }

    # hijack the interval notation code
    interv <- levels(cut(0.5, breaks = breaks, right = right))
    RET <- data.frame(
        imbalance = interv,
        observedFrequency = fuscoDist$density / nBins
    )

    if (correction == TRUE) {
        # find out the distribution of possible values for
        # each node given the number of classes
        allPossI <- function(S, I.prime) {
            m <- ceiling(S / 2)
            RET <- (seq(from = m, to = S - 1) - m) / ((S - 1) - m)
            if (I.prime && (S %% 2) == 1) {
                RET <- RET * (S - 1) / S
            }
            return(RET)
        }

        distrib <- sapply(x$observed$S, FUN = allPossI, I.prime = I.prime)
        distrib <- sapply(
            distrib,
            function(x) {
                hist(x, breaks = breaks, plot = FALSE, right = right)$density
            }
        )
        distrib <- distrib / nBins
        correction <- 1 / nBins - distrib
        correction <- rowMeans(correction)

        RET$correction <- correction
        RET$correctedFrequency <- RET$observedFrequency + RET$correction

        # hijack the histogram plotting
        if (plot) {
            plot(
                structure(
                    list(
                        density = RET$correctedFrequency,
                        breaks = breaks
                    ),
                    class = "histogram"
                ),
                freq = FALSE, ylab = "Corrected Frequency", main = "",
                xlab = xLab
            )
        }
    } else {
        # hijack the histogram plotting
        if (plot) {
            plot(
                structure(
                    list(
                        density = RET$observedFrequency,
                        breaks = breaks
                    ),
                    class = "histogram"
                ),
                freq = FALSE, ylab = "Observed Frequency", main = "",
                xlab = xLab
            )
        }
    }

    if (plot) {
        if (I.prime) {
            abline(v = x$mean.Iprime)
            if (!is.null(x$rand.mean)) abline(v = x$rand.mean, col = "red")
        } else {
            abline(v = median(x$observed$I))
            if (!is.null(x$sim.median)) abline(v = x$sim.median, col = "red")
        }
    }

    invisible(RET)
}
