#' Comparative analysis using the crunch algorithm.
#' 
#' Calculate a linear model using the crunch algorithm.
#' 
#' This function implements the 'crunch' algorithm for modelling the
#' relationship between variables that are phylogenetically non-independent.
#' The method was first described by Felsenstein (1985) and subsequently
#' extended to permit the use of phylogenies with polytomies by Pagel (1992).
#' This method was previously implemented in the Mac Classic computer programs
#' CAIC, written by Andy Purvis, Andy Rambaut (Purvis and Rambaut, 1995) and
#' updated by Nick Isaac and Paul-Michael Agapow.
#' 
#' The user provides a comparative dataset. The formula specifies the model to
#' be fitted and contrasts are calculated in those variables. The specified
#' reference variable is used to ensure that contrasts for multivariate models
#' are calculated in a consistent direction at each node. The function
#' \code{crunch()} acts as a data preparation wrapper for the function
#' \code{contrCalc()}, which is not intended to be directly called by users.
#' Missing data can be present in the explanatory variables: the algorithm
#' makes use of the complete data available at each node as was the case with
#' CAIC.
#' 
#' The resulting table of contrasts is then used to fit the specified model -
#' note that the intercept is automatically dropped from the model if present,
#' following REF HERE.
#' 
#' Contrasts at polytomies are calculated following Pagel (1992). The
#' descendants from the node are split into two groups based on whether they
#' are above or below the group mean in the reference variable. If there is no
#' variation in the reference variable, then a 1:(N-1) split is used. Weighted
#' means in the variables are then calculated for each subgroup and a contrast
#' is calculated between these values using an arbitrary internal branch
#' length.
#' 
#' @aliases crunch contrCalc
#' @param formula A model formula.
#' @param data An 'comparative.data' object. Alternatively, a data frame.
#' @param phy An object of class 'phylo', required when data is not a
#' 'comparative.data' object.
#' @param names.col A name specifying the column in 'data' that matches rows to
#' tips in 'phy', required when data is not a 'comparative.data' object.
#' @param stand.contr A logical flag indicating whether or not to standardize
#' contrasts
#' @param robust A threshold value of studentized residuals to exclude from the
#' model.
#' @param ref.var A reference variable present in the model that is used to
#' specify the direction of calculation of contrasts. If null, this is assumed
#' to be the first explanatory variable.
#' @param node.depth A positive integer greater than 1 used to restrict the
#' model to contrasts with a node depth less than or equal to the specified
#' depth. Tips have a depth of 1.
#' @param polytomy.brlen The internal branch length used for calculating
#' contrasts at a polytomy, following Pagel's (1992) method.
#' @param equal.branch.length If set to 'TRUE' then all branch lengths are set
#' to 2.
#' @param factor.action One of "abort", "warn" or "allow", describing whether
#' to stop if the formula contains a factor ("abort"), or continue after
#' converting the factor to a numeric variable, either with ("warn") or without
#' ("allow") a warning.
#' @return A object of class 'caic'.
#' @section Warning: At a polytomy, subtracting the internal branch length from
#' the real branch lengths can lead to negative branch lengths. CAIC used a
#' hard-coded internal branch length of 1 for calculating crunch contrasts at
#' polytomies. From version 2.6.9, CAIC issued a warning if this lead to
#' negative branch lengths but allowed contrast calculation to continue. In
#' contrast, the implementation in \code{crunch()} uses a default internal
#' branch length (\code{polytomy.brlen}) of 0 and also treats a negative branch
#' length in a polytomy calculation as an error. In either case, contrast
#' calculation on negative branch lengths is not a desirable outcome.
#' Duplication of CAIC results therefore requires \code{polytomy.brlen} to be
#' set to 1 and an analyis \emph{cannot} be duplicated precisely if the
#' phylogeny contains polytomies with descending branches shorter than 1. The
#' method used by \code{pic.lm} to handle polytomies avoids such problems.
#' @author David Orme
#' @seealso \code{\link{caic-class}} for 'caic' object structure and methods.
#' @references Felsenstein, J.  (1985).  Phylogenies and the comparative
#' method.  Am. Nat.  125, 1-15
#' 
#' Pagel, M. D. (1992). A method for the analysis of comparative data.  J.
#' theor. Biol. 156, 431-442.
#' 
#' Purvis, A. and Rambaut, A. (1995) Comparative analysis by independent
#' contrasts (CAIC): an Apple Macintosh application for analysing comparative
#' data.  Computer Appl. Biosciences 11, 247-251.
#' @keywords models regression
#' @examples
#' 
#' data(shorebird)
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)
#' crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data=shorebird)
#' summary(crunchMod)
#' # plot the contrasts
#' crunchTab <- caic.table(crunchMod)
#' plot(Egg.Mass ~ F.Mass, crunchTab)
#' # for the actual model diagnostics
#' par(mfrow=c(3,2))
#' caic.diagnostics(crunchMod)
#' 
crunch <- function(formula, data, phy, names.col, stand.contr = TRUE,
                   robust = Inf, ref.var = NULL, node.depth = NULL,
                   polytomy.brlen = 0, equal.branch.length = FALSE,
                   factor.action = "abort") {
    # Program Flow:
    #   1) setup - check arguments,
    #   2) use model functions to get design and response matrices,
    #      including all NA data
    #   3) feed the model matrices into a function to calculate nodal
    #      values and contrasts
    #   4) feed the returned contrast versions of the design and
    #      response matrices into lm.fit

    # TODO - return node age/height
    # TODO - allow caic to be used as a contrast calculator
    # TODO - explicit check for polytomy.brlen problems

    # CHECKS AND SETUP

    # - test to see if there is a comparative data object and if not then
    #   retrofit the remaining arguments into a comparative data object.
    if (!missing(data)) {
        if (!inherits(data, "comparative.data")) {
            if (missing(names.col)) stop("names column is missing")
            names.col <- deparse(substitute(names.col))
            data <- caicStyleArgs(
                data = data, phy = phy,
                names.col = names.col, warn.dropped = TRUE
            )
        }
    }

    # extract the data and phylogeny
    cdata <- data # in case the original is needed later
    phy <- data$phy
    data <- data$data

    # check node.depth is sensible
    if (!is.null(node.depth)) {
        if (node.depth %% 1 != 0 || node.depth < 1) {
            stop("node.depth must be a positive integer greater than 1.")
        }
    }

    # set branch lengths doesn't get evaluated if FALSE or zero
    if (as.logical(equal.branch.length)) {
        phy$edge.length <- rep(2, nrow(phy$edge))
    } else {
        if (is.null(phy$edge.length)) {
            stop(
                "The phylogeny does not contain branch lengths and crunch ",
                "has not been set to use equal branch lengths."
            )
        }
        if (any(phy$edge.length <= 0)) {
            stop(
                "The phylogeny contains either negative or zero branch ",
                "lengths and crunch has not been set to use equal branch ",
                "lengths."
            )
        }
    }

    # check for factor.action
    factor.action <- match.arg(factor.action, c("abort", "warn", "allow"))

    # CALCULATE MODEL
    # GET THE MODEL MATRIX and Model Response

    # reduce to just the variables used in the formula so they can
    # all be treated as numeric but hang on to tip labels for subsetting
    # the phylogeny down to complete tips
    data <- subset(data, select = all.vars(formula))

    # HANDLE CATEGORICAL VARIABLES:
    # find the factors - (number of levels > 0)
    varLevels <- sapply(data, function(x) length(levels(x)))
    varIsOrdered <- sapply(data, is.ordered)

    if (any(varLevels > 0)) {
        # check for unordered multi states...
        if (any(varLevels > 2 & !varIsOrdered)) {
            stop("Unordered non-binary factors included in model formula.")
        }

        # otherwise check for action on viable factors...

        if (factor.action == "abort") {
            stop(
                "The formula includes factors. Change the factor.action ",
                "argument to allow these to be fitted as numeric variables."
            )
        } else if (factor.action == "warn") {
            warning(
                "The formula includes factors, which have been treated ",
                "as continuous variables "
            )
        }

        data <- as.data.frame(lapply(data, as.numeric))
    }


    # ditch the intercept, if present
    formula <- update(formula, . ~ . - 1)

    # now we have the union of the phylogeny and data
    # get the model frame, matrix and response
    # these show the values at the tips for each term
    mf <- model.frame(formula, data, na.action = na.pass)

    # is there enough data in the model
    # TODO - think whether this check is always sufficient.
    mfComplete <- complete.cases(mf)
    if (sum(mfComplete) < 2) {
        stop("Fewer than two taxa contain complete data for this analysis")
    }

    # get the design matrix
    md <- model.matrix(formula, mf)

    # sort out the reference variable
    if (is.null(ref.var)) {
        ref.var <- colnames(md)[1] # first column in the design matrix
    } else {
        ref.var <- deparse(substitute(ref.var))
        if (is.na(match(ref.var, colnames(md)))) {
            stop("Reference variable not found in design matrix")
        }
    }

    # MODEL RESPONSE
    mr <- model.response(mf)
    # turn into a column matrix
    mr <- as.matrix(mr)
    colnames(mr) <- as.character(formula[2])
    # now that we have the model response for CAIC style contrasts
    # we can substitute the reference variable
    # for empty models (i.e. models specified as X~1)
    if (is.empty.model(formula)) ref.var <- colnames(mr)
    # add to the design matrix - this strips the assign and
    # contrast attributes so save...
    attrMD <- attributes(md)
    md <- cbind(mr, md)

    # NOW SETUP TO GET CONTRASTS AND NODAL VALUES
    # We know the tip values and have the tree
    contr <- contrCalc(md, phy, ref.var, "crunch", polytomy.brlen)

    # GET RESPONSE MATRIX
    # standardize the contrasts if required
    if (stand.contr) contr$contr <- contr$contr / sqrt(contr$var.contr)

    # FEED THE RETURNED DESIGN AND RESPONSE MATRICES INTO THE
    # MODELLING FUNCTIONS
    # assemble the data into a finished contrast object
    ContrObj <- list()
    ContrObj$contr$response <- contr$contr[, 1, drop = FALSE]

    ContrObj$contr$explanatory <- contr$contr[, -1, drop = FALSE]

    ContrObj$nodalVals$response <- contr$nodVal[, 1, drop = FALSE]
    ContrObj$nodalVals$explanatory <- contr$nodVal[, -1, drop = FALSE]
    ContrObj$contrVar <- contr$var.contr
    ContrObj$nChild <- contr$nChild
    ContrObj$nodeDepth <- contr$nodeDepth

    ## need to keep the assign and contrasts attributes from the model
    ## matrix with the contrast object in order to get anova() methods to work
    ## can't store assign permanently with explanatory contrasts because
    # validNode subsetting strips attributes
    attr(ContrObj, "assign") <- attrMD$assign
    if (!is.null(attrMD$contrasts)) {
        attr(ContrObj, "contrasts") <- attrMD$contrasts
    }

    # gather the row ids of NA nodes to drop from the model
    validNodes <- with(
        ContrObj$contr,
        complete.cases(explanatory) & complete.cases(response)
    )

    # enforce any node depth requirement
    if (!is.null(node.depth)) {
        validNodes[ContrObj$nodeDepth > node.depth] <- FALSE
    }

    # save for the user
    ContrObj$validNodes <- validNodes

    # feed the contr.model.response and contr.model.matrix
    # into lm.fit to get the model and then set up the lm object
    # - need to use lm.fit here rather than calling the model on
    #   data=contrData because any functions in the formula are now
    #   set in the column names - don't want lm to try and reinterpret
    #   them in parsing the formula.
    # - the problem then becomes how to get the model to refer to the dataset

    contrMD <- ContrObj$contr$explanatory[validNodes, , drop = FALSE]
    contrRS <- ContrObj$contr$response[validNodes, , drop = FALSE]
    attr(contrMD, "assign") <- attr(ContrObj, "assign") ## replace attributes
    if (!is.null(attr(ContrObj, "contrasts"))) {
        attr(contrMD, "contrasts") <- attr(ContrObj, "contrasts")
    }

    mod <- lm.fit(contrMD, contrRS)
    class(mod) <- "lm"

    # assemble the output
    # return fitted model and contrasts
    RET <- list(contrast.data = ContrObj, mod = mod, data = cdata)
    class(RET) <- c("caic")

    ## THIS CAN BE DONE WITH MUCH MORE FINESSE!
    # convert the ContrObj into a data frame...
    # - removed the OTT call to caic.table
    # remove invalid Nodes to keep the residual and predict lengths the same
    contrData <- with(
        ContrObj$contr,
        as.data.frame(cbind(response, explanatory))
    )
    contrData <- contrData[validNodes, , drop = FALSE]
    RET$mod$call <- substitute(lm(FORM, data = contrData), list(FORM = formula))
    RET$mod$terms <- attr(mf, "terms")

    # put the model.frame in to the lm object so that predict, etc. calls work
    RET$mod$model <- contrData
    attr(RET$mod$model, "terms") <- attr(mf, "terms")

    ## Add studentized residuals: need to use matching in case of invalid nodes
    stRes <- rstudent(mod)
    SRallNodes <- rep(NA, length(RET$contrast.data$validNodes))
    names(SRallNodes) <- names(RET$contrast.data$contrVar)
    SRallNodes[match(names(stRes), names(SRallNodes))] <- stRes
    RET$contrast.data$studentResid <- SRallNodes

    ## add some attributes
    attr(RET, "contr.method") <- "crunch"
    attr(RET, "macro.method") <- ""
    attr(RET, "stand.contr") <- stand.contr
    attr(RET, "robust") <- robust

    # lastly, test for studentised outliers
    if (any(stRes > robust)) {
        RET <- caic.robust(RET, robust)
    }

    return(RET)
}
