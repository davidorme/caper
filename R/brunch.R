#' Comparative analysis using the brunch algorithm.
#'
#' Calculate a linear model using the brunch algorithm.
#'
#' This function implements the 'brunch' algorithm for modelling the
#' relationship between variables that are phylogenetically non-independent.
#' This method was described and previously implemented in the Mac Classic
#' computer programs CAIC, written by Purvis and Rambaut (1995) and updated by
#' Nick Isaac and Paul-Michael Agapow.
#'
#' The 'brunch' algorithm calculates contrasts for models that include binary
#' categorical variables. Contrasts are identified and calculated for all
#' variables in the model for a set of nodes where each side can be
#' unequivocally attributed to one or other of the categories. Unlike 'crunch',
#' nested contrasts are not calculated and each row of data at the tips is used
#' only once. This follows Burt (1989): contrasts whose paths do not meet or
#' cross at any point will be phylogenetically independent.
#'
#' Factors with more than two levels are supported but *must* be ordered to
#' allow sensible contrasts to be drawn. In addition, there is no single best
#' compromise set of contrasts with non-binary factors and implementations may
#' differ in the set chosen.
#'
#' The user provides a comparative dataset. The formula specifies the model to
#' be fitted and contrasts are calculated in those variables. The specified
#' reference variable is used to ensure that contrasts for multivariate models
#' are calculated in a consistent direction at each node. The function
#' \code{brunch} acts as a data preparation wrapper for the function
#' \code{contrCalc}, which is not intended to be directly called by users.
#' Missing data can be present in the explanatory variables: the algorithm
#' makes use of the complete data available at each node as was the case with
#' CAIC.
#'
#' Polytomies - more detail here The Mac Classic program CAIC used 1 for both
#' 'Brunch' and 'Crunch' analyses and this the default.
#'
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
#' @param equal.branch.length If set to 'TRUE' then all branch lengths are set
#' to 2.
#' @return A object of class 'caic'.
#' @author David Orme
#' @seealso \code{\link{caic-class}} for 'caic' object structure and methods.
#' @references Purvis, A. and Rambaut, A. (1995) Comparative analysis by
#' independent contrasts (CAIC): an Apple Macintosh application for analysing
#' comparative data.  Computer Appl. Biosciences 11, 247-251.
#'
#' Burt, A. (1989). Comparative methods using phylogenetically independent
#' contrasts. Oxford Surveys in Evolutionary Biology, 6:33-53.
#' @keywords models regression
#' @examples
#'
#' data(perissodactyla)
#' perisso <- comparative.data(perissodactyla.tree, perissodactyla.data, Binomial)
#' brunchMod <- brunch(log.female.wt ~ Territoriality, data = perisso)
#' summary(brunchMod)
#'
#' # plot the contrasts
#' brunchTab <- caic.table(brunchMod)
#' plot(log.female.wt ~ Territoriality, brunchTab)
#'
#' # for the actual model diagnostics
#' par(mfrow = c(3, 1))
#' caic.diagnostics(brunchMod)
#' @export
#'
brunch <- function(formula, data, phy, names.col, stand.contr = TRUE,
                   robust = Inf, ref.var = NULL, node.depth = NULL,
                   equal.branch.length = FALSE) {
    # Program Flow:
    #   1) setup - check arguments,
    #   2) use model functions to get design and response matrices, including
    #      all NA data
    #   3) feed the model matrices into a function to calculate nodal values
    #      and contrasts
    #   4) feed the returned contrast versions of the design and response
    #      matrices into lm.fit

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
                data = data, phy = phy, names.col = names.col,
                warn.dropped = TRUE
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
                "The phylogeny does not contain branch lengths and brunch ",
                "has not been set to use equal branch lengths."
            )
        }
        if (any(phy$edge.length <= 0)) {
            stop(
                "The phylogeny contains either negative or zero branch ",
                "lengths and brunch has not been set to use equal branch ",
                "lengths."
            )
        }
    }

    # CALCULATE MODEL
    # GET THE MODEL MATRIX and Model Response

    # ditch the intercept, if present
    formula <- update(formula, . ~ . - 1)

    # now we have the union of the phylogeny and data
    # get the model frame, matrix and response
    # these show the values at the tips for each term
    mf <- model.frame(formula, data, na.action = na.pass)

    # TODO - think whether this check is always sufficient...
    mfComplete <- complete.cases(mf)
    if (sum(mfComplete) < 2) {
        stop("Fewer than two taxa contain complete data for this analysis")
    }


    # HANDLE CATEGORICAL VARIABLES:

    # currently, want to exclude the possibility of interactions in a factor
    # in Brunch because I don't have a clue how that should be handled
    # and also need a vector showing which terms are categorical and numeric
    # in order to allow correct standardization of contrasts
    # - do this step before turning the factors into numbers

    varClass <- attributes(attributes(mf)$terms)$dataClasses
    termFactors <- attributes(attributes(mf)$terms)$factors

    if (any(varClass %in% c("ordered", "factor") & rowSums(termFactors) > 1)) {
        stop(
            "Interactions using categorical variables not supported ",
            "in brunch analyses"
        )
    }

    termClass <- apply(
        termFactors, 2, function(X) unique(varClass[as.logical(X)])
    )

    # now modify the variables to be numeric for calculation
    varLevels <- sapply(mf, function(x) length(levels(x)))
    varIsOrdered <- sapply(mf, is.ordered)

    # check for unordered multi states...
    if (any(varLevels > 2 & !varIsOrdered)) {
        stop("Unordered non-binary factors included in model formula.")
    }

    # refit the model frame with numericized data
    data <- as.data.frame(lapply(mf, as.numeric))
    mf <- model.frame(formula, data, na.action = na.pass)

    # get the design matrix
    md <- model.matrix(formula, mf)

    # sort out the reference variable and check for multiple or no factors...
    if (sum(varLevels > 0) == 0) {
        stop("No factors specified in the model formula")
    }
    if (sum(varLevels > 0) > 1) {
        warning("Multiple factors specified in the model formula")
    }

    if (is.null(ref.var)) {
        # the (first) factor in the design matrix
        ref.var <- colnames(data)[which(varLevels > 0)[1]]
    } else {
        ref.var <- deparse(substitute(ref.var))
        if (is.na(match(ref.var, colnames(md)))) {
            stop("Reference variable not found in design matrix")
        }
        if (!ref.var %in% colnames(data)[which(varLevels == 0)]) {
            stop("The reference variable is not a factor")
        }
    }

    # MODEL RESPONSE
    mr <- model.response(mf)
    # turn into a column matrix
    mr <- as.matrix(mr)
    colnames(mr) <- as.character(formula[2])
    # now that we have the model response for CAIC style contrasts we can
    # substitute the reference variable for empty models
    # (i.e. models specified as X~1)
    if (is.empty.model(formula)) ref.var <- colnames(mr)

    # add to the design matrix - this strips the assign and contrast
    # attributes so save...
    attrMD <- attributes(md)
    md <- cbind(mr, md)

    # NOW SETUP TO GET CONTRASTS AND NODAL VALUES
    # We know the tip values, the analysis tree
    # Note that brunch used an internal branch length of 0
    contr <- contrCalc(md, phy, ref.var, "brunch", 0)

    # GET RESPONSE MATRIX
    # first column of contrasts is response
    mrC <- contr$contr[, 1, drop = FALSE]
    mdC <- contr$contr[, -1, drop = FALSE]

    # standardize the contrasts (but not the categorical) if required

    if (stand.contr) {
        notCateg <- !termClass %in% c("factor", "ordered")
        mdC[, notCateg] <- mdC[, notCateg, drop = FALSE] / sqrt(contr$var.contr)
        mrC <- mrC / sqrt(contr$var.contr)
    }

    # FEED THE RETURNED DESIGN AND RESPONSE MATRICES INTO THE
    # MODELLING FUNCTIONS
    # assemble the data into a finished contrast object

    ContrObj <- list()
    ContrObj$contr$response <- mrC
    ContrObj$contr$explanatory <- mdC
    ContrObj$nodalVals$response <- contr$nodVal[, 1, drop = FALSE]
    ContrObj$nodalVals$explanatory <- contr$nodVal[, -1, drop = FALSE]
    ContrObj$contrVar <- contr$var.contr
    ContrObj$nChild <- contr$nChild
    ContrObj$nodeDepth <- contr$nodeDepth

    # need to keep the assign and contrasts attributes from the model
    # matrix with the contrast object in order to get anova() methods to work
    # can't store assign permanently with explanatory contrasts because
    # validNode subsetting strips attributes
    attr(ContrObj, "assign") <- attrMD$assign
    if (!is.null(attrMD$contrasts)) {
        attr(ContrObj, "contrasts") <- attrMD$contrasts
    }

    # gather the row ids of NA nodes to drop from the model
    validNodes <- with(
        ContrObj$contr, complete.cases(explanatory) & complete.cases(response)
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

    # need to pass the assign and contrasts attributes over from the model
    # matrix in order to get anova() methods to work
    contrMD <- ContrObj$contr$explanatory[validNodes, , drop = FALSE]
    contrRS <- ContrObj$contr$response[validNodes, , drop = FALSE]

    # replace attributes
    attr(contrMD, "assign") <- attr(ContrObj, "assign")
    if (!is.null(attr(ContrObj, "contrasts"))) {
        attr(contrMD, "contrasts") <- attr(ContrObj, "contrasts")
    }

    mod <- with(ContrObj$contr, lm.fit(contrMD, contrRS))
    class(mod) <- "lm"

    # assemble the output
    # return fitted model and contrasts
    RET <- list(contrast.data = ContrObj, mod = mod, data = cdata)
    class(RET) <- c("caic")

    # convert the ContrObj into a data frame...
    contrData <- with(
        ContrObj$contr, as.data.frame(cbind(response, explanatory))
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
    attr(RET, "contr.method") <- "brunch"
    attr(RET, "macro.method") <- ""
    attr(RET, "stand.contr") <- stand.contr
    attr(RET, "robust") <- robust

    # lastly, test for studentised outliers
    if (any(stRes > robust)) {
        RET <- caic.robust(RET, robust)
    }

    return(RET)
}
