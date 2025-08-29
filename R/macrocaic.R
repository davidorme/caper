#' Comparative analysis using independent contrasts on species richness data.
#'
#' Macroevolutionary hypotheses about correlates of species richness require
#' testing in a phylogenetic framework in order to avoid phylogenetic
#' autocorrelation. Independent contrasts as described by Felsenstein (1985)
#' are appropriate for explanatory variables in such models but not for species
#' richness as the response variable. This function implements two methods for
#' calculating species richness constrasts described by Agapow and Isaac (2002)
#' and originally implemented in the program MacroCAIC.
#'
#' The 'macrocaic' function fits a regression to the formula provided using
#' 'crunch' contrasts for continuous explanatory variables and species richness
#' contrasts for the response. The species richness contrasts are either the
#' relative rate difference (RRD) or proportion dominance index (PDI):
#'
#' \deqn{RRD = \ln\left(\frac{N_1}{N_2}\right)}
#'
#' \deqn{PDI = \left(\frac{N_1}{N_1+N_2}\right)-0.5}
#'
#' The values \eqn{N_1} and \eqn{N_2} are the species richness of the two
#' daughter nodes and \eqn{N_1} is the species richness of the clade with the
#' larger value of the reference variable. Species richness contrasts are not
#' calculated at polytomies. Nodal values for species richness are calculated
#' as the sum of the richness of the daughter nodes.
#'
#' @param formula A formula describing a linear model predicting species
#' richness.
#' @param data A data frame containing the variables to be used in the model.
#' @param phy An object of class 'phylo'.
#' @param names.col A name identifying a column in \code{data} that contains
#' the tip labels from \code{phy}.
#' @param macroMethod One of either "RRD" or "PDI" (see Details).
#' @param stand.contr A logical flag indicating whether to standardize the
#' contrasts.
#' @param robust A threshold value of studentized residuals to exclude from the
#' model.
#' @param ref.var Identifies a predictor variable used for determining the
#' direction of contrasts.
#' @param node.depth A positive integer greater than 1 used to restrict the
#' model to contrasts with a node depth less than or equal to the specified
#' depth. Tips have a depth of 1.
#' @param macroMinSize A positive integer giving the minimum species richness
#' at a node for contrasts to be included in the model.
#' @param equal.branch.length If set to 'TRUE' then all branch lengths are set
#' to 2.
#' @return A object of class 'caic'.
#' @author David Orme
#' @seealso \code{\link{caic-class}} for 'caic' object structure and methods.
#' @references Felsenstein, J.  (1985).  Phylogenies and the comparative
#' method.  Am. Nat.  125, 1-15 Agapow, P.-M. and Isaac, N. J. B. (2002)
#' MacroCAIC: correlates of species richness. Diversity & Distributions, 8,
#' 41-43 Isaac, N., Agapow, P., Harvey, P., and Purvis, A. (2003).
#' Phylogenetically nested com- parisons for testing correlates of species
#' richness: A simulation study of continuous variables. Evolution,
#' 57(1):18-26.
#' @keywords models regression
#' @examples
#'
#' data(IsaacEtAl)
#' primates <- comparative.data(
#'     primates.tree, primates.data, binomial,
#'     na.omit = FALSE
#' )
#' primatesBodySize <- macrocaic(species.rich ~ body.mass, data = primates)
#' summary(primatesBodySize)
#'
#' @export
macrocaic <- function(formula, data, phy, names.col, macroMethod = "RRD",
                      stand.contr = TRUE, robust = Inf, ref.var = NULL,
                      node.depth = NULL, macroMinSize = 3,
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

    # check branch lengths
    if (as.logical(equal.branch.length)) {
        # doesn't get evaluated if FALSE or zero
        phy$edge.length <- rep(2, nrow(phy$edge))
    } else {
        if (is.null(phy$edge.length)) {
            stop(
                "The phylogeny does not contain branch lengths and ",
                "macrocaic has not been set to use equal branch lengths."
            )
        }
        if (any(phy$edge.length < 0)) {
            stop(
                "The phylogeny contains negative branch lengths and ",
                "macrocaic has not been set to use equal branch lengths."
            )
        }
    }

    # MACROCAIC SPECIFIC
    # set intermediate branch length to use at polytomies - CAIC used 1, whereas
    # MacroCAIC requires 0 in order to give subnode contrast calculation that is
    # equivalent to a weighted mean
    crunch.brlen <- 0

    # GET THE BASIC MODEL MATRICES

    # drop any intercept from the formula, which has no effect if the interecept
    # is already omitted
    formula <- stats::update(formula, . ~ . - 1)
    if (stats::is.empty.model(formula)) {
        stop(
            "Macrocaic requires an explanatory variable to determine the ",
            "direction of species richness contrasts.\n",
            "Models of the form nSpp ~ 1 are not meaningful."
        )
    }

    # Get the model frame including missing data
    # and check the number of complete cases in the model frame
    initMf <- stats::model.frame(formula, data, na.action = "stats::na.pass")
    initMfComplete <- stats::complete.cases(initMf)
    # TODO - think whether this check is always sufficient...
    if (sum(initMfComplete) < 2) {
        stop("Fewer than two taxa contain complete data for this analysis")
    }

    # macro analyses with missing data on species richness are a bad thing
    macroMf <- as.matrix(stats::model.response(initMf))
    colnames(macroMf) <- with(
        attributes(attr(initMf, "terms")), rownames(factors)[response]
    )
    if (any(is.na(macroMf))) {
        stop("MacroCAIC analyses cannot have missing species richness values")
    }
    if (any(macroMf <= 0)) {
        stop("Species richness values cannot be negative or zero")
    }
    if (any((macroMf %% 1) > 0)) {
        stop("Non-integer species richness values present")
    }

    # CALCULATE MODEL
    # GET THE MODEL MATRIX and Model Response

    # get the model frame, matrix and response
    # these show the values at the tips for each term
    mf <- stats::model.frame(formula, data, na.action = stats::na.pass)

    # macroCAIC is not intended to handle categorical data, but could be used
    #  to control for ordered factors where intervals are feasible equal

    # HANDLE CATEGORICAL VARIABLES:
    # find the factors
    varClass <- attributes(attributes(mf)$terms)$dataClasses
    termFactors <- attributes(attributes(mf)$terms)$factors
    factorCols <- names(varClass)[varClass %in% c("ordered", "factor")]

    if (any(varClass %in% c("ordered", "factor") & rowSums(termFactors) > 1)) {
        stop(
            "Interactions using categorical variables not ",
            "supported in macrocaic analyses"
        )
    }
    termClass <- apply(
        termFactors, 2,
        function(X) unique(varClass[as.logical(X)])
    )

    for (fact in factorCols) {
        # - check whether all factors are ordered or binary
        currFact <- with(mf, get(fact))
        lev <- levels(currFact)
        ord <- is.ordered(currFact)
        if (length(lev) > 2 && !ord) {
            stop("Unordered non-binary factor included in model formula.")
        }
        # - modify the model frame object to make the factors numeric
        # - quote the names of the variables to assign to
        eval(parse(
            text = paste("mf$'", fact, "'<- as.numeric(currFact)", sep = "")
        ))
        attr(mf, "dataClasses") <- rep("numeric", dim(termFactors)[2])
    }

    # MODEL RESPONSE
    mr <- stats::model.response(mf)
    # turn into a column matrix
    mr <- as.matrix(mr)
    colnames(mr) <- as.character(formula[2])

    # get the design matrix
    md <- stats::model.matrix(formula, mf)

    # sort out the reference variable
    if (!is.null(ref.var)) {
        ref.var <- deparse(substitute(ref.var))
        if (is.na(match(ref.var, colnames(md)))) {
            stop("Reference variable not found in the model predictors")
        }
    } else {
        ref.var <- colnames(md)[1]
    }

    # add to the design matrix - this strips the assign and contrast attributes
    # so save...
    attrMD <- attributes(md)
    md <- cbind(mr, md)

    # NOW GET CONTRASTS AND NODAL VALUES

    contr <- contrCalc(md, phy, ref.var,
        picMethod = "crunch",
        crunch.brlen, macro = macroMethod
    )

    # the first column has been calculated as a macro contrast
    mrC <- contr$contr[, 1, drop = FALSE]
    mdC <- contr$contr[, -1, drop = FALSE]

    # standardize the contrasts if required
    # (don't standardize macro contrasts or categorical variables in Brunch)

    if (stand.contr) {
        notCateg <- !termClass %in% c("factor", "ordered")
        mdC[, notCateg] <- mdC[, notCateg, drop = FALSE] / sqrt(contr$var.contr)
    }

    # FEED THE RETURNED DESIGN AND RESPONSE MATRICES INTO THE MODELLING
    # FUNCTIONS assemble the data into a finished contrast object

    ContrObj <- list()
    ContrObj$contr$response <- mrC
    ContrObj$contr$explanatory <- mdC

    ContrObj$nodalVals$response <- contr$nodVal[, 1, drop = FALSE]
    ContrObj$nodalVals$explanatory <- contr$nodVal[, -1, drop = FALSE]

    ContrObj$contrVar <- contr$var.contr
    ContrObj$nChild <- contr$nChild

    # need to keep the assign and contrasts attributes from the model matrix
    # with the contrast object in order to get anova() methods to work can't
    # store assign permanently with explanatory contrasts because validNode
    # subsetting strips attributes
    attr(ContrObj, "assign") <- attrMD$assign
    if (!is.null(attrMD$contrasts)) {
        attr(ContrObj, "contrasts") <- attrMD$contrasts
    }

    # gather the row ids of NA nodes to drop from the model (missing data plus
    # polytomies in macro analyses)
    validNodes <- with(
        ContrObj$contr,
        stats::complete.cases(explanatory) & stats::complete.cases(response)
    )

    # macrocaic can exclude species poor nodes (default is to exclude sister
    # species pairs)
    validNodes[which(ContrObj$nodalVal$response < macroMinSize)] <- FALSE

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

    ## need to pass the assign and contrasts attributes over from the model
    ## matrix in order to get anova() methods to work
    contrMD <- ContrObj$contr$explanatory[validNodes, , drop = FALSE]
    contrRS <- ContrObj$contr$response[validNodes, , drop = FALSE]

    attr(contrMD, "assign") <- attr(ContrObj, "assign") ## replace attributes
    if (!is.null(attr(ContrObj, "contrasts"))) {
        attr(contrMD, "contrasts") <- attr(ContrObj, "contrasts")
    }

    mod <- with(ContrObj$contr, stats::lm.fit(contrMD, contrRS))
    class(mod) <- "lm"

    # assemble the output
    # return fitted model and contrasts
    RET <- list(contrast.data = ContrObj, data = cdata, mod = mod)
    class(RET) <- c("caic")

    # convert the ContrObj into a data frame...
    contrData <- with(
        ContrObj$contr,
        as.data.frame(cbind(response, explanatory))
    )
    contrData <- contrData[validNodes, , drop = FALSE]
    RET$mod$call <- substitute(stats::lm(FORM, data = contrData), list(FORM = formula))
    RET$mod$terms <- attr(mf, "terms")

    # put the model.frame in to the lm object so that predict, etc. calls work
    RET$mod$model <- contrData
    attr(RET$mod$model, "terms") <- attr(mf, "terms")

    ## Add studentized residuals: need to use matching in case of invalid nodes
    stRes <- stats::rstudent(mod)
    SRallNodes <- rep(NA, length(RET$contrast.data$validNodes))
    names(SRallNodes) <- names(RET$contrast.data$contrVar)
    SRallNodes[match(names(stRes), names(SRallNodes))] <- stRes
    RET$contrast.data$studentResid <- SRallNodes

    # add some attributes
    attr(RET, "contr.method") <- "crunch"
    attr(RET, "macro.method") <- macroMethod
    attr(RET, "stand.contr") <- stand.contr
    attr(RET, "robust") <- robust

    # lastly, test for studentised outliers
    if (any(stRes > robust)) {
        RET <- caic.robust(RET, robust)
    }

    return(RET)
}
