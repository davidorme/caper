#' The 'caic' S3 object class and methods
#'
#' The functions 'crunch', 'brunch', 'macrocaic' and 'piclm' all return an
#' object containing an independent contrast model. The structure of the object
#' and the available methods are described here.
#'
#' @name caic
#' @format A 'caic' object is a list containing the following:
#' \describe{
#'   \item{contrast.data}{ A list of the following:
#'     \describe{
#'       \item{contr}{
#'          A list containing matrices of the contrasts in the response variable
#'          (contr\$response) and explanatory variables (contr\$explanatory).
#'       }
#'       \item{nodalVals}{
#'          A list containing matrices of the nodal values in the response
#'          variable (contr\$response) and explanatory variables
#'          (contr\$explanatory).
#'       }
#'       \item{contrVar}{
#'          A numeric vector of the expected variance for each contrast.
#'       }
#'       \item{nChild}{
#'          A vector showing the number of nodes descending from each
#'          internal node
#'       }
#'       \item{nodeDepth}{
#'          A vector showing the maximum number of nodes between each internal
#'          node and the tips of the phylogeny (including both the node in
#'          question and the tip and hence always >=2)
#'       }
#'       \item{validNodes}{
#'          A logical vector showing which internal nodes on the tree have valid
#'          contrasts, given the available data and any user constraints.
#'       }
#'     }
#'   }
#'   \item{data}{
#'     A 'comparative.data' object containing the phylogeny used to calculate
#'     contrasts and the original data.
#'   }
#'   \item{lm}{
#'     An 'lm' object containing a regression model through the origin for the
#'    calculated contrast
#'   }
#' }
#'
#' In addition, the object may have the following attributes:
#' \describe{
#'   \item{contr.method}{One of 'crunch', 'brunch' or 'piclm'.}
#'   \item{macro}{
#'     Either 'RRD' or 'PDI' if the response contrasts are calculated as species
#'     richness contrasts using \code{\link{macrocaic}}
#'   }
#'   \item{stand.cont}{
#'     A logical value showing whether the contrasts in the object have been
#'     standardized.
#'   }
#' }
#' @keywords class
#' @export
NULL



#' Summarize a crunch, brunch or macrocaic analysis
#'
#' The summary method simply returns the linear model summary from the 'caic'
#' object. The print method prints some basic information about the analysis
#' followed by the model summary.
#'
#'
#' @aliases summary.caic print.caic
#' @param object An object of class 'caic'.
#' @param x An object of class 'caic'.
#' @param \dots Arguments to be passed to 'summary.lm'.
#' @return The summary method returns an object of class 'summary.lm'.
#' @author David Orme
#' @seealso \code{\link{crunch}},\code{\link{brunch}}, \code{link{macrocaic}}
#' @keywords methods
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)
#' crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data = shorebird)
#' print(crunchMod)
#' summary(crunchMod)
#'
summary.caic <- function(object, ...) {
    summary(object$mod, ...)
}

#' @describeIn summary.caic Print CAIC model details
print.caic <- function(x, ...) {
    cat("Phylogenetic Independent Contrasts analysis using:",
        attr(x, "contr.method"), ".\n",
        sep = ""
    )
    if (!is.null(attr(x, "macro.method"))) {
        cat(
            "Response values are species rich contrasts using: ",
            attr(x, "macro.method"), "\n"
        )
    }
    cat("\nPhylogeny: ", attr(x, "phyName"),
        " (", length(x$data$phy$tip.label), " tips)\n",
        sep = ""
    )
    cat("Data: ", attr(x, "dataName"),
        " (", nrow(x$data$data), " rows)\n",
        sep = ""
    )
    cat("Number of valid contrasts: ",
        sum(x$contrast.data$validNodes), "\n",
        sep = ""
    )

    stres <- na.omit(x$contrast.data$studentResid)
    robust <- attr(x, "robust")
    if (any(abs(stres) > robust)) {
        nNonrobust <- sum(abs(stres) > robust)
        cat("Excluding ", nNonrobust,
            ifelse(nNonrobust > 1, " contrasts", " contrast"),
            " with absolute studentised residuals > ", robust, "\n",
            sep = ""
        )
    }

    print(summary(x))
}


#' Anova and model checking methods for independent contrast models.
#'
#' These functions provide ANOVA tables and model comparison using ANOVA and
#' AIC, along with standard model diagnostic plots and accessor functions for
#' phylogenetic independent contrast objects.
#'
#' The 'anova' method provides access to single anova tables for a model and to
#' comparison of lists of models. The 'logLik' method provides access to the
#' log likelihood of the 'caic' model and hence to AIC comparison of models.
#'
#' The 'plot' method uses the standard set of model diagnostic plots for linear
#' models. It is also wise to check the evolutionary assumptions of independent
#' contrast models using the 'caic' specific diagnostic plots. The 'predict'
#' and 'residuals' functions provide access to these parts of the 'caic'
#' object.
#' @name anova.caic
#' @aliases anova.caic anova.caiclist logLik.caic plot.caic predict.caic
#' residuals.caic coef.caic
#' @param object An object of class 'caic'.
#' @param scale A character string specifying the test statistic to be used.
#' Can be one of "F", "Chisq" or "Cp", with partial matching allowed, or NULL
#' for no test.
#' @param test numeric. An estimate of the noise variance sigma^2. If zero this
#' will be estimated from the largest model considered.
#' @param x An object of class 'caic'.
#' @param \dots Further argument to be passed to methods.
#' @author David Orme
#' @seealso \code{\link{crunch}},
#' \code{\link{brunch}},\code{\link{macrocaic}},\code{\link{caic.diagnostics}}
#' @keywords utils stats
#' @examples
#'
#' data(shorebird)
#' shorebird.data$lgEgg.Mass <- log(shorebird.data$Egg.Mass)
#' shorebird.data$lgM.Mass <- log(shorebird.data$M.Mass)
#' shorebird.data$lgF.Mass <- log(shorebird.data$F.Mass)
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)
#'
#' cMod1 <- crunch(lgEgg.Mass ~ lgM.Mass * lgF.Mass, data = shorebird)
#' cMod2 <- crunch(lgEgg.Mass ~ lgM.Mass + lgF.Mass, data = shorebird)
#' cMod3 <- crunch(lgEgg.Mass ~ lgM.Mass, data = shorebird)
#'
#' anova(cMod1, cMod2, cMod3)
#' AIC(cMod1, cMod2, cMod3)
#'
#' plot(cMod3)
#'
anova.caic <- function(object, ...) {
    ## borrowing from anova.lm
    if (length(list(object, ...)) == 1L) {
        # no other objects, no other args (scale and test only
        # make sense for multiple models)
        anova(object$mod)
    } else {
        # pass on - having a second function allows the easy interception
        # of test and scale arguments out of the list of objects
        return(anova.caiclist(object, ...))
    }
}

#' @describeIn anova.caic Calculate an ANOVA table for a list of CAIC models
anova.caiclist <- function(object, ..., scale = 0, test = "F") {
    ## ANOVA cares about model types - need to check crunch
    # need to check that the contrast methods are the same
    objects <- list(object, ...)

    objectsClass <- sapply(objects, class)
    if (!all(objectsClass == "caic")) {
        stop("anova() on mix of 'caic' and non-'caic' objects.")
    }

    objectsContrMethod <- sapply(objects, attr, "contr.method")
    if (length(unique(objectsContrMethod)) > 1L) {
        stop("anova() on mixed contrast methods")
    }

    objectsMacroMethod <- sapply(objects, attr, "macro.method")
    if (length(unique(objectsMacroMethod)) > 1L) {
        stop("anova() on mixed macrocaic methods")
    }

    ## OK - now pass the mod parts of those object into anova.lmlist()
    mods <- lapply(objects, "[[", "mod")
    args <- c(mods, list(scale = scale, test = test))
    anv <- do.call("anova", args)
    return(anv)
}

#' @describeIn anova.caic Extract model predictions from a CAIC model
predict.caic <- function(object, ...) {
    # need to force the model to get predictions using the contrast table
    # rather than the original data table...
    # don't completely hijack the newdata argument...

    dots <- list(...)
    newdataProv <- pmatch(names(dots), "newdata")
    if (all(is.na(newdataProv))) {
        nD <- caic.table(object)
    } else {
        nD <- dots[[newdataProv]]
    }
    predict(object$mod, newdata = nD)
}

#' @describeIn anova.caic Extract the log likelihood from a CAIC model
logLik.caic <- function(object, ...) {
    logLik(object$mod, ...)
}

#' @describeIn anova.caic Plot a CAIC model
plot.caic <- function(x, ...) {
    plot(x$mod, ...)
}

#' @describeIn anova.caic Extract residuals from a CAIC model
residuals.caic <- function(object, ...) {
    residuals(object$mod, ...)
}

#' @describeIn anova.caic Extract model coefficients from a CAIC model
coef.caic <- function(object, ...) {
    coef(object$mod, ...)
}
