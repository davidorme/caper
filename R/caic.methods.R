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
#' crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data=shorebird)
#' print(crunchMod)
#' summary(crunchMod)
#' 
summary.caic <- function(object, ...) {
    summary(object$mod, ...)
}

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

logLik.caic <- function(object, ...) {
    logLik(object$mod, ...)
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
#' 
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
#' cMod1 <- crunch(lgEgg.Mass ~ lgM.Mass * lgF.Mass, data=shorebird)
#' cMod2 <- crunch(lgEgg.Mass ~ lgM.Mass + lgF.Mass, data=shorebird)
#' cMod3 <- crunch(lgEgg.Mass ~ lgM.Mass , data=shorebird)
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

plot.caic <- function(x, ...) {
    plot(x$mod, ...)
}

residuals.caic <- function(object, ...) {
    residuals(object$mod, ...)
}

coef.caic <- function(object, ...) {
    coef(object$mod, ...)
}
