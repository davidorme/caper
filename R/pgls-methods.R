#' Generic model methods for 'pgls' models.
#'
#' These are simple summary methods, accessor functions and summary and print
#' methods for 'pgls' models.
#'
#' Phylogenetically corrected residuals from 'pgls' models [TODO].
#'
#' Note that the r^2 values reported by \code{summary.pgls} have a specific
#' interpretation. \code{pgls} fits the intercept-only model for the data using
#' _exactly_ the same covariance matrix (phylogeny plugged through any branch
#' length transformations) as the fitted model to get a null model. The
#' r-squared and adjusted r-squared that are reported therefore hold the
#' covariance matrix constant, so show percentage of variance explained between
#' a null model and the actual model given that precise model of trait change.
#'
#' The actual ML null model for the data (optimising the BL transformation
#' independently) might be different from this - but then the r squared values
#' confound change in explanatory power from changing the model parameters and
#' from changing the trait model.
#'
#' @name pgls-methods
#' @aliases pgls-methods coef.pgls residuals.pgls fitted.pgls predict.pgls
#' print.pgls summary.pgls print.summary.pgls nobs.pgls
#' @param object An object of class 'pgls'.
#' @param x An object of class 'pgls'.
#' @param phylo Return phylogenetically corrected residuals or ordinary
#' residuals (see details).
#' @param newdata Alternative data for predicting from 'pgls' models.
#' @param digits Number of digits to show in summary methods.
#' @param ... Further arguments to methods.
#' @return The 'summary' method returns an object of class 'summary.pgls'
#' containing: \item{call}{The original function call creating the model.}
#' \item{df}{A vector of the degrees of freedom used to estimate parameters and
#' the residual degrees of freedom.} \item{sigma}{The square root of the
#' estimated variance of the random error.} \item{residuals}{The
#' phylogenetically corrected residuals.} \item{coefficients}{A table of model
#' coefficient, standard errors and t values.} \item{param}{A vector of branch
#' length parameters used in the model.} \item{mlVals}{A vector showing which
#' branch length parameters have been optimised.} \item{param.CI}{A list of
#' length three containing confidence intervals and p values on parameter
#' bounds for each parameter.} \item{fstatistic}{A vector of the F value,
#' numerator and denominator degrees of freedom for the model.}
#' \item{r.squared}{The r^2 for the model.} \item{adj.r.squared}{The adjusted
#' r^2 for the model.}
#' @author Rob Freckleton, David Orme
#' @seealso \code{\link{pgls}}
#' @keywords utils stats
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(
#'     shorebird.tree, shorebird.data, Species,
#'     vcv = TRUE, vcv.dim = 3
#' )
#' mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird)
#' print(mod1)
#'
#' mod1.sum <- summary(mod1)
#' print(mod1.sum)
NULL



#' @describeIn pgls-methods Extract a summary object from a pgls model
#' @export
summary.pgls <- function(object, ...) {
    ## call and return object
    ans <- list(call = object$call)
    class(ans) <- "summary.pgls"

    ## model size
    p <- object$k
    n <- object$n
    rdf <- n - p
    ans$df <- c(p, rdf)

    ## residuals and residual standard error
    r <- object$phyres
    rss <- object$RSSQ
    resvar <- rss / rdf
    ans$sigma <- sqrt(resvar)
    ans$residuals <- r

    ## coefficient matrix
    cf <- object$model$coef
    se <- object$sterr
    t <- cf / se

    coef <- cbind(cf, se, t, 2 * (1 - pt(abs(t), rdf)))
    colnames(coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    ans$coefficients <- coef

    ## parameter matrix
    ans$param <- object$param
    ans$mlVals <- object$mlVals

    if (!is.null(object$param.CI)) ans$param.CI <- object$param.CI

    if (!is.null(object$na.action)) ans$na.action <- object$na.action

    # model statistics: p includes the intercept - it is the number of columns
    # of the design matrix
    ans$fstatistic <- c(
        value = ((object$NSSQ - object$RSSQ) / object$RMS) / (object$k - 1),
        numdf = p - 1,
        dendf = rdf
    )
    ans$r.squared <- (object$NSSQ - object$RSSQ) / object$NSSQ
    ans$adj.r.squared <- (object$NMS - object$RMS) / object$NMS

    return(ans)
}

#' @describeIn pgls-methods Print a pgls model summary object.
#' @export
print.summary.pgls <- function(x, digits = max(3, getOption("digits") - 3),
                               ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n",
        sep = ""
    )

    r <- zapsmall(quantile(x$resid), digits + 1)
    names(r) <- c("Min", "1Q", "Median", "3Q", "Max")
    cat("Residuals:\n")
    print(r, digits = digits)


    cat("\nBranch length transformations:\n\n")
    for (p in names(x$param)) {
        cat(sprintf(
            "%-6s [%s]  : %0.3f\n", p,
            ifelse(x$mlVals[p], " ML", "Fix"), x$param[p]
        ))
        if (!is.null(x$param.CI[[p]])) {
            blopt <- x$param.CI[[p]]
            cat(sprintf(
                "   lower bound : %0.3f, p = %-5s\n",
                blopt$bounds.val[1], format.pval(blopt$bounds.p[1])
            ))
            cat(sprintf(
                "   upper bound : %0.3f, p = %-5s\n",
                blopt$bounds.val[2], format.pval(blopt$bounds.p[2])
            ))
            cat(sprintf(
                "   %2.1f%% CI   : (%0.3f, %0.3f)\n",
                blopt$ci * 100, blopt$ci.val[1], blopt$ci.val[2]
            ))
        }
    }


    cat("\nCoefficients:\n")
    printCoefmat(x$coef)

    cat("\nResidual standard error:", format(signif(
        x$sigma,
        digits
    )), "on", x$df[2L], "degrees of freedom\n")
    if (nzchar(mess <- naprint(x$na.action))) {
        cat("  (", mess, ")\n", sep = "")
    }
    cat("Multiple R-squared:", formatC(x$r.squared, digits = digits))
    cat(
        ",\tAdjusted R-squared:", formatC(x$adj.r.squared,
            digits = digits
        ), "\nF-statistic:", formatC(x$fstatistic[1L],
            digits = digits
        ), "on", x$fstatistic[2L], "and",
        x$fstatistic[3L], "DF,  p-value:", format.pval(
            pf(x$fstatistic[1L],
                x$fstatistic[2L], x$fstatistic[3L],
                lower.tail = FALSE
            ),
            digits = digits
        ), "\n"
    )
}

#' @describeIn pgls-methods Print a pgls model
#' @export
print.pgls <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n",
        sep = ""
    )

    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
        print.gap = 2,
        quote = FALSE
    )
    cat("\n")
}

#' @describeIn pgls-methods Extract model coefficients from a pgls model
#' @export
coef.pgls <- function(object, ...) {
    cf <- object$model$coef
    nm <- rownames(cf)
    cf <- structure(as.vector(cf), names = nm)
    return(cf)
}

#' @describeIn pgls Extract residuals from a pgls model
#' @export
residuals.pgls <- function(object, phylo = FALSE, ...) {
    ret <- NULL
    if (phylo == FALSE) {
        ret <- object$res
    } else {
        ret <- object$phyres
    }
    return(ret)
}

#' @describeIn pgls-methods Extract fitted values from a pgls model
#' @export
fitted.pgls <- function(object, ...) {
    ret <- object$fitted
    return(ret)
}

#' @describeIn pgls-methods Extract predicted values from a pgls model
#' @export
predict.pgls <- function(object, newdata = NULL, ...) {
    # pull the data from the model if no new data is provided
    if (is.null(newdata)) {
        newdata <- object$data$data
    }

    # turn that into a design matrix
    # need to drop the response from the formula
    dmat <- model.matrix(delete.response(terms(formula(object))),
        data = newdata
    )

    # multiply through by the coefficients
    mu <- as.matrix(coef(object))
    ret <- dmat %*% mu
    return(ret)
}

#' @describeIn pgls-methods Extract the log likelihood from a pgls model
#' @export
logLik.pgls <- function(object, REML = FALSE, ...) {
    val <- object$model$log.lik

    attr(val, "nall") <- object$n
    attr(val, "nobs") <- object$n
    attr(val, "df") <- object$k
    class(val) <- "logLik"
    val
}

#' @describeIn pgls-methods Extract the number of observations from a pgls model
#' @export
nobs.pgls <- function(object, ...) length(resid(object))



#' Diagnostic plots for 'pgls' models.
#'
#' The function generates four diagnostics plots for 'pgls' models.
#'
#' The first two plots show the fit of the phylogenetic residuals from the
#' model to a normal distribution: a density plot of the residuals and a normal
#' Q-Q plot. The second two plots scatterplots show pattern in the distribution
#' of the fitted values against the observed and residual values.
#'
#' @param x An object of class 'pgls'.
#' @param ... Additional arguments to plot functions.
#' @author Rob Freckleton, David Orme
#' @seealso \code{\link{pgls}}
#' @keywords utils graphics
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(
#'     shorebird.tree, shorebird.data, Species,
#'     vcv = TRUE, vcv.dim = 3
#' )
#' mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird)
#' par(mfrow = c(2, 2))
#' plot(mod1)
#'
plot.pgls <- function(x, ...) {
    res <- residuals(x, phylo = TRUE)
    res <- res / sqrt(var(res))[1]
    plot(density(res))
    qqnorm(res)
    qqline(res)
    plot(fitted(x), res,
        xlab = "Fitted value",
        ylab = "Residual value (corrected for phylogeny)"
    )
    plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value")
}



#' Anova and AIC tables for 'pgls' models.
#'
#' The 'anova' function creates ANOVA tables for a 'pgls' models using
#' sequential sums of squares.
#'
#' The sequential sums of squares are calculated by refitting the model in the
#' order of the terms of the formula and so can take a little time to
#' calculate. Branch length transformations are held at the values of the
#' initial object. The 'logLik.pgls' provides a simple accessor function that
#' allows the use of AIC model comparisons. Note that the generic AIC methods
#' do no checking to ensure that sensible models are being compared.
#'
#' @aliases anova.pgls anova.pglslist logLik.pgls
#' @param object A 'pgls' model object.
#' @param \dots Additional 'pgls' models.
#' @param scale A character string specifying the test statistic to be used.
#' Can be one of "F", "Chisq" or "Cp", with partial matching allowed, or NULL
#' for no test.
#' @param test numeric. An estimate of the noise variance sigma^2. If zero this
#' will be estimated from the largest model considered.
#' @return A table of class 'anova' and 'data.frame' that employs the generic
#' plot methods for 'anova' tables.
#' @note The functions build heavily on the generic methods 'anova.lm' and
#' 'anova.lmlist'.
#' @author Rob Freckleton, David Orme
#' @seealso \code{\link{pgls}}
#' @keywords utils stats
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(
#'     shorebird.tree, shorebird.data, Species,
#'     vcv = TRUE, vcv.dim = 3
#' )
#'
#' mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird)
#' anova(mod1)
#'
#' mod2 <- pgls(log(Egg.Mass) ~ log(M.Mass) + log(F.Mass), shorebird)
#' mod3 <- pgls(log(Egg.Mass) ~ log(M.Mass), shorebird)
#' mod4 <- pgls(log(Egg.Mass) ~ 1, shorebird)
#'
#' anova(mod1, mod2, mod3, mod4)
#' AIC(mod1, mod2, mod3, mod4)
#'
anova.pgls <- function(object, ...) {
    ## SEQUENTIAL SUMS OF SQUARES.
    ## ASSUMES ORDER OF TERMS PRESERVE MARGINALITY

    if (length(list(object, ...)) > 1L) {
        return(anova.pglslist(object, ...))
    } else {
        data <- object$data
        tlabels <- attr(terms(object$formula), "term.labels")
        k <- object$k
        n <- object$n
        NR <- length(tlabels) + 1

        # track residual ss and residual df and get residuals
        # and df of null model
        rss <- resdf <- rep(NA, NR)
        rss[1] <- object$NSSQ
        resdf[1] <- n - 1

        lm <- object$param["lambda"]
        dl <- object$param["delta"]
        kp <- object$param["kappa"]

        # fit the sequential models
        for (i in seq_along(tlabels)) {
            fmla <- as.formula(
                paste(object$namey, " ~ ", paste(tlabels[1:i], collapse = "+"))
            )
            plm <- pgls(fmla, data, lambda = lm, delta = dl, kappa = kp)
            rss[i + 1] <- plm$RSSQ
            resdf[i + 1] <- (n - 1) - plm$k + 1
        }

        ss <- c(abs(diff(rss)), object$RSSQ)
        df <- c(abs(diff(resdf)), n - k)
        ms <- ss / df
        fval <- ms / ms[NR]
        P <- pf(fval, df, df[NR], lower.tail = FALSE)

        table <- data.frame(df, ss, ms, f = fval, P)
        table[length(P), 4:5] <- NA
        dimnames(table) <- list(c(tlabels, "Residuals"), c(
            "Df",
            "Sum Sq", "Mean Sq", "F value", "Pr(>F)"
        ))
        structure(table,
            heading = c(
                "Analysis of Variance Table",
                sprintf(
                    paste0(
                        "Sequential SS for pgls: lambda = %0.2f, ",
                        "delta = %0.2f, kappa = %0.2f\n"
                    ), lm, dl, kp
                ),
                paste("Response:", deparse(formula(object)[[2L]]))
            ),
            class = c("anova", "data.frame")
        )
    }
}

#' @describeIn anova.pgls Calculate an ANOVA table for a list of pgls models.
#' @export
anova.pglslist <- function(object, ..., scale = 0, test = "F") {
    objects <- list(object, ...)

    ## check the models use the same response
    responses <- as.character(
        lapply(objects, function(x) deparse(terms(x$formula)[[2L]]))
    )
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
        objects <- objects[sameresp]
        warning(
            "models with response ", deparse(responses[!sameresp]),
            " removed because response differs from ", "model 1"
        )
    }

    ## check the models have the same number of cases (not actually that they
    # are the same values)
    ns <- sapply(objects, function(x) length(x$residuals))
    if (any(ns != ns[1L])) {
        stop("models were not all fitted to the same size of dataset")
    }

    ## check that the model parameters are the same
    param <- sapply(objects, "[[", "param")
    paramChk <- apply(param, 1, function(X) all(X == X[1]))
    if (!all(paramChk)) {
        stop("models were fitted with different branch length transformations.")
    }

    nmodels <- length(objects)
    if (nmodels == 1) {
        return(anova(object))
    }
    resdf <- as.numeric(lapply(objects, function(X) X$n - X$k))
    resdev <- as.numeric(lapply(objects, "[[", "RSSQ"))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(
        NA,
        -diff(resdev)
    ))
    variables <- lapply(objects, function(x) {
        paste(deparse(formula(x)),
            collapse = "\n"
        )
    })
    dimnames(table) <- list(1L:nmodels, c(
        "Res.Df", "RSS", "Df",
        "Sum of Sq"
    ))
    title <- "Analysis of Variance Table"
    subtitle <- sprintf(
        "pgls: lambda = %0.2f, delta = %0.2f, kappa = %0.2f\n",
        param["lambda", 1], param["delta", 1], param["kappa", 1]
    )
    topnote <- paste("Model ", format(1L:nmodels), ": ", variables,
        sep = "", collapse = "\n"
    )
    if (!is.null(test)) {
        bigmodel <- order(resdf)[1L]
        scale <- if (scale > 0) {
            scale
        } else {
            resdev[bigmodel] / resdf[bigmodel]
        }
        table <- stat.anova(
            table = table, test = test, scale = scale,
            df.scale = resdf[bigmodel], n = length(objects[bigmodel$residuals])
        )
    }
    structure(table, heading = c(title, subtitle, topnote), class = c(
        "anova",
        "data.frame"
    ))
}
