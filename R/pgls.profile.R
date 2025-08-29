#' Likelihood profiles and confidence intervals for 'pgls' models.
#'
#' These functions create likelihood profiles for branch length transformations
#' in phylogenetic generalised least squares models and fit confidence
#' intervals to estimated branch length parameters.
#'
#' The 'pgls.profile' function calculates the likelihood of a 'pgls' model
#' under different values of branch length transformations. A single parameter
#' is chosen from 'lambda', 'kappa' or 'delta' to be profiled and the model
#' likelihood is calculated at 'N' equally spaced points between the parameter
#' bounds used in the model. If the model contains a maximum likelihood
#' estimate of the parameter (or if param.CI is not null) then the resulting
#' 'pgls.profile' object will contain estimated confidence intervals.
#'
#' Only one parameter is profiled at a time and the other branch length
#' parameters will be held at the fixed or ML estimates used to fit the model.
#' The 'pgls.confint' function is used by either 'pgls' or 'pgls.profile' to
#' find confidence intervals around a maximum likelihood estimate of a given
#' branch length. The model must contain an ML estimate of the parameter for
#' confidence intervals to be calculated.
#'
#' The plot method simply draws an annotated profile plot, showing the location
#' of the ML estimate and confidence intervals if present.
#'
#' @aliases pgls.profile plot.pgls.profile pgls.confint
#' @param pgls A \code{pgls} object.
#' @param which A choice of which branch length transformation ('lambda',
#' 'kappa' or 'delta') to use.
#' @param N The number of points used to profile the likelihood
#' @param param.CI A p value used to add confidence intervals to a likelihood
#' profile for a parameter.
#' @param x A 'pgls.profile' object to plot.
#' @param ... Further arguments to plot functions.
#' @return The 'pgls.profile' function returns a list containing:
#'
#' \item{x}{Parameter values at which the likelihood has been calculated.}
#' \item{logLik}{The likelihood value at each value.}
#' \item{which}{The parameter being profiled.}
#' \item{pars}{The value of the other fixed parameters.}
#' \item{dname}{
#'      The name of the 'comparative.data' object used to fit the model.
#' }
#' \item{formula}{The formula of the model being profiled}
#'
#' If the model contains an ML estimate of the parameter being profiled, then
#' the 'pgls.profile' object will also contain the output of 'pgls.confint':
#'
#' \item{opt}{The maximum likelihood value of the parameter.}
#' \item{bounds.val}{The values of the bounds on the parameter.}
#' \item{bounds.p}{
#'      The p value of the likelihood at the bounds, given the ML value.
#' }
#' \item{ci.val}{The values of the parameter at the confidence intervals.}
#' \item{ci}{The confidence interval value used.}
#'
#' @author David Orme
#' @seealso \code{\link{pgls}}
#' @keywords util stats
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(
#'     shorebird.tree, shorebird.data, Species,
#'     vcv = TRUE, vcv.dim = 3
#' )
#' mod <- pgls(log(Egg.Mass) ~ log(M.Mass), shorebird, lambda = "ML")
#' mod.l <- pgls.profile(mod, "lambda")
#' plot(mod.l)
#' pgls.confint(mod, "lambda")
#'
pgls.profile <- function(pgls, which = c("lambda", "kappa", "delta"),
                         N = 50, param.CI = NULL) {
    ## takes a pgls model and profiles one of the branch length transformations

    # get the x sequence for the parameter
    which <- match.arg(which)
    bnds <- pgls$bounds[[which]]
    x <- seq(from = bnds[1], to = bnds[2], length = N)

    # get a matrix of parameter values
    pars <- matrix(pgls$param, nrow = N, ncol = 3, byrow = TRUE)
    colnames(pars) <- names(pgls$param)
    pars[, which] <- x

    ## now get the sequence of likelihoods for the parameter in question
    logLik <- sapply(seq_along(x), function(X) {
        pgls.likelihood(
            optimPar = NULL, fixedPar = pars[X, ], y = pgls$y, x = pgls$x,
            V = pgls$data$vcv, optim.output = TRUE
        )
    })

    RET <- list(
        x = x, logLik = logLik, which = which, pars = pgls$param,
        dname = pgls$dname, formula = pgls$formula
    )
    class(RET) <- "pgls.profile"

    # test for existing parameter ci otherwise create if asked
    if (!is.null(pgls$param.CI[which])) {
        RET$ci <- pgls$param.CI[[which]]
    } else if (!is.null(param.CI)) {
        RET$ci <- pgls.confint(pgls, which, param.CI)
    }

    return(RET)
}

#' @describeIn pgls.profile Plot a pgls model profile object
#' @export
plot.pgls.profile <- function(x, ...) {
    xlab <- as.expression(x$which)
    xsub <- sprintf(
        "Data: %s; Model: %s\nkappa %0.2f; lambda %0.2f; delta %0.2f",
        x$dname, deparse(x$formula), x$pars["kappa"],
        x$pars["lambda"], x$pars["delta"]
    )

    with(x, plot(logLik ~ x, type = "l", xlab = xlab, ...))
    graphics::title(sub = xsub, cex.sub = 0.7, line = graphics::par("mgp")[1] + 1.5)


    if (!is.null(x$ci)) {
        graphics::abline(v = x$ci$opt, col = "red", ...)
        graphics::abline(v = x$ci$ci.val, lty = 2, col = "red", ...)
    }
}

#' @describeIn pgls.profile Confidence intervals on pgls branch transformations
#' @export
pgls.confint <- function(pgls, which = c("lambda", "kappa", "delta"),
                         param.CI = 0.95) {
    # Are we dealing with a same confidence interval
    # ha! first planned use of sequential 'or'  operator
    if (!is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) {
        stop("ci is not a number between 0 and 1.")
    }

    # find the parameter being checked
    which <- match.arg(which)

    # is the value in the object for this parameter an ML value?
    # - if not, then this needs to be estimated in order
    #   to get confidence intervals.
    # - currently, bail out but could refit to model to get this.
    if (pgls$mlVals[which] == FALSE) {
        stop("The pgls object contains a fixed, not ML, estimate of ", which)
    }
    ML <- pgls$model$log.lik

    # separate out the values held constant and varied
    fix <- pgls$param
    whichNum <- which(names(fix) == which)
    opt <- fix[whichNum]
    fix <- fix[-whichNum]

    # only one optimPar so get bounds and two intervals
    bounds <- pgls$bounds[[which]]
    belowML <- c(bounds[1], opt)
    aboveML <- c(opt, bounds[2])

    # the offset needed to get the root of the ML surface
    # at zero is  - (observed ML) + a chisq component

    MLdelta <- (stats::qchisq(param.CI, 1) / 2)
    offset <- (-ML) + MLdelta

    ## get the model components
    y <- pgls$y
    x <- pgls$x
    V <- pgls$data$vcv

    ## find the confidence intervals on the parameter
    ## - first need to find the logLik at the bounds
    ## - as long as the bound is outside the CI, can then use uniroot
    ##   to find the actual confidence interval.

    lowerBound.ll <- pgls.likelihood(
        structure(bounds[1], names = which), fix, y, x, V,
        optim.output = TRUE
    )
    upperBound.ll <- pgls.likelihood(
        structure(bounds[2], names = which), fix, y, x, V,
        optim.output = TRUE
    )

    lrt0 <- 2 * (ML - lowerBound.ll)
    lrt1 <- 2 * (ML - upperBound.ll)
    lowerBound.p <- 1 - stats::pchisq(lrt0, 1)
    upperBound.p <- 1 - stats::pchisq(lrt1, 1)

    # A problem with uniroot is that the identity of the variables gets stripped
    #   which is why pgls.likelihood now has an optim.names option used here.
    ll.fun <- function(opt) {
        pg <- pgls.likelihood(
            opt, fix, y, x, V,
            optim.output = TRUE, names.optim = which
        )
        ll <- pg + offset
        return(ll)
    }

    lowerCI <- if (lowerBound.ll < (ML - MLdelta)) {
        stats::uniroot(ll.fun, interval = belowML)$root
    } else {
        NA
    }
    upperCI <- if (upperBound.ll < (ML - MLdelta)) {
        stats::uniroot(ll.fun, interval = aboveML)$root
    } else {
        NA
    }

    return(list(
        opt = opt,
        bounds.val = bounds,
        bounds.p = c(lowerBound.p, upperBound.p),
        ci.val = c(lowerCI, upperCI),
        ci = param.CI
    ))
}
