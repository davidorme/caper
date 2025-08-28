#' Phylogenetic generalized linear models
#'
#' Fits a linear model, taking into account phylogenetic non-independence
#' between data points. The strength and type of the phylogenetic signal in the
#' data matrix can also be accounted for by adjusting branch length
#' transformations (lambda, delta and kappa). These transformations can also be
#' optimised to find the maximum likelihood transformation given the data and
#' the model.
#'
#' This function fits a linear model controlling for the non-independence
#' between cases resulting from phylogenetic structure in the data. The
#' stucture of the phylogenetic signal can be controlled by altering the
#' parameters lambda, delta and kappa (see the 'caper' vignette for details).
#' The implementation of the method is currently as described in Freckleton et
#' al (2002).
#'
#' The branch length transformations can be optimised between bounds using
#' maximum likelihood by setting the value for a transformation to 'ML'. The
#' default bounds are: lambda = c(1e-6,1), kappa = c(1e-6,3) and
#' delta=c(1e-6,3). These defaults may be overridden by passing a named list
#' with new elements to the bounds argument - only the bounds to be changed
#' need to be provided (e.g. bounds=list(lambda=c(0,3))).
#'
#' The 'pgls.likelihood' and 'pgls.blenTransform' methods are not primarily
#' intended to be called by users. The 'pgls.likelihood' function provides a
#' general method to calculate the likelihood of a model, given the covariance
#' matrix, response, design matrix and branch length parameters.
#'
#' @aliases pgls pgls.likelihood pgls.blenTransform
#' @param formula A model formula
#' @param data A 'comparative.data' object containing the covariance matrix and
#' data to be used in the model.
#' @param lambda A value for the lambda transformation.
#' @param kappa A value for the kappa transformation.
#' @param delta A value for the delta transformation.
#' @param param.CI A p value used to calculate confidence intervals.
#' @param control A list of control parameters for the optim function.
#' @param bounds A list of bounds to use for branch length transformations (see
#' Details).
#' @param optimPar A named vector of branch length parameters to be optimised
#' to find the maximum likelihood value.
#' @param fixedPar A named vector of fixed values for branch length parameters.
#' @param y A column matrix of the model response.
#' @param x The design matrix of the model.
#' @param V A phylogenetic covariance matrix.
#' @param optim.output A logical value. If true then 'pgls.likelihood' returns
#' only the likelihood value for use in the 'optim' function.
#' @param names.optim The name of a single parameter being optimised. This is
#' only required for estimating parameter confidence intervals, where the
#' function 'uniroot' strips names from vectors.
#' @return The 'pgls' function returns an object of class \code{pgls}
#' containing the following:
#'
#' "na.action" "param.CI" \item{call}{The original call to the 'pgls' function}
#' \item{model }{A summary of the fitted model containing:} \item{formula}{The
#' model formula supplied.} \item{data}{The comparative data object provided.}
#' \item{dname}{The name of the comparative data object.} \item{logLikY}{The
#' log likelihood of the response variable given the model.} \item{RMS}{The
#' residual mean square variance in the model.} \item{RSSQ}{The residual sum of
#' squares from the model.} \item{NMS}{The null mean square variance for the
#' model.} \item{NSSQ}{The null sum of squares for the response.}
#' \item{aic}{The AIC score of the model} \item{aicc}{The AICc score of the
#' model, correcting for the number of cases and parameters estimated}
#' \item{n}{The number of rows of data used in fitting the model} \item{k}{The
#' number of parameter estimates} \item{sterr}{The standard errors of the
#' parameter estimates} \item{Vt}{The phylogenetic covariance matrix used in
#' the model, with branch length transformations applied.} \item{fitted}{The
#' predicted values} \item{residuals}{The non-phylogenetic residuals}
#' \item{phyres}{The phylogenetic residuals} \item{x}{The design matrix of the
#' model } \item{varNames}{The variables include in the model.} \item{y}{The
#' response of the model.} \item{namey}{The name of the response variable.}
#' \item{param}{A named numeric vector of length three giving the branch length
#' transformations used in the model.} \item{mlVals}{A named logical vector of
#' length three indicating which branch length values in 'param' are maximum
#' likelihood estimates.} \item{bounds}{The bounds on branch length parameter
#' estimates used in the model.} \item{param.CI}{A named list of length three
#' giving confidence intervals and the p values at the parameter bounds for
#' optimised branch length transformations. Fixed parameters will have a NULL
#' entry in this list.} \item{na.action}{A named vector identifying any rows of
#' missing data excluded from the model.}
#' @section Warning: The model is fitted using a data frame reduced to complete
#' row cases to eliminate missing values. In order to ensure that the models
#' fitted using different subsets of the data are comparable, the whole data
#' frame \code{data} is reduced to complete cases. In the future, a scope
#' argument may be provided to control this but at present the data frame
#' should be reduced to only those variables used in the maximal model in order
#' to avoid prevent redundant variables causing rows to be dropped
#' unnecessarily.
#' @author Rob Freckleton; David Orme
#' @seealso \code{\link{pgls.profile}}, \code{\link{anova.pgls}},
#' \code{\link{summary.pgls}}
#' @references R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic
#' analysis and comparative data: A test and review of evidence. American
#' Naturalist, 160:712-726, 2002.
#' @keywords models regression
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv = TRUE, vcv.dim = 3)
#' mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird, lambda = "ML")
#' mod2 <- pgls(log(Egg.Mass) ~ log(M.Mass), data = shorebird, lambda = "ML", delta = "ML")
#' @export
pgls <- function(formula, data, lambda = 1.0, kappa = 1.0, delta = 1.0,
                 param.CI = 0.95, control = list(fnscale = -1),
                 bounds = NULL) {
    # Local Dfun function
    Dfun <- function(Cmat) {
        iCmat <- solve(Cmat, tol = .Machine$double.eps)
        svdCmat <- La.svd(iCmat)
        D <- svdCmat$u %*% diag(sqrt(svdCmat$d)) %*% t(svdCmat$v)
        return(t(D))
    }

    ## think about this - allow old data + V use?
    if (!inherits(data, "comparative.data")) {
        stop("data is not a 'comparative' data object.")
    }
    dname <- deparse(substitute(data))
    call <- match.call()

    ## check for missing data in the formula and replace the data
    # object with a complete version
    miss <- model.frame(formula, data$data, na.action = na.pass)
    miss.na <- apply(miss, 1, function(X) (any(is.na(X))))
    if (any(miss.na)) {
        miss.names <- data$phy$tip.label[miss.na]
        data <- data[-which(miss.na), ]
    }

    # Get the design matrix, number of parameters and response variable
    m <- model.frame(formula, data$data)
    y <- m[, 1]
    x <- model.matrix(formula, m)
    k <- ncol(x)
    namey <- names(m)[1]

    # test for variables with no variance in model matrices
    # (thx to Sarah Dryhurst)
    # - will cause singularity - lm() filters for this and aliases variables
    #   but here we'll just fail for the time being
    xVar <- apply(x, 2, var)[-1] # drop intercept
    badCols <- xVar < .Machine$double.eps
    if (any(badCols)) {
        stop(
            "Model matrix contains columns with zero variance: ",
            paste(names(xVar)[badCols], collapse = ", ")
        )
    }

    ## if the comparative data doesn't contain a VCV,
    ## then get one and put it in the data object too. Bit wasteful
    if (is.null(data$vcv)) {
        V <- if (kappa == 1) {
            VCV.array(data$phy)
        } else {
            VCV.array(data$phy, dim = 3)
        }
        data$vcv <- V
    } else {
        V <- data$vcv
    }

    ## sort out the data
    nm <- names(data$data)
    n <- nrow(data$data)

    # if a ci is specified, check (early) that it is sensible for
    # use at the end!
    if (!is.null(param.CI)) {
        if (!is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) {
            stop("param.CI is not a number between 0 and 1.")
        }
    }

    # check and insert elements from user bounds
    usrBounds <- bounds
    bounds <- list(kappa = c(1e-6, 3), lambda = c(1e-6, 1), delta = c(1e-6, 3))
    if (!is.null(usrBounds)) {
        if (!is.list(usrBounds)) {
            stop(
                "Bounds must be a list of named bounds for any ",
                "or all of kappa, lambda and delta"
            )
        }

        usrNames <- names(usrBounds)
        badNames <- setdiff(usrNames, c("kappa", "lambda", "delta"))
        if (length(badNames) > 0) {
            stop(
                "The list of bounds contains names other than ",
                "kappa, lambda and delta"
            )
        }

        for (nm in usrNames) {
            bounds[nm] <- usrBounds[nm]
        }
    }

    ## check the branch length transformations to be applied
    ## - gather into a named list: names are used throughout to
    ##   get the right values in the right place.
    parVals <- list(kappa = kappa, lambda = lambda, delta = delta)

    ## - test the bounds and parameter values are sensible
    for (i in seq_along(parVals)) {
        ## is the parameter a single number or 'ML'
        p <- parVals[[i]]
        nm <- names(parVals)[i]

        if (length(p) > 1) stop(nm, " not of length one.")
        if (is.character(p) && p != "ML") {
            stop(nm, " is character and not 'ML'.")
        }

        ## are the bounds of length 2, numeric and positive or zero
        bnds <- bounds[[nm]]
        if (length(bnds) > 2) {
            stop("Bounds specified for ", nm, " not of length one.")
        }
        if (!is.numeric(bnds)) {
            stop("Non-numeric bounds specified for ", nm, ".")
        }
        if (any(bnds < 0)) {
            stop("Negative values in bounds specified for ", nm, ".")
        }
        lb <- bnds[1]
        ub <- bnds[2]
        if (lb > ub) {
            stop("Lower bound greater than upper bound for ", nm, ".")
        }

        ## are specified transforms (not 'ML') in range (also filters
        #  out negative transforms)
        if (is.numeric(p) && (p < lb || p > ub)) {
            stop(sprintf(
                "%s value (%0.2f) is out of specified bounds [%0.2f, %0.2f]",
                nm, p, lb, ub
            ))
        }
    }

    if (kappa != 1 && length(dim(V)) != 3) {
        stop("3D VCV.array needed for kappa transformation.")
    }

    ## which are being optimised
    mlVals <- sapply(parVals, "==", "ML")

    ## if any are being optimised then run pgls.likelihood as a
    # simple optimising function,
    ## returning the logLik for a particular set of transformations
    ##  - start the search for ML estimates from the midpoint of
    # the specified bounds

    if (any(mlVals)) {
        # isolate parameters to be optimized and set to a sensible start.
        parVals[mlVals] <- lapply(bounds, mean)[mlVals]
        # collapse list down to a vector
        parVals <- as.numeric(parVals)
        names(parVals) <- c("kappa", "lambda", "delta")

        # split them up
        optimPar <- parVals[mlVals]
        fixedPar <- parVals[!mlVals]

        # define the optimization bounds
        lower.b <- sapply(bounds, "[", 1)[mlVals]
        upper.b <- sapply(bounds, "[", 2)[mlVals]

        # TODO - could isolate single optimisations here to use optimise()
        # rather than optim() likelihood function swapped out for externally
        # visible one
        optim.param.vals <- optim(optimPar,
            fn = pgls.likelihood, # function and start vals
            method = "L-BFGS-B", control = control,
            upper = upper.b, lower = lower.b,
            V = V, y = y, x = x, fixedPar = fixedPar, optim.output = TRUE
        ) # arguments to function

        if (optim.param.vals$convergence != "0") {
            stop(
                "Problem with optim:", optim.param.vals$convergence,
                optim.param.vals$message
            )
        }

        fixedPar <- c(optim.param.vals$par, fixedPar)
        fixedPar <- fixedPar[c("kappa", "lambda", "delta")]
    } else {
        ## reduce the list of parameters to a vector
        fixedPar <- as.numeric(parVals)
        names(fixedPar) <- c("kappa", "lambda", "delta")
    }

    ## run the likelihood function again with the fixed parameter values
    ll <- pgls.likelihood(
        optimPar = NULL, fixedPar = fixedPar, y, x, V, optim.output = FALSE
    )

    ## store the log likelihood of the optimized solution for use in ci.searchs
    log.lik <- ll$ll

    ## get the transformed vcv matrix for the fitted model for use
    ## in calculating the remaining outputs.
    Vt <- pgls.blenTransform(V, fixedPar)

    ## start collating outputs:

    ## AIC
    aic <- -2 * log.lik + 2 * k
    aicc <- -2 * log.lik + 2 * k + ((2 * k * (k + 1)) / (n - k - 1))

    ## coefficients
    coeffs <- ll$mu
    names(coeffs) <- colnames(x)
    varNames <- names(m)

    ## predicted values
    pred <- x %*% ll$mu

    ## residuals
    res <- y - pred
    D <- Dfun(Vt)
    pres <- D %*% res # TODO - what is this exactly

    ## fitted model
    fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)

    ## various variances
    RMS <- ll$s2
    RSSQ <- ll$s2 * (n - k)

    ## null model
    xdummy <- matrix(rep(1, length(y)))
    nullMod <- pgls.likelihood(
        optimPar = NULL, fixedPar = fixedPar, y, xdummy, V, optim.output = FALSE
    )
    NMS <- nullMod$s2
    NSSQ <- nullMod$s2 * (n - 1)

    # Bits for parameter errors
    errMat <- t(x) %*% solve(Vt) %*% x
    errMat <- solve(errMat) * RMS[1]
    sterr <- diag(errMat)
    sterr <- sqrt(sterr)


    RET <- list(
        model = fm, formula = formula, call = call, RMS = RMS, NMS = NMS,
        NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, aicc = aicc, n = n, k = k,
        sterr = sterr, fitted = pred, residuals = res, phyres = pres,
        x = x, data = data, varNames = varNames, y = y, param = fixedPar,
        mlVals = mlVals,
        namey = namey, bounds = bounds, Vt = Vt, dname = dname
    )

    class(RET) <- "pgls"

    ## missing data
    if (any(miss.na)) {
        RET$na.action <- structure(
            which(miss.na),
            class = "omit", .Names = miss.names
        )
    }
    ## if requested, get the confidence intervals on the optimized parameters
    ## if any are actually optimised
    if (!is.null(param.CI) && any(mlVals)) {
        ## Loop over optimized parameters
        param.CI.list <- list(kappa = NULL, lambda = NULL, delta = NULL)
        mlNames <- names(mlVals)[which(mlVals)]

        for (param in mlNames) {
            param.CI.list[[param]] <- pgls.confint(RET, param, param.CI)
        }

        RET$param.CI <- param.CI.list
    }

    return(RET)
}



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
#' \item{x}{Parameter values at which the likelihood has been calculated.}
#' \item{logLik}{The likelihood value at each value.} \item{which}{The
#' parameter being profiled.} \item{pars}{The value of the other fixed
#' parameters.} \item{dname}{The name of the 'comparative.data' object used to
#' fit the model.} \item{formula}{The formula of the model being profiled}
#'
#' If the model contains an ML estimate of the parameter being profiled, then
#' the 'pgls.profile' object will also contain the output of 'pgls.confint':
#'
#' \item{opt}{The maximum likelihood value of the parameter.}
#' \item{bounds.val}{The values of the bounds on the parameter.}
#' \item{bounds.p}{The p value of the likelihood at the bounds, given the ML
#' value.} \item{ci.val}{The values of the parameter at the confidence
#' intervals.} \item{ci}{The confidence interval value used.}
#' @author David Orme
#' @seealso \code{\link{pgls}}
#' @keywords util stats
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv = TRUE, vcv.dim = 3)
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

plot.pgls.profile <- function(x, ...) {
    xlab <- as.expression(x$which)
    xsub <- sprintf(
        "Data: %s; Model: %s\nkappa %0.2f; lambda %0.2f; delta %0.2f",
        x$dname, deparse(x$formula), x$pars["kappa"],
        x$pars["lambda"], x$pars["delta"]
    )

    with(x, plot(logLik ~ x, type = "l", xlab = xlab, ...))
    title(sub = xsub, cex.sub = 0.7, line = par("mgp")[1] + 1.5)


    if (!is.null(x$ci)) {
        abline(v = x$ci$opt, col = "red", ...)
        abline(v = x$ci$ci.val, lty = 2, col = "red", ...)
    }
}

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

    MLdelta <- (qchisq(param.CI, 1) / 2)
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
    lowerBound.p <- 1 - pchisq(lrt0, 1)
    upperBound.p <- 1 - pchisq(lrt1, 1)

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
        uniroot(ll.fun, interval = belowML)$root
    } else {
        NA
    }
    upperCI <- if (upperBound.ll < (ML - MLdelta)) {
        uniroot(ll.fun, interval = aboveML)$root
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

pgls.likelihood <- function(optimPar, fixedPar, y, x, V,
                            optim.output = TRUE, names.optim = NULL) {
    # Full ML estimation for given x and V: modified to also act as an engine
    # for optim
    # - this is why the branch length  parameters are passed as two chunks, so
    #   that the first acts as the targets for optimisation.
    # - the function is passed named vectors containing kappa, lambda and delta
    #   which might be available for optimization (optimPar) or user defined
    #   (fixedPar)

    # merge the values of KLD from the two parameter vectors if names.optim is
    # provided then add it (uniroot in the ci.search strips it out)

    # Estimates the GLS parameters for given data
    get.coeffs <- function(Y, iV, X) {
        xVix <- crossprod(X, iV %*% X)
        xViy <- crossprod(X, iV %*% Y)
        # This is  a bad thing to do!!!!
        mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy
        return(mu)
    }

    # Estimates the variance of a given trait (accounting for phylogeny)
    est.var <- function(y, iV, x, mu) {
        e <- y - x %*% mu
        s2 <- crossprod(e, iV %*% e)
        n <- length(y)
        k <- length(x[1, ])
        return(s2 / (n - k))
    }

    if (!is.null(names.optim)) names(optimPar) <- names.optim
    allPar <- c(optimPar, fixedPar)

    # get the transformed VCV matrix and its inverse
    V <- pgls.blenTransform(V, allPar)
    iV <- solve(V, tol = .Machine$double.eps)

    mu <- get.coeffs(y, iV, x)
    s2 <- est.var(y, iV, x, mu)
    n <- nrow(x)
    k <- ncol(x)
    logDetV <- determinant(V, logarithm = TRUE)$modulus[1]

    ## Likelihood calculation
    ll <- (-n / 2.0 * log(2 * pi) -
        n / 2.0 * log((n - k) * s2 / n) -
        logDetV / 2.0 - n / 2.0)

    # if being used for optimization, only return the log likelihood
    if (optim.output) {
        return(ll)
    } else {
        return(list(ll = ll, mu = mu, s2 = s2))
    }
}

pgls.blenTransform <- function(V, fixedPar) {
    ## applies the three branch length transformations to a VCV matrix

    # apply transformations
    if (!is.null(fixedPar["kappa"]) && fixedPar["kappa"] != 1) {
        if (length(dim(V)) < 3) {
            stop("Kappa transformation requires a 3 dimensional VCV array.")
        }
    }

    if (fixedPar["kappa"] == 0) {
        # kappa catching NA^0=1
        V <- (V > 0)
    } else {
        V <- V^fixedPar["kappa"]
    }

    V <- apply(V, c(1, 2), sum, na.rm = TRUE) # collapse 3D array
    V <- ifelse(upper.tri(V) + lower.tri(V), V * fixedPar["lambda"], V) # lambda
    if (fixedPar["delta"] == 0) {
        # delta catching NA^0=1
        V <- (V > 0)
    } else {
        V <- V^fixedPar["delta"]
    }

    attr(V, "blenTransform") <- fixedPar
    return(V)
}



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
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv = TRUE, vcv.dim = 3)
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

coef.pgls <- function(object, ...) {
    cf <- object$model$coef
    nm <- rownames(cf)
    cf <- structure(as.vector(cf), names = nm)
    return(cf)
}

residuals.pgls <- function(object, phylo = FALSE, ...) {
    ret <- NULL
    if (phylo == FALSE) {
        ret <- object$res
    } else {
        ret <- object$phyres
    }
    return(ret)
}

fitted.pgls <- function(object, ...) {
    ret <- object$fitted
    return(ret)
}


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


logLik.pgls <- function(object, REML = FALSE, ...) {
    val <- object$model$log.lik

    attr(val, "nall") <- object$n
    attr(val, "nobs") <- object$n
    attr(val, "df") <- object$k
    class(val) <- "logLik"
    val
}

nobs.pgls <- function(object, ...) length(resid(object))




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
#' shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv = TRUE, vcv.dim = 3)
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
