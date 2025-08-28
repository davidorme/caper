#' Tree simulation with traits.
#'
#' This function provides a very general environment in which to simulate
#' trees. The basic philosophy is that the user provides a series of
#' expressions that define speciation rates, extinction rates and trait
#' evolution. These expressions can make use of information about the internal
#' state of the tree, allowing for very flexible definitions of rules for tree
#' growth.
#'
#' The main idea behind this function (which is still in development) is to
#' provide a flexible framework for simulating tree growth and trait evolution.
#' The user provides expressions for the main arguments (\code{b}, \code{d} and
#' \code{halt}) which act as rules defining speciation and extinction and the
#' ending of the simulation. These can be simple constants, but can also make
#' use of the properties of the environment of the evolving tree. This includes
#' both lineage specific properties (as described in the \code{lineages}
#' section of the returned value) or properties of the clade as a whole (as
#' described in the \code{clade} section of the returned value). For example, a
#' extinction rate might increase with lineage age (\code{d=0.01*lin.age}) or a
#' speciation rate might decrease according to a density dependent process
#' (\code{b=1 - (nExtantTip/500)}). Halt expressions will typically use clade
#' properties (\code{halt=clade.age >= 5} or \code{nTips >= 50}) but could use
#' lineage properties, for example stopping when a trait value hits a certain
#' value (\code{halt=any(ct1 >= 10)}). It is not permitted to use '==' in a
#' halt function of clade.age because it will allow the simulation to run away
#' if the actual value steps over the test value.
#'
#' Discrete traits are defined using a matrices of rates for transitions
#' between states for each trait. At present, these are fixed for the duration
#' of a simulation epoch and cannot be set as expressions of tree variables.
#'
#' Continuous trait evolution currently employs a simple Brownian model, given
#' a starting value and variance per unit time. The traits can have defined
#' co-variance (the simulation uses mvrnorm at present) and can also have a
#' defined mean change, allowing for a directional walk in the trait values. At
#' present, it is not possible for the trait variance to vary according to the
#' internal state of the tree; continuous characters retain the same variance
#' and covariance for the whole of the simulation epoch.
#'
#' Whilst none of the \code{halt} rules are TRUE, then the function evaluates
#' the birth, death and discrete trait rates and converts these to waiting
#' times using random variates from a exponential distribution with the
#' calculated rates. These competing waiting times are compared both to each
#' other and the \code{grain} of the simulation, the shortest waiting time is
#' found and the relevant event is then triggered. The winning event is
#' identified in the character vector \code{winnerName} in order to allow
#' inheritance rules to differentiate events.
#'
#' @aliases growTree as.comparative.data.growTree
#' @param b A speciation rate. This can be a numeric constant, as in the
#' default, which specifies a single speciation rate for the simulation.
#' Alternatively, this can be an expression, or a list of expressions which
#' define speciation rate in terms of the properties of the tree. See details
#' for discussion of those properties.
#' @param d An extinction rate, described as above.
#' @param halt A rule use to halt the simulation. The default is the number of
#' tips in the simulation, specified as a single integer, but this can also be
#' an expression or list of expressions on the properties of the tree. The
#' simulation is halted when any of these expressions becomes true.
#' @param grain Where rates depend on time or trait values, it becomes
#' necessary to allow time to pass discretely in order to re-evaluate waiting
#' times under the changing values. This sets the amount of time that is
#' allowed to pass before re-evaluation. If rates do not depend on such
#' changing parameters, it is sensible to set this to infinity - this will
#' ensure that the flow of the simulation is not slowed by checking.
#' @param linObj This can be used to supply an existing simulation object,
#' which will then continue to grow under the provided rules. This allows the
#' user to simulate trees with different sets of rules operating in different
#' epochs. The function \code{linToApe} will convert such an object to a
#' 'phylo' object, retaining additional trait data as extra components of the
#' 'phylo' object list.
#' @param ct.start A numeric vector specifying the starting values for
#' continuous traits. If unnamed these will be sequentially named as 'ct1',
#' 'ct2' etc. The names of traits may be used in expressions governing tree
#' growth rules.
#' @param ct.change A numeric vector describing the mean change per unit time
#' in continuous trait values, used to simulate a directional bias in character
#' evolution. If ct.change is NULL, then this is assumed to be zero for each
#' species.
#' @param ct.var Either a vector of variances for each trait or a square matrix
#' describing the variances and covariances amongst the continuous traits. If
#' this is NULL, then uncorrelated traits with a variance of 1 are assumed.
#' @param dt.rates A list of matrices describing the rate of transition between
#' discrete character traits. Each matrix defines a trait and, as with
#' ct.start, the list names are used to identify the traits in the simulation
#' and default to 'dt1', 'dt2', etc. The dimnames of the matrix are used to
#' identify the states of the trait and default to 'st1', 'st2', etc. The
#' matrix need not be symmetrical: the rates are defined from the states in the
#' columns to the states in the rows, hence the diagonal should probably be
#' zero. Each trait is assumed to start the simulation in the first state in
#' the matrix.
#' @param inheritance A list of rules that are applied after a speciation and
#' can be used to modify trait values for the descendent lineages. The names of
#' the list specify which traits are to be modified and, for each trait
#' specified, should return a vector of length two which replaces the existing
#' values.
#' @param trace.events A logical flag, indicating whether or not report
#' speciation, extinction and discrete character evolution events.
#' @param trace.cladesize A positive integer giving an increment size for the
#' simulation to report clade size if required.
#' @param output.lineages A logical flag indicating whether to return the
#' internal lineages object.
#' @param neg.rates One of 'abort', 'warn' or 'quiet', defining the behaviour
#' when a rate calculation produces a negative number. With 'warn' and 'quiet',
#' negative rates are set to zero and the simulation continues.
#' @param inf.rates One of 'abort', 'warn' or 'quiet', as for\code{neg.rates}.
#' With 'warn' or 'quiet', infinite rates are left in place, resulting in
#' events happening instantly. This may, in some cases, be desirable!
#' @param stall.time If the all rates within the simulation are zero then only
#' this length of time is allowed to pass before the simulation exits with a
#' 'stalled' status. If \code{grain} is infinite, then the simulation stalls
#' immediately when all rates are zero.
#' @param extend.proportion This option allows the simulation to continue
#' running for a given proportion of the time to the next speciation. This
#' makes sense when growing a clade to a given number of extant taxa; with the
#' default setting of zero, the resulting tree ends at a bifurcation with zero
#' branch lengths and this option allows the tree to grow (and taxa to go
#' extinct and traits to evolve).
#' @param x A lineage table output from growTree
#' @param ... Further arguments to as.
#' @return Depending on the value of \code{output.phylo}, either an object of
#' class 'phylo' or an object of class 'growTree' with the following structure:
#' \item{lineages}{A data frame with a row for each lineage in the tree. Each
#' row identifies the \code{parent.id} and \code{id} of the row along with the
#' total age of the lineage (\code{lin.age}) and the times at which the lineage
#' was born (\code{birth.time}). If the species went extinct (or speciated)
#' then the \code{death.time} is recorded and \code{extinct} is set to TRUE.
#' Speciating lineages have \code{tip} set to FALSE. Each row also records the
#' \code{caic.code} of the lineage - this is used as a sorting code for
#' conversion to a 'phylo' object and is a kludge. If traits are defined in the
#' simulation then the values or states are recorded in this table. These are
#' the current values for extant tips and the values at extinction for extinct
#' tips and internal nodes.} \item{clade}{A list containing:\code{clade.age},
#' the total age of the simulation; \code{nLin}, the total number of lineages;
#' \code{nTip}, the total number of tips, differentated into \code{nExtantTip}
#' and \code{nExtinctTip}.} \item{rules}{A list reporting the birth (\code{b}),
#' death (\code{d}) and stopping (\code{halt}) rules and any inheritance
#' rules.} \item{ct.set}{If continuous characters were simulated, a list of the
#' \code{ct.start}, \code{ct.change} and \code{ct.var} details provided.}
#' \item{dt.rates}{If discrete characters were simulated, a list containing the
#' \code{dt.rates} details provided.}
#' @author David Orme, drawing heavily on discussions with Paul-Michael Agapow.
#' @keywords utilities
#' @examples
#'
#'
#' ## see the package vignette for a much fuller discussion of examples.
#'
#' # A basic 200 tip tree, output as a 'comparative.data' object
#' tree <- growTree(halt = 200, grain = Inf)
#' plot(tree$phy)
#'
#' # A basic tree of age 4 time units, output as a 'comparative.data' object
#' tree <- growTree(halt = expression(clade.age >= 4), grain = Inf)
#' plot(tree$phy)
#' @export
growTree <- function(b = 1, d = 0, halt = 20, grain = 0.1, linObj = NULL,
                     ct.start = NULL, ct.change = NULL, ct.var = NULL,
                     dt.rates = NULL, inheritance = NULL, trace.events = FALSE,
                     trace.cladesize = FALSE, output.lineages = FALSE,
                     neg.rates = "abort", inf.rates = "abort",
                     stall.time = 10, extend.proportion = 0) {
    waitTime <- function(rates) {
        # rexp can't handle rates of zero or Inf without a warning so handle
        # that in producing waiting times
        zero <- rates == 0
        inf <- rates == Inf
        rates[zero | inf] <- 1
        wait <- rexp(length(rates), rates)
        wait[zero] <- Inf
        wait[inf] <- 0

        # handle matrix of rates, such as the discrete trait rates uses...
        if (is.matrix(rates)) {
            wait <- array(wait, dim = dim(rates), dimnames = dimnames(rates))
        }

        return(wait)
    }

    ## The simulation has the following clade level properties:
    ## - clade.age, nLin, nTip, nExtantTip, nExtinctTip
    ## and the following lineage properties:
    ## - id, parent.id, lin.age, birth.time, death.time, extinct, tip

    ## CREATE A LINEAGE OBJECT OR USE THE PROVIDED ONE
    if (is.null(linObj)) {
        ## IF NO EXISTING LINEAGES ARE PROVIDED, THEN INTIALIZE A NEW ONE
        lineages <- data.frame(
            parent.id = 0, id = 1, lin.age = 0, birth.time = 0,
            death.time = as.numeric(NA), extinct = FALSE, tip = TRUE,
            caic.code = "", stringsAsFactors = FALSE
        )
        clade <- list(
            clade.age = 0, nLin = 1, nTip = 1,
            nExtantTip = 1, nExtinctTip = 0
        )
        linStartAge <- 0
        epochRules <- NULL

        ## CONTINUOUS TRAIT SETUP
        if (!is.null(ct.start)) {
            ctFlag <- TRUE

            ## check we have the right sort of info to pass to
            if (!is.numeric(ct.start)) {
                stop(
                    "'ct.start' must be a numeric vector providing the ",
                    "starting values for continuous traits."
                )
            }
            if (is.null(ct.change)) {
                ct.change <- rep(0, length(ct.start))
            }
            if (!is.numeric(ct.change)) {
                stop(
                    "'ct.change' must be a numeric vector providing the mean ",
                    "change per unit time for continuous traits."
                )
            }
            if (is.null(ct.var)) {
                ct.var <- diag(rep(1, length(ct.start)))
            }
            if (!is.numeric(ct.var)) {
                stop(
                    "'ct.var' must be a numeric vector  or matrix providing ",
                    "the (co-)variances for continuous traits."
                )
            }

            if (length(ct.start) != length(ct.change)) {
                stop("ct.start and ct.change must be the same length")
            }

            if (is.matrix(ct.var)) {
                if (dim(ct.var)[1] != dim(ct.var)[2]) {
                    stop("ct.var must be a square matrix")
                }
                if (dim(ct.var)[1] != length(ct.start)) {
                    stop(
                        "The dimensions of ct.var must match the ",
                        "length of ct.start"
                    )
                }
            } else {
                if (length(ct.start) != length(ct.var)) {
                    stop("ct.start and ct.var must be the same length")
                }
                ct.varMat <- array(0, dim = rep(length(ct.var), 2))
                diag(ct.varMat) <- ct.var
                ct.var <- ct.varMat
            }

            if (is.null(names(ct.start))) {
                names(ct.start) <- paste("ct", seq(along = ct.start), sep = "")
            }

            ## add the traits to the lineages
            lineages <- cbind(lineages, t(ct.start))
        } else {
            ctFlag <- FALSE
        }

        ## DISCRETE TRAIT SETUP
        if (!is.null(dt.rates)) {
            ## check that the transition matrix matches the number of levels
            if (!is.list(dt.rates) || !all(sapply(dt.rates, is.matrix))) {
                stop(
                    "dt.rates must be a list of matrices giving ",
                    "the rates of change between states"
                )
            }
            rateDims <- sapply(dt.rates, function(X) dim(X))
            if (!all(rateDims[1, ] - rateDims[2, ] == 0)) {
                stop(
                    "dt.rates must be a list of _square_ matrices ",
                    "giving the rates of change between states"
                )
            }
            if (any(unlist(dt.rates) < 0)) {
                stop("Negative values in dt.rates matrix")
            }

            dtFlag <- TRUE

            ## name the traits
            if (is.null(names(dt.rates))) {
                names(dt.rates) <- paste("dt", seq(along = dt.rates), sep = "")
            }
            ## name the states
            dt.rates <- lapply(dt.rates, function(X) {
                if (is.null(dimnames(X))) {
                    dimnmX <- paste("st", seq(to = dim(X)[1]), sep = "")
                    dimnames(X) <- list(dimnmX, dimnmX)
                }
                return(X)
            })

            ## add the traits to the lineages
            ## - get a data.frame of one row with first state in each matrix
            ##   as a factor with levels set by the trait states...
            dt.start <- lapply(
                dt.rates,
                function(X) {
                    factor(dimnames(X)[[1]][1],
                        levels = dimnames(X)[[1]]
                    )
                }
            )
            lineages <- cbind(lineages, dt.start)
        } else {
            dtFlag <- FALSE
        }
    } else {
        if (!inherits(linObj, "growTree")) {
            stop("Lineage object is not of class 'growTree'.")
        }
        clade <- linObj$clade
        lineages <- linObj$lineage
        epochRules <- linObj$epochRules
        linStartAge <- linObj$clade$clade.age

        # TODO - some checking to make sure that possible traits are catered for
        #        or assume that they are with the ct slot on the lineages object
        #      - current setup is to use the rules from the last epoch unless
        #        overridden
        #      - continuous - only looks for a change in the ct.change for
        #        overriding

        lastRules <- linObj$epochRules[[length(linObj$epochRules)]]

        if (is.null(ct.change)) {
            if (is.null(lastRules$ct.set)) {
                ## no continuous trait rules in previous epoch
                ctFlag <- FALSE
            } else {
                ## copy in previous rule set
                # think about this - ct.start makes no sense beyond the first
                # epoch but currently carries the trait names
                ct.start <- lastRules$ct.set$ct.start
                ct.change <- lastRules$ct.set$ct.change
                ct.var <- lastRules$ct.set$ct.var
                ctFlag <- TRUE
            }
        } else {
            ## TODO -  check that new rules make sense
            ctFlag <- TRUE
        }

        if (is.null(dt.rates)) {
            if (is.null(lastRules$dt.rates)) {
                ## no discrete trait rules in previous epoch
                dtFlag <- FALSE
            } else {
                ## copy in previous rule set
                dt.rates <- lastRules$dt.rates
                dtFlag <- TRUE
            }
        } else {
            ## TODO -  check that new rules make sense
            dtFlag <- TRUE
        }
    }

    # CHECKING RATE INPUTS
    # - can be a numeric constant, an expression, a list of expressions
    #   or a named object containing one of those things...
    # - expressions can only use the available properties of lineages and clades
    validVar <- c(names(lineages), names(clade))
    validExpr <- function(X) {
        if (any(is.na(match(all.vars(X), validVar)))) FALSE else TRUE
    }

    ## convert them all into a list of 'rules'
    ## SPECIATION RULE(S)
    switch(mode(b),
        "numeric" = if (length(b) > 1 | b < 0) {
            stop("Speciation rate 'b' is numeric but not a positive scalar")
        } else {
            b <- list(as.expression(b))
        },
        "expression" = if (!validExpr(b)) {
            stop("Speciation expression 'b' contains unknown variables.")
        } else {
            b <- list(b)
        },
        "list" = {
            if (!all(sapply(b, mode) == "expression")) {
                stop("The list 'b' must be a list of expressions")
            }
            if (!all(sapply(b, validExpr))) {
                stop(
                    "One or more expressions in the list 'b' contain ",
                    "unknown variables"
                )
            }
        },
        stop("Speciation rate 'b' not in a recognized format.")
    )

    ## assign names to the list
    if (is.null(names(b))) names(b) <- paste("b", seq(along = b), sep = "")

    ## EXTINCTION RULE(S)
    switch(mode(d),
        "numeric" = if (length(d) > 1 | d < 0) {
            stop("Extinction rate 'd' is numeric but not a positive scalar")
        } else {
            d <- list(as.expression(d))
        },
        "expression" = if (!validExpr(d)) {
            stop("Extinction expression 'd' contains unknown variables.")
        } else {
            d <- list(d)
        },
        "list" = {
            if (!all(sapply(d, mode) == "expression")) {
                stop("The list 'd' must be a list of expressions")
            }
            if (!all(sapply(d, validExpr))) {
                stop(
                    "One or more expressions in the list 'd' contain ",
                    "unknown variables"
                )
            }
        },
        stop("Extinction rate 'd' not in a recognized format.")
    )

    ## assign names to the list
    if (is.null(names(d))) names(d) <- paste("d", seq(along = d), sep = "")

    ## HALT RULE(S)
    ## default is to convert a single numeric value into a number of extant tups
    switch(mode(halt),
        "numeric" = if (length(halt) > 1 | halt <= 1) {
            stop(
                "If numeric, 'halt' must be a single number greater than ",
                "1 giving the number of extant tips to create."
            )
        } else {
            halt <- list(as.expression(
                substitute(nExtantTip >= XXX, list(XXX = halt))
            ))
        },
        "expression" = if (!validExpr(halt)) {
            stop("Stopping expression 'halt' contains unknown variables.")
        } else {
            halt <- list(halt)
        },
        "list" = {
            if (!all(sapply(halt, mode) == "expression")) {
                stop("The list 'halt' must be a list of expressions")
            }
            if (!all(sapply(halt, validExpr))) {
                stop(
                    "One or more expressions in the list 'halt' contain ",
                    "unknown variables"
                )
            }
        },
        stop("Stopping expression 'halt' not in a recognized format.")
    )

    ## test for clade.age ==
    haltCheck <- deparse(halt)
    if (length(grep("clade.age ?==", haltCheck)) > 0) {
        stop("'clade.age == ' used in halt expression. Use 'clade.age >='.")
    }

    ## assign names to the list
    if (is.null(names(halt))) {
        names(halt) <- paste("halt", seq(along = halt), sep = "")
    }

    ## check for behaviour on encountering negative or infinite rates
    neg.rates <- match.arg(neg.rates, c("abort", "warn", "quiet"))
    inf.rates <- match.arg(inf.rates, c("abort", "warn", "quiet"))


    ## now start simulation - keep the simulation running whilst all of the halt
    ## expressions are FALSE and while something is alive and , approximately,
    ## whilst anything is actually happening
    status <- "complete"
    lastRealEvent <- linStartAge

    ## A function to get a rate for each lineage from either b or d
    ## - is also used to extend the simulation after the last speciation
    ratesCheck <- function(rateExp, lineages, inf.rates, neg.rates) {
        ## - force multiplication by unity to extend constants across the clade
        unitConst <- rep(1, length(lineages$id))
        rates <- lapply(
            rateExp,
            function(X) eval(X, envir = c(lineages, clade)) * unitConst
        )

        ## negative rates?
        if (any(unlist(rates) < 0)) {
            switch(neg.rates,
                abort = stop("Negative rates produced"),
                warn = warning("Negative rates produced - setting to zero")
            )
        }

        rates <- lapply(rates, function(X) ifelse(X < 0, 0, X))

        ## infinite rates
        if (any(is.infinite(unlist(rates)))) {
            switch(inf.rates,
                abort = stop("Infinite rates produced"),
                warn = warning("Infinite rates produced")
            )
        }

        ## ensure extinct stay dead
        rates <- lapply(rates,
            function(X, ext) ifelse(ext, 0, X),
            ext = lineages$extinct
        )

        return(rates)
    }

    while (!any(haltStatus <- sapply(halt, eval, envir = c(lineages, clade)))) {
        ## evaluate the birth and death rate expressions
        bRates <- ratesCheck(b, lineages, inf.rates, neg.rates)
        dRates <- ratesCheck(d, lineages, inf.rates, neg.rates)

        allRates <- c(unlist(bRates), unlist(dRates))

        ## check to see if the stall criteria are met...
        if (all(allRates == 0)) {
            if ((clade$clade.age - lastRealEvent) > stall.time) {
                status <- "stalled"
                warning(
                    "Rates are all zero and stall.time is exceeded: ",
                    "exiting stalled simulation"
                )
                break # so end the simulation
            }
            if (grain == Inf) {
                warning(
                    "All rates are zero and grain is set to infinity ",
                    "giving no finite waiting times: exiting stalled simulation"
                )
                status <- "stalled"
                break
            }
        }

        ## make sure something is alive...
        if (all(lineages$extinct)) {
            status <- "extinct"
            warning("All lineages extinct: exiting simuation")
            break # everything is dead so end the simulation
        }


        ## turn rates into waiting times...
        bWait <- lapply(bRates, waitTime)
        dWait <- lapply(dRates, waitTime)

        # look at discrete trait changes TODO - rethink this along the lapply
        # lines - not sure it can be easily done...
        if (dtFlag) {
            # for each trait, take the relevant column in lineages and use it to
            # sample columns from the appropriate rate matrix
            #  - vector to hold fastest time to change
            dtWait <- numeric(length(dt.rates))
            names(dtWait) <- names(dt.rates)
            #  - matrix to hold which lineage is to change to what state for
            #    those fastest events
            dtWhich <- rbind(changeToState = dtWait, lineage = dtWait)

            for (dt in names(dt.rates)) {
                currDtRate <- dt.rates[[dt]][, lineages[, dt], drop = FALSE]
                colnames(currDtRate) <- lineages$id
                currDtRate[, lineages$extinct] <- 0
                currDtWait <- waitTime(currDtRate)
                dtWait[dt] <- min(currDtWait)
                dtWhich[, dt] <- which(
                    currDtWait == dtWait[dt],
                    arr.ind = TRUE
                )[1, ]
                # TODO - think about that [1,] - removes ties but these are
                # always likely to be between Inf so this is reasonable
            }
        } else {
            dtWait <- Inf
            dtWhich <- 0
        }

        ## ... look for the winning event...
        firstB_ID <- sapply(bWait, which.min)
        firstD_ID <- sapply(dWait, which.min)
        firstB_Time <- sapply(bWait, min)
        firstD_Time <- sapply(dWait, min)

        # order meaningful here - breaks ties in this order
        competWait <- c(firstB_Time, firstD_Time, dtWait, grain)
        # if grain wins the race then no row will be selected...
        competID <- c(firstB_ID, firstD_ID, dtWhich, 0)
        competType <- c(
            rep("Spec", length(bWait)),
            rep("Ext", length(dWait)),
            rep("Discrete", length(dtWait)), "Grain"
        )
        competName <- c(names(b), names(d), names(dtWait), "Grain")

        winner <- which.min(competWait)
        winnerWait <- competWait[winner]
        winnerID <- competID[winner]
        winnerType <- competType[winner]
        winnerName <- competName[winner]

        ## ... allow time to pass for the clade and extant lineages...
        clade$clade.age <- clade$clade.age + winnerWait
        lineages$lin.age[
            !lineages$extinct
        ] <- lineages$lin.age[!lineages$extinct] + winnerWait

        ## trait changes?
        if (ctFlag) {
            ## need some code in here
            delta <- as.matrix(
                MASS::mvrnorm(
                    n = dim(lineages)[1],
                    mu = ct.change,
                    Sigma = ct.var
                ) * winnerWait
            )
            delta[lineages$extinct, ] <- 0 # the dead don't evolve
            lineages[, names(ct.start)] <- lineages[, names(ct.start)] + delta
        }

        ## ... and maybe something happens...
        if (winnerType == "Spec") { # a birth!

            # copy parent into the daughters (incidentally inheriting any
            # traits...)
            daughterID <- clade$nLin + (1:2)
            parent <- lineages[winnerID, ]
            daughterInfo <- rbind(parent, parent)

            ## inheritance rules go in here
            if (!is.null(inheritance)) {
                for (x in seq(along = inheritance)) {
                    daughterInfo[
                        , names(inheritance)[x]
                    ] <- eval(inheritance[[x]], envir = daughterInfo)
                }
            }

            daughterInfo$id <- daughterID
            daughterInfo$parent.id <- parent$id
            daughterInfo$birth.time <- clade$clade.age
            daughterInfo$lin.age <- 0
            daughterInfo$caic.code <- paste(daughterInfo$caic.code,
                c("A", "B"),
                sep = ""
            )
            rownames(daughterInfo) <- daughterID
            lineages <- rbind(lineages, daughterInfo)

            ## kill the parent
            lineages$extinct[winnerID] <- TRUE
            lineages$tip[winnerID] <- FALSE
            lineages$death.time[winnerID] <- clade$clade.age

            ## update the clade
            clade$nLin <- clade$nLin + 2
            clade$nTip <- clade$nTip + 1
            clade$nExtantTip <- clade$nExtantTip + 1
            ## record that something happened here
            lastRealEvent <- clade$clade.age
            if (trace.events) {
                cat(clade$clade.age, ": Lineage ", winnerID, " speciated\n")
            }
        }

        if (winnerType == "Ext") { # an extinction
            ## kill the lineage
            lineages$extinct[winnerID] <- TRUE
            lineages$death.time[winnerID] <- clade$clade.age
            ## update the clade
            clade$nExtantTip <- clade$nExtantTip - 1
            clade$nExtinctTip <- clade$nExtinctTip + 1
            ## record that something happened here
            lastRealEvent <- clade$clade.age
            if (trace.events) {
                cat(clade$clade.age, ": Lineage ", winnerID, " went extinct\n")
            }
        }

        if (winnerType == "Discrete") { # a discrete trait changes

            lin <- dtWhich["lineage", winnerName]
            # as a number
            state <- dtWhich["changeToState", winnerName]
            # as the state name
            state <- dimnames(dt.rates[[winnerName]])[[1]][state]
            lineages[lin, winnerName] <- state
            lastRealEvent <- clade$clade.age
            if (trace.events) {
                cat(
                    clade$clade.age, ": Discrete trait ", winnerName,
                    " on lineage", lin, " changed to state ", state, "\n"
                )
            }
        }
        if (trace.cladesize) {
            if (clade$nTip %% trace.cladesize == 0) {
                cat(clade$clade.age, ": Clade contains ", clade$nTip, "tips\n")
            }
        }
    }

    if (trace.events) {
        cat(
            clade$clade.age, ":Simulation halted by halt rule '",
            names(halt)[haltStatus], "'\n"
        )
    }

    ## create a lineage structure
    RET <- list(lineages = lineages, clade = clade, status = status)

    ## get the settings for the epoch
    newEpochRules <- list(
        epochStart = linStartAge, b = b, d = d,
        halt = halt, inheritance = inheritance
    )
    if (ctFlag) {
        newEpochRules <- c(newEpochRules, list(ct.set = list(
            ct.start = ct.start, ct.change = ct.change, ct.var = ct.var
        )))
    }
    if (dtFlag) {
        newEpochRules <- c(newEpochRules, list(dt.rates = dt.rates))
    }
    RET$epochRules <- c(epochRules, list(newEpochRules))

    class(RET) <- "growTree"

    # if the user wants to allow clade growth to continue for a further
    # proportion (p) of the waiting time to the next clade. Generate a waiting
    # time (W) and then run again with speciation set to zero, allowing the tree
    # to grow to the current clade age + W*p, with the same extinction and
    # character change rules TODO? - build this into a helper function
    # extendTree() rather than doing it internally...
    if (extend.proportion > 0) {
        if (trace.events) {
            cat(
                clade$clade.age,
                ": Simulation extended beyond last speciation.\n"
            )
        }

        ## generate a waiting time
        bRates <- ratesCheck(b, lineages, inf.rates, neg.rates)
        bWait <- lapply(bRates, waitTime)
        timeToStop <- clade$clade.age + min(unlist(bWait)) * extend.proportion
        haltExpr <- as.expression(
            substitute(clade.age >= XXX, list(XXX = timeToStop))
        )

        # grow the tree a bit more with no births, but the other processess
        # running as usual but set the grain to a non infinite value to avoid
        # stalls on infinite waiting times
        # - TODO: this grain choice is arbitrary
        # - note that the trait evolution information is carried over via the
        #   epoch rules recording
        RET <- growTree(
            b = 0, d = d, halt = haltExpr,
            grain = if (!is.finite(grain)) {
                0.001
            } else {
                grain
            },
            linObj = RET,
            ct.start = NULL, ct.change = NULL, ct.var = NULL, dt.rates = NULL,
            inheritance = NULL, trace.events = trace.events,
            output.lineages = TRUE,
            neg.rates = neg.rates, inf.rates = inf.rates,
            stall.time = stall.time, extend.proportion = 0
        )
    }

    class(RET) <- "growTree"
    if (!output.lineages) RET <- as.comparative.data(RET)

    return(RET)
}


#' @describeIn growTree Convert a growTree object to a comparative data object
as.comparative.data.growTree <- function(x, ...) {
    lineages <- x$lineages
    clade <- x$clade

    # get a phylo object

    if (dim(lineages)[1] > 1) {
        # extract information (excluding root node)
        linNoR <- lineages[-1, ]

        # now need to coerce the numbering into ape style
        parentMap <- data.frame(
            linPar = unique(linNoR$parent.id),
            apePar = with(clade, (nTip + 1):nLin)
        )
        childMap <- data.frame(
            linPar = with(linNoR, id[tip]),
            apePar = with(clade, 1:nTip)
        )
        nodeMap <- rbind(parentMap, childMap)
        linNoR$pnode <- nodeMap$apePar[match(linNoR$parent.id, nodeMap$linPar)]
        linNoR$node <- nodeMap$apePar[match(linNoR$id, nodeMap$linPar)]

        edge <- as.matrix(linNoR[, c("pnode", "node")])
        dimnames(edge) <- NULL

        # lets be honest... the caic.code column is really only there as a cheap
        # route to a cladewise sorting of the edge matrix!
        ord <- order(linNoR$caic.code)
        edge <- edge[ord, ]
        edge.length <- linNoR$lin.age[ord]

        phy <- list(
            edge = edge, edge.length = edge.length, tip.label = 1:clade$nTip,
            Nnode = with(clade, nLin - nTip), root.edge = lineages$lin.age[1]
        )
        class(phy) <- "phylo"
    } else {
        # have an unspeciated root so put in slightly  as a single tip with a
        # root edge of zero
        phy <- list(
            edge = matrix(c(1, 2), ncol = 2), edge.length = lineages$lin.age[1],
            tip.label = 1, root.edge = 0, Nnode = 1
        )
        class(phy) <- "phylo"
        # extract data for use in comparative data objects.
        linNoR <- lineages
        names(linNoR)[2] <- "node"
    }

    # sort out comparative data
    lastRules <- x$epochRules[[length(x$epochRules)]]
    datCol <- c(
        "node", "lin.age", "birth.time", "death.time", "extinct", "tip",
        names(lastRules$ct.set$ct.start), names(lastRules$dt.rates)
    )
    dat <- linNoR[, datCol]

    # get tip data set without tips flag
    tipData <- dat[dat$tip, -6]
    nodeData <- dat[!dat$tip, -5:-6]
    com <- comparative.data(phy, tipData, "node", na.omit = FALSE)
    com$node.data <- nodeData
    com$epochRules <- x$epochRules
    attr(com, "growTree") <- TRUE

    return(com)
}
