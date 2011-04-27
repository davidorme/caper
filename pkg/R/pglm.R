## RESTRUCTURE AND EXPANSION/MERGING OF PGLM CODE

pglm <- function(formula, data, V = NULL, lambda = 1.0, kappa = 1.0,  delta= 1.0, 
	             param.CI = NULL, control = list(fnscale=-1), 
                 bounds = list(lambda=c(1e-6,1), kappa=c(1e-6,3), delta=c(1e-6,3))) {

	## bounds go singular: bounds = list(delta = c(1e-04, 3), lambda = c(1e-04,  0.99999), kappa = c(1e-04, 3))
	
	## pglm replaces lik.lambda - exactly the same as a null model
	
	## all the internal functions that were here are now farmed out to externally accessible functions
	## - except because I don't know what it does!
	
	## What does this do?
	Dfun <- function(Cmat) {
		iCmat <- solve(Cmat,  tol = .Machine$double.eps)
		svdCmat <- La.svd(iCmat)
		D <- svdCmat$u %*% diag(sqrt( svdCmat$d )) %*% t(svdCmat$v)
		return( t(D) )
	}

	## END OF INTERNAL FUNCTIONS, START OF FUNCTION CODE

	## think about this - allow old data + V use?
	if(! inherits(data, 'comparative.data')) stop("data is not a 'comparative' data object.")
	
	dname <- deparse(substitute(data))
	
	## if the comparative data doesn't contain a VCV,
	## then add one.
	if(is.null(data$vcv)){
		V <- vcv.array(data$phy)
	} else {
		V <- data$vcv
	}

	## sort out the data
	nm <- names(data$data)
	n <- nrow(data$data)
	
	# Get the design matrix, number of parameters and response variable
	m <- model.frame(formula, data$data)
	y <- m[,1]
	x <- model.matrix(formula, m)
	k <- ncol(x)
	namey <- names(m)[1]

	# if a ci is specified, check (early) that it is sensible for use at the end!
	# ha! first planned use of sequential 'or'  operator
	if(! is.null(param.CI)){
		if(! is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) 
			stop('param.CI is not a number between 0 and 1.')
	}
	
	# check and sort elements of bounds
	if(! setequal(names(bounds), c('kappa', 'lambda', 'delta'))){
		stop("Bounds does not contain elements labelled 'kappa','lambda' and 'delta'")
	}
	bounds <- bounds[c('kappa','lambda','delta')]
	
	## check the branch length transformations to be applied
	## - gather into a named list: names are used throughout to 
	##   get the right values in the right place.
	parVals <- list(kappa=kappa, lambda=lambda, delta=delta)
	
	## - test the bounds and parameter values are sensible
	for(i in seq_along(parVals)){

		## is the parameter a single number or 'ML'
		p <- parVals[[i]]
		nm <- names(parVals)[i]

		if(length(p) > 1) stop(nm, " not of length one.")
		if(is.character(p) & p != "ML") stop(nm, " is character and not 'ML'.")

		## are the bounds of length 2, numeric and positive or zero
		bnds <- bounds[[nm]]
		if(length(bnds) > 2) stop("Bounds specified for ",nm, " not of length one.")
		if(! is.numeric(bnds)) stop("Non-numeric bounds specified for ",nm, ".")
		if(any(bnds < 0)) stop("Negative values in bounds specified for ",nm, ".")
		lb <- bnds[1]
		ub <- bnds[2]
		if(lb > ub) stop("Lower bound greater than upper bound for ",nm, ".")
		
		## are specified transforms (not 'ML') in range (also filters out negative transforms) 
		if(is.numeric(p) & ( p < lb | p > ub))
			stop(sprintf("%s value (%0.2f) is out of specified bounds [%0.2f, %0.2f]", nm, p, lb, ub))
	}
	
	if(kappa != 1 && length(dim(V)) != 3) stop("3D vcv.array needed for kappa transformation.")

	## which are being optimised
	mlVals <- sapply(parVals,  "==", "ML")

	## if any are being optimised then run pglm.logLik as a simple optimising function,
	## returning the logLik for a particular set of transformations
	##  - start the search for ML estimates from the midpoint of the specified bounds
	
	if(any(mlVals)){
	    
    	# isolate parameters to be optimized and set to a sensible start.
    	parVals[mlVals] <- lapply(bounds, mean)[mlVals]
		# collapse list down to a vector
    	parVals <- as.numeric(parVals)
    	names(parVals) <- c("kappa", "lambda","delta")
		
		# split them up
    	optimPar <- parVals[mlVals]
    	fixedPar <- parVals[!mlVals]
    	
    	# define the optimization bounds
    	lower.b <- sapply(bounds,  "[", 1)[mlVals]
    	upper.b <- sapply(bounds,  "[", 2)[mlVals]
    	
		## TODO - could isolate single optimisations here to use optimise() rather than optim()
		## likelihood function swapped out for externally visible one
    	optim.param.vals <- optim(optimPar, fn = pglm.loglik, # function and start vals
    	    method="L-BFGS-B", control=control, upper=upper.b, lower=lower.b, # optim control
    	    V = V, y=y, x=x, fixedPar = fixedPar, optim.output=TRUE) # arguments to function
	    
    	if(optim.param.vals$convergence != "0"){
    		stop("Problem with optim:", optim.param.vals$convergence, 
    		    optim.param.vals$message)}
	    
    	fixedPar <- c(optim.param.vals$par, fixedPar)
    	fixedPar <- fixedPar[c("kappa","lambda","delta")]
    } else {
		## reduce the list of parameters to a vector
        fixedPar <- as.numeric(parVals)
		names(fixedPar) <- c("kappa", "lambda","delta")
    }
    
	## run the likelihood function again with the fixed parameter values
	## ll <- log.likelihood(optimPar=NULL, fixedPar=fixedPar, y, x, V, optim=FALSE)
	ll <- pglm.loglik(optimPar=NULL, fixedPar=fixedPar, y, x, V, optim=FALSE)
	
	## store the log likelihood of the optimized solution for use in ci.searchs
	log.lik <- ll$ll
	
	## get the transformed vcv matrix for the fitted model for use
	## in calculating the remaining outputs.
	Vt <- blenTransform(V, fixedPar)
	
	## start collating outputs:
	
	## AIC
	aic <- -2 * log.lik + 2 * k
	aicc <- -2 * log.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
	
	## coefficients
	coeffs <- ll$mu
	names(coeffs) <- colnames(x)
	varNames <- names(m)

	## predicted values
	pred <- x %*% ll$mu 
	
	##residuals
	res <- y - pred
	D <- Dfun(Vt)
	pres <- D %*% res # TODO - what is this exactly
	
	## fitted model
	fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
	
	## log likelihood of the data given the transformed vcv matrix
	logDetV <- determinant(Vt, logarithm = TRUE)$modulus[1]
	logLikY <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log( (n - k) * ll$s2 / n) - logDetV / 2.0  - n / 2.0
	
	## various variances
	RMS <- ll$s2
	RSSQ <- ll$s2 * (n - k)
	
	## null model
	xdummy <- matrix(rep(1, length(y)))
	nullMod <- pglm.loglik(optimPar=NULL, fixedPar=fixedPar, y, xdummy, V, optim.output=FALSE)
	NMS <- nullMod$s2
	NSSQ <- nullMod$s2 * (n -1) 
	
	# Bits for parameter errors	
	errMat <- t(x)%*% solve(Vt) %*% x  
	errMat <- solve(errMat) * RMS[1] 
	sterr <- diag(errMat)
	sterr <- sqrt(sterr)
	
	
	RET <- list(model = fm, formula = formula, logLikY = logLikY, RMS = RMS, NMS = NMS,
	            NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, aicc = aicc, n = n, k = k,
	            sterr = sterr, vcv = errMat, fitted = pred, residuals = res, phyres = pres,
	            x = x, data = data,  varNames = varNames, y = y, param = fixedPar, mlVals=mlVals,
	            namey = namey, bounds=bounds, Vt=Vt, dname=dname)
	
	class(RET) <- "pglm"
	
	## if requested, get the confidence intervals on the optimized parameters
	## if any are actually optimised
	if(! is.null(param.CI) && any(mlVals)){
		
		## Loop over optimized parameters
		param.CI.list <- list(kappa=NULL, lambda=NULL, delta=NULL)
		mlNames <- names(mlVals)[which(mlVals)]
		
		for(param in mlNames){
			param.CI.list[[param]] <- pglm.confint(RET, param, param.CI)
		}
		
		RET$param.CI <- param.CI.list
	}
	
	return(RET)
	
}

pglm.profile <- function(pglm, which=c('lambda','kappa','delta'), N=50, param.CI=NULL){
	
	## takes a pglm model and profiles one of the branch length transformations
	
	# get the x sequence for the parameter
	which <- match.arg(which)
	bnds <- pglm$bounds[[which]]
	x <- seq(from=bnds[1], to=bnds[2], length=N)
	
	# get a matrix of parameter values
	pars <- matrix(pglm$param, nrow=N, ncol=3, byrow=TRUE)
	colnames(pars) <- names(pglm$param)
	pars[,which] <- x
	
	## now get the sequence of likelihoods for the parameter in question
	logLik <- sapply(seq_along(x), function(X){ pglm.loglik(optimPar=NULL, fixedPar=pars[X,], y=pglm$y, x=pglm$x, V=pglm$data$vcv, optim.output=TRUE)})
	
	RET <- list(x=x,logLik=logLik, which=which, pars=pglm$param, dname=pglm$dname, formula=pglm$formula)
	class(RET) <- 'pglm.profile'
	
	# test for existing parameter ci otherwise create if asked
	if(! is.null(pglm$param.CI[which])){
		RET$ci <- pglm$param.CI[[which]]
	} else if(! is.null(param.CI)){
		RET$ci <- pglm.confint(pglm, which, param.CI)
	} 

	return(RET)
}

plot.pglm.profile <- function(x, ...){
	
	xlab <- as.expression(x$which)
	xsub <- sprintf('Data: %s; Model: %s\nkappa %0.2f; lambda %0.2f; delta %0.2f', 
	                x$dname, deparse(x$formula), x$pars['kappa'], x$pars['lambda'], x$pars['delta'])
	
	with(x, plot(logLik ~ x, type='l', xlab=xlab, ...))
	title(sub=xsub, cex.sub=0.7, line=par('mgp')[1]+1.5)
	
	
	if(! is.null(x$ci)){
		abline(v=x$ci$opt, col='red', ...)
		abline(v=x$ci$ci.val, lty=2, col='red', ...)
	}

}

pglm.confint <- function(pglm, which=c('lambda','kappa','delta'), param.CI=0.95){
	
	# Are we dealing with a same confidence interval
	# ha! first planned use of sequential 'or'  operator
	if(! is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) 
		stop('ci is not a number between 0 and 1.')
	
	# find the parameter being checked
	which <- match.arg(which)
	
	# is the value in the object for this parameter an ML value?
	# - if not, then this needs to be estimated in order 
	#   to get confidence intervals.
	# - currently, bail out but could refit to model to get this.
	if(pglm$mlVals[which] == FALSE) stop('The pglm object contains a fixed, not ML, estimate of ', which)
	ML <- pglm$model$log.lik
	
	# separate out the values held constant and varied
	fix <- pglm$param
	whichNum <- which(names(fix) == which)
	opt <- fix[whichNum]
	fix <- fix[-whichNum]
	
	# only one optimPar so get bounds and two intervals
	bounds  <- pglm$bounds[[which]]
	belowML <- c(bounds[1], opt)
	aboveML <- c(opt, bounds[2])
	
	# the offset needed to get the root of the ML surface
	# at zero is  - (observed ML) + a chisq component
	
	MLdelta <- (qchisq(param.CI, 1)/2)
	offset <- (- ML) + MLdelta

	## get the model components
	y <- pglm$y
	x <- pglm$x
	V <- pglm$data$vcv
	
	## find the confidence intervals on the parameter
	## - first need to find the logLik at the bounds
	## - as long as the bound is outside the CI, can then use uniroot 
	##   to find the actual confidence interval.

	lowerBound.ll <- pglm.loglik(structure(bounds[1], names=which), fix, y, x, V, optim.output=TRUE)
	upperBound.ll <- pglm.loglik(structure(bounds[2], names=which), fix, y, x, V, optim.output=TRUE)
	
	lrt0 <- 2 * (ML - lowerBound.ll)
	lrt1 <- 2 * (ML - upperBound.ll)
	lowerBound.p <- 1 - pchisq(lrt0, 1)
	upperBound.p <- 1 - pchisq(lrt1, 1)
	
	## - a problem with uniroot is that the identity of the variables gets stripped
	##   which is why pglm.logLik now has an optim.names option used here.
	ll.fun <- function(opt){
        pg <- pglm.loglik(opt, fix, y, x, V, optim.output=TRUE, names.optim=which)
        ll <- pg + offset
        return(ll)
    }
	
	lowerCI <- if(lowerBound.ll < (ML -MLdelta)) uniroot(ll.fun , interval=belowML)$root else NA
	upperCI <- if(upperBound.ll < (ML -MLdelta)) uniroot(ll.fun , interval=aboveML)$root else NA

	return(list(opt=opt, bounds.val=bounds, bounds.p=c(lowerBound.p, upperBound.p), ci.val=c(lowerCI, upperCI), ci=param.CI))

}

pglm.loglik <- function(optimPar, fixedPar, y, x, V, optim.output=TRUE, names.optim=NULL) {
    
	# Full ML estimation for given x and V: modified to also act as an engine for optim
	# - this is why the branch length  parameters are passed as two chunks, so that
	#   the first acts as the targets for optimisation.
	# - the function is passed named vectors containing kappa, lambda and delta
	#   which might be available for optimization (optimPar) or user defined (fixedPar)

    # merge the values of KLD from the two parameter vectors
	# if names.optim is provided then add it (uniroot in the ci.search strips it out)
	
	# Estimates the GLS parameters for given data
	get.coeffs <- function(Y, iV, X) {
		xVix <- crossprod(X, iV %*% X)
		xViy <- crossprod(X, iV %*% Y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy 	#This is  a bad thing to do!!!!
		return(mu)
	}

	# Estimates the variance of a given trait (accounting for phylogeny)
	est.var <- function(y, iV, x, mu ) {
		e <- y - x %*% mu
		s2 <- crossprod(e, iV %*% e)
		n <- length(y) 
		k <- length(x[1,])
		return( s2 / (n- k) )
	}

	## applies the three branch length transformations to a VCV matrix
    blenTransform <- function(V, par){
		
        # apply transformations
        if(par["kappa"] == 0) V <- (V > 0) else V <-  V ^ par["kappa"] # kappa catching NA^0=1
        V <- apply(V, c(1,2), sum, na.rm=TRUE) # collapse 3D array
        V <- ifelse(upper.tri(V)+lower.tri(V), V * par["lambda"], V) # lambda
        if(par["delta"] == 0) V <- (V > 0) else V <-  V ^ par["delta"] # delta catching NA^0=1

		attr(V, 'blenTransform') <- par
    	return(V)
	}
	
	if(! is.null(names.optim)) names(optimPar) <- names.optim
    allPar <- c(optimPar, fixedPar)
    
	# get the transformed VCV matrix and its inverse
    V <- blenTransform(V, allPar)
	iV <- solve(V, tol = .Machine$double.eps)
	
	mu <- get.coeffs(y, iV, x)
	s2 <- est.var(y, iV, x, mu)
	n <- length(x[,1])
	logDetV <- determinant(V, logarithm = TRUE)$modulus[1]
	ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - 1)/2.0

	# if being used for optimization, only return the log likelihood
	if(optim.output) return(ll)  else return( list(ll = ll, mu = mu, s2 = s2) )
}

plot.pglm <- function(x, ...) {
	
	# layout(matrix(c(1,2,3,4), 2, 2, byrow = FALSE))
	res <- residuals(x, phylo = TRUE)
	res <- res / sqrt(var(res))[1]
	truehist(res, xlab = "Residual value (corrected for phylogeny)")
	qqnorm(res)
	abline(0, 1)
	plot(fitted(x), res, xlab = "Fitted value", ylab = "Residual value (corrected for phylogeny)"  )
	plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value")
}

summary.pglm <- function(object,...) {
		
	## coefficient matrix
	
	cf <- object$model$coef
	se <- object$sterr
	t  <- cf/se
	p  <- 2 * ( 1 - pt( abs(t), object$n - object$k) )

	coef <- cbind(cf,se,t,p)
	colnames(coef) <- c('Estimate','Std. Error','t value','Pr(>|t|)')
	
	
		testLambda <- function(pobj) {
			
			lrt0 <- 2 * (pobj$logLikY - pobj$L0)
			lrt1 <- 2 * (pobj$logLikY - pobj$L1)
			
			p0 <- 1 - pchisq(lrt0, 1)
			p1 <- 1 - pchisq(lrt1, 1)
			
			cat("     Test of Lambda = 0: chisq = ", lrt0, " P = ", p0, "\n")
			cat("     Test of Lambda = 1: chisq = ", lrt1, " P = ", p1, "\n")
			}
			
	
	cat("\n\n--------------------------------------------------------\n")
	cat("Summary of Generalised Least Squares, correcting for \n")
	cat("Phylogeny:\n\n")
	cat("Number of parameters = ", object$k,"\n")
	cat("Number of data points = ", object$n,"\n\n")
	cat("Branch length transformations:\n")
    cat(sprintf("%-6s : %0.4f [%s]", names(object$param), object$param, ifelse(object$mlVals, "ML", "Fixed")), sep="\n")
	if(object$LamOptimised == TRUE) { testLambda(object)}
	cat("Maximised log-likelihood = ", object$logLikY,"\n\n")
	cat("Model AIC = ", object$aic, "\n")
	cat("Model AICc = ", object$aicc, "\n\n")

	cat("Null Mean Square = ", object$NSSQ,"\n")
	cat("Residual Mean Square = ", object$RSSQ, "\n\n")
	cat("Raw R^2 = ", (object$NSSQ - object$RSSQ) / object$NSSQ, "\n")
	cat("Adjusted R^2 = ", (object$NMS - object$RMS) / object$NMS, "\n")
	
	Fstat <- ((object$NSSQ - object$RSSQ) / object$RMS) / (object$k - 1)
	
	cat("F statistic = ",  ((object$NSSQ - object$RSSQ) / object$RMS) / (object$k - 1), " ")
	cat("P model = ", pf(Fstat, object$k - 1, object$n - object$k,  ncp=0, lower.tail = FALSE, log.p = FALSE), "\n\n")
	cat("Summary of coefficients:\n\n")
	coeffs <- coef(object)
	errs <- object$sterr
	cat("Term\t\tEstimate\t\tStd Err\t\tT-value\tP\n")
	storet<-c()
	for(i in 1:length(coeffs) ) {
		est <- coeffs[1,i]
		nm <- names(coeffs)[i]
		se <- errs[i]
		Tstat <- est / se
		storet<-c(storet,Tstat)	
		Pval <- 2 * ( 1 - pt( abs(Tstat), object$n - object$k) )
		cat(nm,"\t")
		cat(est, "\t", se, "\t", Tstat, "\t", Pval, "\n")
		}
	
	cat("\n\n--------------------------------------------------------\n")
}

