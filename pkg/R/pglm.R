
## TODO - rename phylomat to V throughout to maintain consistenct with profile functions...

pglm <- function(formula, data, V = NULL, lambda = 1.0, kappa = 1.0,  delta= 1.0, 
	             param.CI = NULL, control = list(fnscale=-1), 
                 bounds = list(lambda=c(0,1), kappa=c(0,3), delta=c(0,3))) {

	# bounds go singular: bounds = list(delta = c(1e-04, 3), lambda = c(1e-05,  0.99999), kappa = c(0, 3))
	
	## think about this - allow old data + V use?
	if(! inherits(data, 'comparative.data')) stop("data is not a 'comparative' data object.")
	
	## if the comparative data doesn't contain a VCV,
	## then add one.
	if(is.null(data$vcv)){
		V <- vcv.array(data$phy)
	} else {
		V <- data$vcv
	}
	
	## What does this do?
	Dfun <- function(Cmat) {
		iCmat <- solve(Cmat,  tol = .Machine$double.eps)
		svdCmat <- La.svd(iCmat)
		D <- svdCmat$u %*% diag(sqrt( svdCmat$d )) %*% t(svdCmat$v)
		return( t(D) )
	}

	# Estimates the GLS parameters for given data
	get.coeffs <- function(Y, V, X) {
		iV <- solve(V, tol = .Machine$double.eps)
		xVix <- crossprod(X, iV %*% X)
		xViy <- crossprod(X, iV %*% Y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy 	#This is  a bad thing to do!!!!
		return(mu)
	}

	# Estimates the variance of a given trait (accounting for phylogeny)
	est.var <- function(y, V, x, mu ) {
		iV <- solve(V, tol = .Machine$double.eps)
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

	# Full ML estimation for given x and V
	# modified to also act as an engine for optim
	# - the function is passed named vectors containing kappa, lambda and delta
	#   which might be available for optimization (optimPar) or user defined (fixedPar)
	log.likelihood <- function(optimPar, fixedPar, y, x, V, optim.output=TRUE, names.optim=NULL) {
	    
	    # merge the values of KLD from the two parameter vectors
		# if names.optim is provided then add it (uniroot in the ci.search strips it out)
		if(! is.null(names.optim)) names(optimPar) <- names.optim
        allPar <- c(optimPar, fixedPar)
	    
        V <- blenTransform(V, allPar)

		mu <- get.coeffs(y, V, x)
		s2 <- est.var(y, V, x, mu)
		n <- length(x[,1])
		logDetV <- determinant(V, logarithm = TRUE)$modulus[1]
		ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - 1)/2.0

		ypred <- x%*%mu	

		# if being used for optimization, only return the log likelihood
		if(optim.output) return(ll)  else return( list(ll = ll, mu = mu, s2 = s2) )
	}

	null.var <- function(y, V) {
		X <- matrix(1, nrow = length(y))
		mu <- get.coeffs(y, V, X)
		return(est.var(y, V, X, mu))
	}

	## generalised function to search for the  
	## confidence intervals on the parameters
	## optimPar holds the ML estimate for the parameter to search
	## fixedPar holds the others
	
	ci.search <- function(optimPar, fixedPar, y, x, V, ML, ci=0.95, bounds=bounds){
		
		# only one optimPar so get bounds and two intervals
		nm <- names(optimPar)
		bounds  <- bounds[[nm]]
		belowML <- c(bounds[1], optimPar)
		aboveML <- c(optimPar, bounds[2])
		
		# the offset needed to get the root of the ML surface
		# at zero is  - (observed ML) + a chisq component
 		offset <- (- ML) + (qchisq(ci, 1)/2)
		
	    ll.fun <- function(optimPar){
	        pg <- log.likelihood(optimPar, fixedPar, y, x, V, optim.output=TRUE, names.optim=nm)
	        ll <- pg + offset
	        return(ll)
	    }
		
		## need to find out in advance whether the bounds are within the CI, or 
		## uniroot doesn't work
		## a problem with uniroot is that the identity of the variables gets stripped
		## this is why log.likelihood has an optim.names option
		lowerCI <- uniroot(ll.fun , interval=belowML)
		upperCI <- uniroot(ll.fun , interval=aboveML)

		return(cbind(lowerCI, upperCI))

	}

	#prune.dat <- prune(data, V) 
	#V <- prune.dat$V
	#data <- prune.dat$dat
	nm <- names(data)
	
	n <- length(data[,1])
	
	# Get the design matrix
	m <- model.frame(formula, data$data)
	y <- m[,1]
	x <- model.matrix(formula, m)
	k <- length(x[1,])
	
	namey <- names(m)[1]

	# sort out the branch length transformations
	# sensible values?
	parVals <- list(kappa=kappa, lambda=lambda, delta=delta)
	
	# check and sort bounds
	if(! setequal(names(bounds), c('kappa', 'lambda', 'delta'))){
		stop("Bounds does not contain elements labelled 'kappa','lambda' and 'delta'")
	}
	bounds <- bounds[c('kappa','lambda','delta')]
	
	for(i in seq_along(parVals)){
		p <- parVals[[i]]
		nm <- names(parVals)[i]
		lb <- bounds[[nm]][1]
		ub <- bounds[[nm]][2]
		
		if(length(p) > 1) stop(n, " not of length one.")
		if(is.character(p) & p != "ML") stop(nm, " is character and not 'ML'.")
		if(is.numeric(p) & ( p < lb | p > ub))
			stop(sprintf("%s is out of bounds [%0.1f, %0.1f]", nm, lb, ub))
	}
	
	if(kappa != 1 && length(dim(V)) != 3) stop("3D vcv.array needed for kappa transformation.")

	mlVals <- sapply(parVals,  "==", "ML")

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
    	
    	optim.param.vals <- optim(optimPar, fn = log.likelihood, # function and start vals
    	    method="L-BFGS-B", control=control, upper=upper.b, lower=lower.b, # optim control
    	    V = V, y=y, x=x, fixedPar = fixedPar, optim.output=TRUE) # arguments to function
	    
    	if(optim.param.vals$convergence != "0"){
    		stop("Problem with optim:", optim.param.vals$convergence, 
    		    optim.param.vals$message)}
	    
    	fixedPar <- c(optim.param.vals$par, fixedPar)
    	fixedPar <- fixedPar[c("kappa","lambda","delta")]
    } else {
        fixedPar <- as.numeric(parVals)
		names(fixedPar) <- c("kappa", "lambda","delta")
    }
    
	# run the likelihood function again with the fixed parameter values
	V <- blenTransform(V, fixedPar)
	ll <- log.likelihood(optimPar=NULL, fixedPar=fixedPar, y, x, V, optim=FALSE)

	log.lik <- ll$ll
	
	## if requested, get the confidence intervals on the optimized parameters
	## if any are actually optimised
	if(! is.null(param.CI) && any(mlVals)){
		
		# ha! first planned use of sequential 'or'  operator
		if(! is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) 
			stop('param.CI is not a number between 0 and 1.')
		
		browser()
		## Loop over optimized parameters
		param.CI.list <- list(kappa=NULL, lambda=NULL, delta=NULL)
		for(prm in which(mlVals)){
			opt <- fixedPar[prm]
			fix <- fixedPar[-prm]
			param.CI.list[[prm]] <- ci.search(opt,  fix, y, x, V, log.lik, ci=param.CI, bounds=bounds)
		}

	}
	


	aic <- -2 * log.lik + 2 * k
	aicc <- -2 * log.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
	
	coeffs <- ll$mu
	coeffs <- data.frame(t(coeffs))
	names(coeffs) <- colnames(x)
	varNames = names(m)

	
	pred <- x %*% ll$mu 
	
	res <- y - pred
	D <- Dfun(V)
	pres <- D %*% res
	
	fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
	
	logDetV <- determinant(V, logarithm = TRUE)$modulus[1]
	
	logLikY <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log( (n - k) * ll$s2 / n) - logDetV / 2.0  - n / 2.0

	RMS <- ll$s2
	RSSQ <- ll$s2 * (n - k)
	NMS <- RMS
	NSSQ <- RSSQ
	
	if(k > 0) {
		NMS <- null.var(y, V)
		NSSQ <- NMS * (n - 1)
		}

	# Bits for parameter errors	
	errMat <- t(x)%*% solve(V) %*% x  
	errMat <- solve(errMat) * RMS[1] 
	sterr <- diag(errMat)
	sterr <- sqrt(sterr)
	
	
	ret <- list(model = fm, formula = formula, logLikY = logLikY, RMS = RMS, NMS = NMS,
	            NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, aicc = aicc, n = n, k = k,
	            sterr = sterr, vcv = errMat, fitted = pred, residuals = res, phyres = pres,
	            x = x, data = data,  varNames = varNames, y = y, V = V, param = fixedPar, mlVals=mlVals,
	            L0 = NULL, L1 = NULL, LamOptimised = FALSE, namey = namey)

	class(ret) <- "pglm"
	return(ret)
	
}

