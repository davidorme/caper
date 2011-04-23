## ## testing code for the phylo.d method of contrCalc using
## ## the example tree from Fritz and Purvis 2010
## 
## FnPtree <- read.tree(text='((((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1):1,(((I:1,J:1):1,(K:1,L:1):1):1,((M:1,N:1):1,(O:1,P:1):1):1):1):1;')
## FnPtree$node.label <- 17:31
## FnPdat <- data.frame(extra.clumpy  = c(1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0),
## 					 brown.clumpy  = c(0,0,1,1,1,0,1,0,0,0,1,1,0,0,1,1),
## 					 random        = c(1,0,1,1,0,0,1,0,0,0,1,1,1,0,0,1),
## 					 overdispersed = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0))
## 
## FnPdat <- as.matrix(FnPdat)
## rownames(FnPdat) <- FnPtree$tip.label
## 
## diffs <- contrCalc(FnPdat, phy=FnPtree, ref.var='extra.clumpy', picMethod='phylo.d', crunch.brlen=0)
## 
## colSums(diffs$contrMat)
## ## and hey presto
## ## extra.clumpy  brown.clumpy        random overdispersed 
## ##          1.0           5.0           6.5           8.0

phylo.d <- function(data, phy, names.col, binvar, permut=1000) {

    # - test to see if there is a comparative data object and if not then
    #   retrofit the remaining arguments into a comparative data object.
	if(! missing(data)){
		if(! inherits(data, 'comparative.data')){
			if(missing(names.col)) stop('names column is missing')
			names.col <- deparse(substitute(names.col))
			data <- caicStyleArgs(data=data, phy=phy, names.col=names.col)
		}
	}
	
	# look for binary variable
	binvar <- deparse(substitute(binvar))
    bininds <- match(binvar, names(data$data))
    if (is.na(bininds)) (stop("'", binvar, "' is not a variable in data."))

	# get the variable out and do a general test for binarity
	ds <- data$data[ ,bininds]
	if(length(unique(ds)) != 2) stop("'", binvar, "' contains more than two states.")
	if(any(is.na(ds))) stop("'", binvar, "' contains missing values.")
	
	# proportion - applies to any two unique values
	propStates <- unclass(table(ds))
	propState1 <- propStates[1]/sum(propStates)
	
	# check for a number
    if (!is.numeric(permut)) (stop("'", permut, "' is not numeric.")) 
	
	## This is rewritten away from the original version with internal functions
	##  - structure was slowing and the functions aren't externalised ever
	
	## Random Association model random data
		ds.ran <- replicate(permut, sample(ds))
	
	## Brownian Threshold model random data
	
		## there was a call to lambdaTree(phy,1) - why???	
		## - get the variance covariance for the tree
		if(is.null(data$vcv)){
			vcv <- vcv.array(data$phy)
		} else {
			vcv <- data$vcv
		}
		
		# Simulate traits up the tree
		ds.phy <- rmvnorm(permut, sigma=unclass(vcv)) # class of 'vcv.array' throws the method dispatch
		ds.phy <- as.data.frame(t(ds.phy))
	
		## - find the threshold in each variable. 
		## - quantile interpolates between values
		ds.phy.thresh <- apply(ds.phy, 2, quantile, propState1)
		
		## sweep out the thresholds
		ds.phy <- sweep(ds.phy, 2, ds.phy.thresh, '<')
		ds.phy <- as.numeric(ds.phy) ## bah! kills dims so reinstate
		dim(ds.phy) <- dim(ds.ran)

	## Get change along edges
		
		## ## It is very slow to use crunch for big formulae because there
		## ## is a massive overhead (~ 95% of crunch run time) in using the model
		## ## formula apparatus for such large formulae. Although the code
		## ## below works it is a huge performance hit compared to just 
		## ## running through contrCalc directly. Advantage of comparative.data!
		## ds.ran <- cbind(Obs=ds, ds.ran)
		## ds.ran <- as.data.frame(ds.ran)
		## ## get default formulae
		## ds.ran.formula <- formula(ds.ran)
		## ## would be too paranoid to use the comparative data function rather than hacking the object!
		## ds.ran.CD <- data
		## ds.ran.CD$data <- ds.ran
		## ds.phy.caic <- crunch(ds.ran.formula, ds.ran.CD)
		
		## insert observed and set dimnames for contrCalc
		ds.ran <- cbind(Obs=ds, ds.ran)
		ds.phy <- cbind(Obs=ds, ds.phy)
		dimnames(ds.ran) <- dimnames(ds.phy) <- list(data$phy$tip.label, c('Obs', paste('V',1:permut, sep='')))
		
		## being careful with the edge order - pre-reorder the phylogeny
		## because the method won't reorder an already matching order.
		## Plus we need the pruningwise order later.
		phy <- reorder(data$phy, 'pruningwise')
		
		## now run that through the contrast engine 
		## - in fact, the change calculation requires a tree traversal to compare 
		##   change along the edges from the nodal values of the daughters to the parent
		##   and this traversal is what contrCalc does. So create a new contrCalc method.
		ds.ran.cc <- contrCalc(vals=ds.ran, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
		ds.phy.cc <- contrCalc(vals=ds.phy, phy=phy, ref.var='V1', picMethod='phylo.d', crunch.brlen=0)
		
	## get sums of change and distributions
	
		ransocc <- colSums(ds.ran.cc$contrMat)
		physocc <- colSums(ds.phy.cc$contrMat)
		# double check the observed
		if(ransocc[1] != physocc[1]) stop('Problem with character change calculation in phylo.d')
		obssocc <- ransocc[1]
		ransocc <- ransocc[-1]
		physocc <- physocc[-1]
		
		soccratio <- (obssocc - mean(physocc)) / (mean(ransocc) - mean(physocc))
		soccpval1 <- sum(ransocc < obssocc) / permut
		soccpval0 <- sum(physocc > obssocc) / permut
		
	
	dvals <- list(DEstimate=soccratio, Pval1=soccpval1, Pval0=soccpval0,
		        Parameters=list(Observed=obssocc, 
		        MeanRandom=mean(ransocc), MeanBrownian=mean(physocc)), 
		        Permutations=list(random=ransocc, brownian=physocc), 
		        NodalVals=list(observed = ds.ran.cc$nodVals[, 1,drop=FALSE], 
			                   random   = ds.ran.cc$nodVals[,-1,drop=FALSE], 
			                   brownian = ds.phy.cc$nodVals[,-1,drop=FALSE]),
				binvar = binvar,  data=data, nPermut = permut)
	
	class(dvals) <- 'phylo.d'
	return(dvals)
	
}

print.phylo.d <- function(x, ...){
    summary(x)
}

summary.phylo.d <- function(object, ...){
    cat('\nCalculation of D statistic for the phylogenetic structure of a binary variable\n')
    cat('\n  Data : ', object$data$data.name)
    cat('\n  Binary variable : ', object$binvar)
    cat('\n  Phylogeny : ', object$data$phy.name)
    cat('\n  Number of permutations : ', object$nPermut)
    
    cat("\n\nEstimated D : ", object$DEstimate)
    cat("\nProbability of E(D) resulting from no (random) phylogenetic structure : ", object$Pval1)
    cat("\nProbability of E(D) resulting from Brownian phylogenetic structure    : ", object$Pval0)
    cat("\n\n")
}