
phylo.d <- function(data, phy, names.col, binvar, permut=1000, rnd.bias=NULL) {

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

	# get the variable out and convert to ranks
	ds <- data$data[ ,bininds]
	if(any(is.na(ds))) stop("'", binvar, "' contains missing values.")
	dsRank <- rank(ds)
	
	# check for a number
    if (!is.numeric(permut)) (stop("'", permut, "' is not numeric.")) 
	
	# check tree branch lengths
	el    <- data$phy$edge.length
	elTip <- data$phy$edge[,2] <= length(data$phy$tip.label)
	
	if(any(el[elTip] == 0)) 
		stop('Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate')
	if(any(el[! elTip] == 0)) 
		stop('Phylogeny contains zero length internal branches. Use di2multi.')
	
	## This is rewritten away from the original version with internal functions
	##  - structure was slowing and the functions aren't externalised ever
	
	## Random Association model random data
	##  - with weighted shuffling if weights are given
	ds.ran <- replicate(permut, sample(dsRank))
	
	## Brownian Threshold model random data
	
		## there was a call to lambdaTree(phy,1) - why???	
		## - get the variance covariance for the tree
		if(is.null(data$vcv)){
			vcv <- VCV.array(data$phy)
		} else {
			vcv <- data$vcv
		}
		
		# Simulate traits up the tree
		ds.phy <- rmvnorm(permut, sigma=unclass(vcv)) # class of 'VCV.array' throws the method dispatch
		ds.phy <- as.data.frame(t(ds.phy))
		
		# turn those into rank values
		ds.phy <- apply(ds.phy, 2, rank)

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
		ds.ran <- cbind(Obs=dsRank, ds.ran)
		ds.phy <- cbind(Obs=dsRank, ds.phy)
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
		# double check the observed, but only to six decimal places or you can get floating point errors
		if(round(ransocc[1], digits=6) != round(physocc[1], digits=6)) stop('Problem with character change calculation in phylo.d')
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
				binvar = binvar,  data=data, nPermut = permut, rnd.bias=rnd.bias)
	
	class(dvals) <- 'phylo.d'
	return(dvals)
	
}