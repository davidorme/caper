caic.label <- function(phy, charset=NULL, action="insert", style="CAIC"){
    
    # OLD2NEW STATUS: CONVERTED...

    if(! inherits(phy, "phylo")) 
         stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
        
    match.arg(action, c("insert", "replace", "append"))
    match.arg(style, c("RLE", "CAIC"))
    
    contrGp <- split(phy$edge[,2], f=phy$edge[,1]) # handily, split retains numeric order not alphabetic...
    caicLab <- character(max(phy$edge)) 
    names(caicLab) <- 1:max(phy$edge)


    if(is.null(charset)) charset <- LETTERS

    # loop the nodes
    for(nd in seq(along=contrGp)){

        parent <- names(contrGp)[nd]
        children <- contrGp[[nd]]
        if(length(children) > length(charset)) stop("Insufficient characters to code polytomies")
        caicLab[children] <- paste(caicLab[parent], charset[1:length(children)], sep="")

    }

    if(style=="RLE"){
        caicLab <- strsplit(caicLab, split="")
        caicLab <- sapply(caicLab, function(X) with(rle(X), paste(ifelse(lengths > 1, lengths, ""), values, sep="", collapse="")))
    }

    # put in the root label
    caicLab[caicLab == ""] <- "@Root"

    # OLD2NEW: intBool <- as.numeric(names(caicLab)) < 0 # internal nodes now from max(phy$edge)-phy$Nnode +1 to  max(phy$edge)
    intBool <- with(phy, 1:max(edge) > (max(edge) - Nnode))

    # insert option changed from an ordered match to edge[,2] in order to preserve the root label
    switch(action, 
        "replace" = { phy$tip.label <- caicLab[! intBool]
                      phy$node.label <- caicLab[intBool]},
        "append"  = { if(is.null(phy$node.label)) phy$node.label <- with(phy, (max(edge)-Nnode +1):max(edge))
                      phy$tip.label <- paste(phy$tip.label, caicLab[! intBool])
                      phy$node.label <- paste(phy$node.label, caicLab[intBool])},
        "insert"  =   phy$edge.caic.code <- caicLab) #[match(phy$edge[,2], names(caicLab))])

    return(phy)
}

caic.table <- function(caicObj, validNodes=TRUE, nodalValues=FALSE, ultrametric.tol=0.0001, CAIC.codes=FALSE){
    
		# simple code to create a data frame of the contrasts from the caic object
		rowID <- names(caicObj$contrast.data$contrVar)
		
        contr <- with(caicObj$contrast.data$contr, cbind(response, explanatory))

        if(nodalValues){
            nv <- with(caicObj$contrast.data$nodalVals, cbind(response, explanatory))
            colnames(nv) <- paste("nodal.", colnames(nv), sep="")
            tab <- as.data.frame(cbind(contr, nv), row.names=rowID)
        } else {
            tab <- as.data.frame(cbind(contr), row.names=rowID)
        }
        
        tab$contrVar <- caicObj$contrast.data$contrVar
        tab$validNodes <- caicObj$contrast.data$validNodes
        tab$nChild <- caicObj$contrast.data$nChild
        tab$nodeDepth <- caicObj$contrast.data$nodeDepth
        if(is.ultrametric(caicObj$data$phy, tol=ultrametric.tol)) {
			tab$nodeAge <- branching.times(caicObj$data$phy) 
		} else { 
			tab$nodeAge <- NA
		} 
		
		# not sure the match() calls here are necessary
        stRes <- rstudent(caicObj$mod)
        tab$studResid <- NA
        tab$studResid[match(names(stRes), rowID)] <- stRes

        if(validNodes) tab <- subset(tab, validNodes, select=-validNodes)
       
       if(CAIC.codes){
           Cphy <- caic.label(caicObj$data$phy)
           tab$CAIC.code <- Cphy$edge.caic.code[match(rowID, names(Cphy$edge.caic.code))]
       }
       
       return(tab)

}

summary.caic <- function(object, ...){

    summary(object$mod, ...)

}

print.caic <- function(x, ...){

    cat("Phylogenetic Independent Contrasts analysis using ",  attr(x, "contr.method"), ".\n\n", sep="")

    cat("Phylogeny: ", attr(x, "phyName"), " (",  length(x$data$phy$tip.label)  ," tips)\n", sep="")
    cat("Data: ",  attr(x, "dataName"), " (",  nrow(x$data$data)  ," rows)\n", sep="")
    cat("Number of valid contrasts: ", sum(x$contrast.data$validNodes), "\n", sep="")

    print(summary(x))

}

predict.caic <- function(object, ...){
    
    # need to force the model to get predictions using the contrast table rather than the original data table...
    # don't completely hijack the newdata argument...
    
    dots <- list(...)
    newdataProv <- pmatch(names(dots), "newdata")
    if(all(is.na(newdataProv))) nD <- caic.table(object) else nD <- dots[[newdataProv]]
    predict(object$mod, newdata=nD)
    
    
}

## AIC is completely agnostic about the model types fed into it.
AIC.caic <- function(object, ..., k=2){
	
	# borrowing heavily from AIC.default
    ll <- if ("stats4" %in% loadedNamespaces()) stats4:::logLik else logLik

	# look inside CAIC objects
    if (length(list(...))) {
	        val <- lapply(list(object, ...), function(X) ll(X$mod))
	        val <- as.data.frame(t(sapply(val, function(el) c(attr(el, 
	            "df"), AIC(el, k = k)))))
	        names(val) <- c("df", "AIC")
	        Call <- match.call()
	        Call$k <- NULL
	        row.names(val) <- as.character(Call[-1L])
	        val
	  }
	    else AIC(ll(object$mod), k = k)
}

## ANOVA cares about model types - need to check crunch
anova.caic <- function(object, ...){

	## borrowing from anova.lm
	if(length(list(object, ...)) == 1L){
		# no other objects, no other args (scale and test only make sense for multiple models)
		anova(object$mod)
	} else {
		## pass on - having a second function allows the easy interception of test and scale
		## arguments out of the list of objects
		return(anova.caiclist(object, ...))	
	}
	
	
}

anova.caiclist <- function(object, ..., scale=0, test='F'){

	# need to check that the contrast methods are the same
	objects <- list(object, ...)

	objectsClass <- sapply(objects, class)
	if(! all(objectsClass == 'caic')) stop("anova() on mix of 'caic' and non-'caic' objects.")

	objectsContrMethod <- sapply(objects, attr, 'contr.method')
	if(length(unique(objectsContrMethod)) > 1L) stop("anova() on mixed contrast methods")

	objectsMacroMethod <- sapply(objects, attr, 'macro.method')
	if(length(unique(objectsMacroMethod)) > 1L) stop("anova() on mixed macrocaic methods")

	## OK - now pass the mod parts of those object into anova.lmlist()
	mods <- lapply(objects, '[[', 'mod')
	args <- c(mods, list(scale=scale, test=test))
	anv  <- do.call('anova.lmlist', args)
	# attr(anv, 'heading')[1] <- "Analysis of Variance Table from 'caic' objects\n"
	return(anv)
}
	
	
