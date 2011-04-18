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
    # simple code to create a table of the contrasts from the caic object
    
        nodeNum <- matrix(as.numeric(names(caicObj$contrast.data$contrVar)),
                          ncol=1, dimnames=list(NULL, "nodeNumber"))
        contr <- with(caicObj$contrast.data$contr, cbind(response, explanatory))
        # colnames(contr) <- paste("C_", colnames(contr), sep="")
        if(nodalValues){
            nv <- with(caicObj$contrast.data$nodalVals, cbind(response, explanatory))
            colnames(nv) <- paste("nodal.", colnames(nv), sep="")
            tab <- as.data.frame(cbind(nodeNum, contr, nv))
        } else {
            tab <- as.data.frame(cbind(nodeNum, contr))
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
        stRes <- rstudent(caicObj$mod)
        tab$studResid <- NA
        tab$studResid[match(as.numeric(names(stRes)), tab$nodeNumber)] <- stRes
        if(validNodes) tab <- subset(tab, validNodes, select=-validNodes)
       
       if(CAIC.codes){
           Cphy <- caic.label(caicObj$data$phy)
           tab$CAIC.code <- Cphy$edge.caic.code[match(tab$nodeNumber, names(Cphy$edge.caic.code))]
       }
       
       return(tab)

}

caicStyleArgs <- function(data, phy, names.col){
	
	# general function to handle old style non-'comparative.data' use of 
	# crunch, brunch, macrocaic functions 
	
	if(missing(data)) stop('data object is missing')
	if(missing(phy)) stop('phy object is missing')

	# check the classes (and allow for them being in the wrong order)
	args <- list(data, phy)
	argClass <- sapply(args, class)
	
	# bail back to calling function if we don't have targets
	# i.e. precisely a phylogeny and a data.frame
	targets <- c('data.frame', 'phylo')
	if(! identical( sort(intersect( argClass, targets)), targets)){
		return(NULL)
	}
	
	# try and build a comparative data object
    if(argClass[1] == 'data.frame'){
		data <- eval(substitute(comparative.data(phy = phy, data = data, names.col = XXX), list(XXX=names.col)))
	} else {
		data <- eval(substitute(comparative.data(phy = data, data = phy, names.col = XXX), list(XXX=names.col)))
	}
	return(data)
	
}

## ## THIS NEEDS SOME WORK TO PASS ALL THESE
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, phy=shorebird.tree, data=shorebird.data, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.tree, shorebird.data, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.data, shorebird.tree, names.col=Species)
## ## BREAKS WITH THE FOLLOWING MISSING ARGUMENTS BUT THESE WOULD HAVE BROKEN ANYWAY
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, phy=shorebird.tree, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, data=shorebird.data, names.col=Species) 
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.tree, names.col=Species)
## crunchMod <- crunch(Egg.Mass ~ F.Mass + M.Mass, shorebird.data, names.col=Species)

summary.caic <- function(object, ...){

    summary(object$mod, ...)

}

