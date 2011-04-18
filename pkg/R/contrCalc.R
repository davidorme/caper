`contrCalc` <-
function(vals, phy, ref.var, picMethod, crunch.brlen){

    # Takes the tip values and analysis tree and calculates matrices 
    # of contrasts and nodal values.

    # if picMethod is a macro analysis, then ref.var is a vector of nodal values,
    # otherwise it is a reference to a column in vals. Either way, these are used to
    # determine the direction of the calculation of contrasts
  
    # DESIGN THOUGHTS - simpler to maintain with a single function for each type of contrast
    #                 - create an initialize function to do these shared first steps?
    
    Root <-  with(phy, (max(edge)-Nnode)+1)
    IntNd <- Root:max(phy$edge) 
    nIntNd <- phy$Nnode

    contrMat <- matrix(NA, ncol=dim(vals)[2], nrow=nIntNd, dimnames=list(IntNd, dimnames(vals)[[2]]))
    nodVal <- rbind(vals, contrMat)
    
    # node depth vector
    nodeDepth <- rep(c(1,NA), times=c(length(phy$tip.label), nIntNd))
    names(nodeDepth) <- rownames(nodVal)
    
    # vector of number of children and variance used in calculation
    nChild <- numeric(nIntNd)
    names(nChild) <- IntNd

    contrVar <- numeric(nIntNd)
    names(contrVar) <- IntNd

    # for use in brunch analyses, a vector to note whether an
    # internal node has been used or not
    brunchUsed <- logical(dim(nodVal)[1])
    names(brunchUsed) <- rownames(nodVal)

    # identify the groups contributing to each contrast - split() sorts into numeric order
	contrGp <- split(phy$edge[, 2], as.numeric(phy$edge[, 1]))
	# put into a pruningwise order
	phy <- reorder(phy, 'pruningwise')
	contrGp <- contrGp[as.character(unique(phy$edge[,1]))]
    
    # loop the nodes
    for(nd in seq_along(contrGp)){
        
        # ID the parents and children
        parent <- names(contrGp)[nd]
        children <- contrGp[[nd]]
        bl <- with(phy, edge.length[match(children, edge[,2])])
        vals <- nodVal[children,, drop=FALSE]
        
        # find complete cases at node... NB: Here, need to use a single set
        # of children nodes in order to have a consistent set of branch lengths
        # used in calculations - must use complete cases rather than dealing with
        # the available data for each variable
        compChild <- complete.cases(vals)
        nChild[parent] <- sum(compChild)
        
        switch(picMethod,  
            "PDI" = { 
                     # allow for missing data in the reference variable (can't make a contrast)
                     # - don't need to check complete cases because complete spp rich data is mandatory
                    if( any(is.na(ref.var[children]))){
                        currContr <-NA
                    } else {
                    if(dim(vals)[1] > 2 ){ # polytomy
                        currContr <- NA
                    } else {
                        vals <- vals[order(ref.var[children], decreasing=TRUE), , drop=FALSE]
                        currContr <- (vals[1,]/colSums(vals)) - 0.5
                    }}
                    currNV <- colSums(vals)
                    currVar <- sum(bl)
                    currBlAdj <- prod(bl)/sum(bl)
                    },
            "RRD" = {
                    # allow for missing data in the reference variable (can't make a contrast)
                    if( any(is.na(ref.var[children]))){
                        currContr <-NA
                    } else {
                    if(dim(vals)[1] > 2 ){
                         currContr <- NA
                    } else {
                        vals <- vals[order(ref.var[children], decreasing=TRUE), , drop=FALSE]
                        currContr <- log(vals[1,]/vals[2,])
                    }}
                    currNV <- colSums(vals)
                    currVar <- sum(bl)
                    currBlAdj <- prod(bl)/sum(bl)
                    },
        "crunch" = {
                    # get the values of the reference variable
                    # - ref.var identifies the column
                    rv <- vals[, ref.var]
                    
                    # can't do anything with _no_ data, but need to do something with
                    # one or more children with information...                     
                    if(any(compChild)){
                        
                        # continue calculation with those complete rows
                        compVals <- vals[compChild,, drop=FALSE]
                        bl <- bl[compChild]
                        rv <- rv[compChild]
                         
                        if(sum(compChild) == 1){
                            # only one complete row... pass through to next node down
                            currNV <- compVals
                            currContr <- NA
                            currBlAdj <- bl
                            currVar <- NA
                        }
                            
                        if(sum(compChild) == 2){
                            # dichotomous node
                            currNV <- colSums(compVals*(1/bl))/sum(1/bl)
                            currContr <- diff(compVals[order(rv),])
                            currBlAdj <- prod(bl)/sum(bl)
                            currVar <- sum(bl)
                        }
                        
                        if(sum(compChild) > 2){
                            # This is the pagel method for polytomies 
                            ## think about this - check out Nick Isaac's comment about only stalling on negative variance at a node.                            
                            bl <- bl - crunch.brlen
                            if(any(bl < 0)) stop("Crunch contrast calculation at a polytomy gives negative branch lengths.")
                            
                            # is there any (meaningful) variance in the reference variable?
                           if(var(rv) < .Machine$double.eps){ 
                               # compare first to the rest (as in CAIC)
                                group <- c(TRUE, rep(FALSE, length(rv) -1))
                           } else {
                                # find groupings a vector indicating whether each row 
                                # is bigger or smaller than the mean of reference variable or is NA
                                group <- (rv > mean(rv, na.rm=TRUE)) # TODO - think >= or >?                               
                           }
                     
                            ProdValBl <- aggregate(compVals * (1/bl), by=list(group), FUN=sum)[,-1]
                            SumWght <- aggregate((1/bl), by=list(group), FUN=sum)[,-1]
                            subNV <- as.matrix(ProdValBl/SumWght)
                            subBL <- crunch.brlen + (1 / SumWght)
                     
                            currContr <- diff(subNV)
                            currNV <- colSums(subNV*(1/subBL))/sum(1/subBL) # weighted means
                            currVar <- sum(subBL)
                            currBlAdj <- 1/(sum(1/subBL))
                        }
                    } else {
                        # no complete cases so pass NA down to the parent node
                        currContr <- NA
                        currNV <- NA
                        currBlAdj <- 0
                        currVar <- NA
                    }},
        "brunch" = {

                    # get the values of the reference variable
                    # - ref.var identifies the column
                    rv <- vals[, ref.var]
                    
                    # further exclude any nodes which have been used to calculate contrasts
                    compChild <- compChild & ! brunchUsed[children]
                    
                    # can't do anything with _no_ data, but need to do something with
                    # one or more children with information...                     
                    if(any(compChild)){
                        
                        # continue calculation with those complete rows
                        compVals <- vals[compChild,, drop=FALSE]
                        bl <- bl[compChild]
                        rv <- rv[compChild]
                         
                        if(sum(compChild) == 1){
                            # only one complete row... pass through to next node down as with crunch
                            currNV <- compVals
                            currContr <- NA
                            currBlAdj <- bl
                            currVar <- NA
                        }
                            
                        if(sum(compChild) == 2){
                            # dichotomous node
                            # are there differences in the reference variable? 
                            # if there are, make a contrast and then pass NA vals to the nodal values 
                            # otherwise don't make a contrast and pass a weighted average...
                            if(var(rv) == 0){
                                currContr <- NA
                            } else {
                                currContr <- diff(compVals[order(rv),] )
                                brunchUsed[parent] <- TRUE
                            }
                            
                            currNV <- colSums(compVals*(1/bl))/sum(1/bl)
                            currBlAdj <- prod(bl)/sum(bl)
                            currVar <- sum(bl)
                            
                        }
                        
                        if(sum(compChild) > 2){
                            if(var(rv) == 0){ # weighted average of data
                                currContr <- NA
                                currNV <- colSums(compVals*(1/bl))/sum(1/bl)
                                currBlAdj <- 1/(sum(1/bl))
                            } else {
                                # This is the pagel method for polytomies 
                                
                                ## think about this - check out Nick Isaac's comment about only stalling on negative variance at a node.
                                bl <- bl - crunch.brlen
                                if(any(bl <= 0)) stop("Brunch contrast calculation at a polytomy gives negative or zero branch lengths.")
                                
                                
                                # find groupings a vector indicating whether each row 
                                # is bigger or smaller than the mean of reference variable or is NA
                                # TODO - Hmm. what to do if there is no variance in the reference variable
                                group <- (rv > mean(rv, na.rm=TRUE)) # From CAIC, this looks like gt not geq?
                        
                                ProdValBl <- aggregate(compVals * (1/bl), by=list(group), FUN=sum)[,-1]
                                SumWght <- aggregate((1/bl), by=list(group), FUN=sum)[,-1]
                                subNV <- as.matrix(ProdValBl/SumWght)
                                subBL <- crunch.brlen + (1 / SumWght)
 
                                currContr <- diff(subNV)
                                currNV <- colSums(subNV*(1/subBL))/sum(1/subBL) # weighted means
                                currVar <- sum(subBL)
                                currBlAdj <- 1/(sum(1/subBL))
                                brunchUsed[parent] <- TRUE
                            }
                        }
                    } else {
                        # no complete cases so pass NA down to the parent node
                        currContr <- NA
                        currNV <- NA
                        currBlAdj <- 0
                        currVar <- NA
                    }},
        "piclm" = {

                    # get the values of the reference variable
                    # - ref.var identifies the column
                    rv <- vals[, ref.var]
                    
                    # further exclude any nodes which have been used to calculate contrasts
                    compChild <- compChild & ! brunchUsed[children]
                    
                    # can't do anything with _no_ data, but need to do something with
                    # one or more children with information...                     
                    if(any(compChild)){
                        
                        # continue calculation with those complete rows
                        compVals <- vals[compChild,, drop=FALSE]
                        bl <- bl[compChild]
                        rv <- rv[compChild]
                         
                        if(sum(compChild) == 1){
                            # only one complete row... pass through to next node down as with crunch
                            currNV <- compVals
                            currContr <- NA
                            currBlAdj <- bl
                            currVar <- NA
                        }
                            
                        if(sum(compChild) == 2){
                            # dichotomous node
                            # are there differences in the reference variable? 
                            # if there are, make a contrast and then pass NA vals to the nodal values 
                            # otherwise don't make a contrast and pass a weighted average...
                            if(var(rv) == 0){
                                currContr <- NA
                            } else {
                                currContr <- diff(compVals[order(rv),] )
                                brunchUsed[parent] <- TRUE
                            }
                            
                            currNV <- colSums(compVals*(1/bl))/sum(1/bl)
                            currBlAdj <- prod(bl)/sum(bl)
                            currVar <- sum(bl)
                            
                        }
                        
                        if(sum(compChild) > 2){
                            if(var(rv) == 0){ # weighted average of data
                                currContr <- NA
                                currNV <- colSums(compVals*(1/bl))/sum(1/bl)
                                currBlAdj <- 1/(sum(1/bl))
                            } else {
                                # This is the pic.lm method for polytomies 
                                # randomly resolve the polytomy with 0 length branches 
                                # calculate contrasts at each node and downweight in the regression
                                # possibly applying a correction in the degrees of freedom.
                                browser()
                                nTips <- sum(compChild)
                                randRes <- rtree(nTips)
                                randOrd <- sample(nTips)
                                vals <- vals[randOrd,]
                                rownames(vals) <- 1:nTips
                                bl <- bl[randOrd]
                                randRes$edge.length <- rep(0, 2*nTips-2)
                                randRes$edge.length[match(1:nTips, randRes$edge[,2])] <- bl
                            }
                        }
                    } else {
                        # no complete cases so pass NA down to the parent node
                        currContr <- NA
                        currNV <- NA
                        currBlAdj <- 0
                        currVar <- NA
                    }})
 
        # Insert nodal values back into the table to be used in more nested nodes
        nodVal[parent,] <- currNV
        contrMat[parent,] <- currContr
        contrVar[parent] <- currVar

       # Adjust the parent branch length
        parInd <- with(phy, match(parent, edge[,2]))
        if(! parent == Root){
                phy$edge.length[parInd] <- phy$edge.length[parInd] + currBlAdj}
                
       # track the node depth of  nodes with contrasts
       if(sum(compChild) < 2){ # i.e. a node with data for one tip being passed through or with no data
           nodeDepth[parent] <- max(nodeDepth[children])
       } else { # a contrast has been calculated and this node has a greater depth
           nodeDepth[parent] <- max(nodeDepth[children]) + 1 
       }
    }

    RET <- list(contrMat=contrMat, nodVal=nodVal, var.contr=contrVar, nChild=nChild, nodeDepth=nodeDepth) 
    attr(RET, "contr.type") <- picMethod

    return(RET)
}

