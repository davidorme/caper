comparative.data <- function(data, phy, names.col, vcv=FALSE, vcv.dim=2, na.omit=TRUE){

    # TODO - is something odd happening with missing arguments?
    
    # record call and dataset names
    phy.name <- deparse(substitute(phy))
    data.name <- deparse(substitute(data))
    names.col <- as.character(substitute(names.col))

    # check inputs are what they should be
    # DATA:
    
        # ...is a dataframe
        if(! is.data.frame(data)) stop("'data' must be an object of class 'data.frame'.")
        # ...contains the name column and make sure it is of mode character
        namesInd <- match(names.col, names(data))
        if(is.na(namesInd)) {
            stop("Names column '",  names.col, "' not found in data frame '", data.name, "'")
        }
        rownames(data) <- as.character(data[,namesInd])
        # drop the names column
        data <- data[,-namesInd, drop=FALSE]
        
    # PHYLOGENY:
        # check the phylogeny is a rooted phylogeny and set branch lengths...
        if(! inherits(phy, "phylo")) 
            stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
        if(! is.rooted(phy))
            stop("'", deparse(substitute(phy)), "' is not rooted.")

    # MERGE
    
        # store original dataset size
        origTips <- with(phy, max(edge) - Nnode)
        origData <- nrow(data)
        
        # find the intersection between tip.labels and names in data frame
        in.both <- intersect(rownames(data), phy$tip.label)
        if(length(in.both) < 2) stop("Fewer than two tips are common to the dataset and phylogeny")
        
        # TODO: work out what is needed to get rid of the node.label backreferencing
        #       or if it is still needed. Ugly.
        # Label the internal nodes by their node number in the original tree to provide a backreference
        phy$node.label <- with(phy, ((max(edge)-Nnode) +1):max(edge)) 

        # i >> ditch rows with no tip
        row.in.tree <- match(rownames(data), in.both)
        row.not.in.tree <- rownames(data)[is.na(row.in.tree)]
        data <- subset(data, !is.na(row.in.tree))
    
        # ii >> ditch tips which have no rows.
        tip.in.data <-  match(phy$tip.label, in.both)
        to.drop <- phy$tip.label[is.na(tip.in.data)]
        
        #  get subset of phylogeny to be used
        if(length(to.drop) > 0) matchedPhy <- drop.tip(phy, to.drop) else matchedPhy <- phy
        
        # useful info...
        root <- with(matchedPhy, (max(edge) - Nnode) + 1)

        # get the data into the same order as the tips
        tip.order <- match(matchedPhy$tip.label, rownames(data))
        if(any(is.na(tip.order))) stop("Problem with sorting data frame: mismatch between tip labels and data frame labels")
        data <- data[tip.order,, drop=FALSE]
        
        # Label the data frame rows by tip number to allow the tree to be traversed
        # TODO: check on this - requirement of crunch functions.
        # TODO: CLASH WITH PGLM!!!
        rownames(data) <- 1:dim(data)[1]
        
        # Size of conjunction of tree and dataset
        unionData <- dim(data)[1]

    # Compile comparative dataset
    RET <- list(phy=matchedPhy, data=data, N=unionData, data.name=data.name, 
                phy.name=phy.name, na.omit=na.omit, 
                dropped=list(unmatched.rows=row.not.in.tree, tips=to.drop))
    class(RET) <- 'comparative.data'
    
    # Add a VCV array if requested
    if(vcv) {
        RET$vcv <- vcv.array(matchedPhy, dim=vcv.dim)
        RET$vcv.dim <- vcv.dim
    }
    
    # NA handling
    if(na.omit){
        RET <- na.omit(RET) 
    }
    

    return(RET)
}

print.comparative.data <- function(x, ...){

    # basic summary data
    cat("Comparative dataset of ", x$N, " taxa:\n")
    cat("Phylogeny: \n")
    cat("   Source:", x$phy.name, "\n")
    cat("   ", length(x$phy$tip.label), " tips, ", x$phy$Nnode, " internal nodes\n", sep='')
    cat("   Tip names:\n  ") # this is a bit of a hack - can't get str for a vector to take an indent
    str(shorebird.tree$tip.label)
    if(! is.null(x$vcv)){
	    cat('VCV matrix present:\n  ')
	    str(x$vcv, give.attr=FALSE)
	}
    cat("Data: \n")
    str(as.list(x$data), no.list=TRUE, indent='   ')

}

na.omit.comparative.data <- function(x, ...){

    # strips data rows, tips and vcv row/columns for a comparative.data object
    to.drop <- which(! complete.cases(x$data))
    
    # guard against ape 'feature' of dropping all tips for empty vectors
    if(length(to.drop) > 0){
        
        # lose bits of tree
        x$phy <- drop.tip(x$phy, to.drop) 
        # lose rows
        x$data <- x$data[-to.drop,]
        
        # lose VCV elements if needed
        if(! is.null(x$vcv)){
            x$vcv <- x$vcv[-to.drop, -to.drop]
        }
    }
    
    return(x)
}

subset.comparative.data <- function(x, subset, select,  ...){

    ## ripping out the innards of subset.data.frame
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, x$data, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    if (missing(select)) 
        vars <- TRUE
    else {
        nl <- as.list(seq_along(x$data))
        names(nl) <- names(x$data)
        vars <- eval(substitute(select), nl, parent.frame())
    }
    
    ## TODO  - institute tip name subset.
    
    ## now know which rows and columns to keep in the data frame 
    ## guard against ape 'feature' of dropping all tips for empty vectors
    if(any(! r)){
        
        to.drop <- which(! r)
        # lose bits of tree
        x$phy <- drop.tip(x$phy, to.drop) 
        # lose rows
        x$data <- x$data[r, vars, drop = FALSE] # CANNOT lose data.frame-ness
        
        # lose VCV elements if needed
        if(! is.null(x$vcv)){
            x$vcv <- x$vcv[r, r]
        }
    }
    
    return(x)
}


# would be good to have a [ function but may need S4 to do this now
# setMethod("[", signature(x="comparative.data", i="ANY",j="ANY"),
# 
#   function(x, i, j, ..., drop=TRUE) {
#       ## do whatever you want here; your class is not derivedfrom a list
#       ## so we cannot use NextMethod
#       print('argh')
#   })
