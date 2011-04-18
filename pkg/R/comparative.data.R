
comparative.data <- function(phy, data, names.col, vcv=FALSE, vcv.dim=2, na.omit=TRUE){

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
        

    # Compile comparative dataset
	# COULD HAVE phy components as first level slots rather than nested within $phy
	# and have comparative.data inherit methods from phylo, but that presupposes
	# that none of the phylo functions you might use strip out contents
	
    RET <- list(phy=matchedPhy, data = data, 
                data.name=data.name, phy.name=phy.name, 
                dropped=list(tips=to.drop, unmatched.rows=row.not.in.tree))
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

# some useful generics
print.comparative.data <- function(x, ...){

    # basic summary data
    cat("Comparative dataset of", nrow(x$data), "taxa:\n")
    cat("Phylogeny:", x$phy.name, "\n")
    cat("   ", length(x$phy$tip.label), " tips, ", x$phy$Nnode, " internal nodes\n  ", sep='')
    # this is a bit of a hack - can't get str for a vector to take an indent
    str(x$phy$tip.label)
    if(! is.null(x$vcv)){
	    cat('VCV matrix present:\n  ')
	    str(x$vcv, give.attr=FALSE)
	}
    cat("Data:" , x$data.name, "\n")
    str(as.list(x$data), no.list=TRUE, indent='   ')

	# report on mismatch on merge
	dropCount <- sapply(x$dropped, length)
    if(any(dropCount)){
	    cat('Dropped taxa:\n')
	    cat('   ', x$phy.name , ' { ', dropCount[1], ' ( ',nrow(x$data), 
	        ' } ', dropCount[2], ' ) ', x$data.name, sep='')
    }
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
    
		# add to dropped list
		x$dropped$tips <- c(x$dropped$tips, to.drop)
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

"[.comparative.data" <- function(x, i, j) {
	
	# how many args?
	# 2 and missing(i) = x[] --> return x untouched
	# 2 and ! missing(i) = x[i] --> column subset 
	# otherwise 3 = x[i,j] or x[i,] or x[,j] or x[,] 
	if( nargs() ==  2) {
		if(missing(i)){
			return(x)
		} else {
			j <- i
			hasI <- FALSE
		}
	} else { hasI <- TRUE }
	
	# no drop argument permitted. can't lose dataframeness
	if(! missing(j)){
		if(is.null(j)) stop('Null indices not permitted on comparative data objects')
		x$data <- x$data[,j, drop=FALSE]
	}

	# no recycling, no out of index rows
	# no simple reordering possible because the tree implies
	# an order to the data frame
	if(! missing(i) & hasI){
		if(is.null(i)) stop('Null indices not permitted on comparative data objects')
		rownames <- x$phy$tip.label
		if(is.character(i)){
			toKeep <- na.omit(match(i, rownames))
			if(length(toKeep) != length(i)) warning('Some tip names were not found')
			toKeep <- rownames[toKeep]
		} else if (is.numeric(i)) {
			# convert to integer (same as [.data.frame)
			i <- as.integer(i)
			if(all(i > 0)){
				toKeep <- intersect(i, seq_along(rownames))
			} else if(all(i < 0)) {
				toKeep <- setdiff(seq_along(rownames), abs(i))
			} else {
				stop("only 0's may be mixed with negative subscripts")
			}
			toKeep <- rownames[toKeep]
			if(! all(abs(i) %in% seq_along(rownames))) warning('Some row numbers were not found')
		} else if (is.logical(i)) {
			if(length(i) != length(rownames)) stop('Logical index does not match number of tips')
			toKeep <- rownames[i]
		}
		
		# Work out which to drop
		rowToKeep <- match(toKeep, rownames)
		if(length(rowToKeep) < 2) stop('Comparative dataset reduced to fewer than 2 taxa')
		
		toDrop    <- setdiff(rownames, toKeep)
		# this assumes that drop.tip preserves the remaining order - testing suggests ok
		x$phy <- drop.tip(x$phy, toDrop)
        # lose VCV elements if needed
        if(! is.null(x$vcv))  x$vcv <- x$vcv[rowToKeep, rowToKeep]
        x$data <- x$data[rowToKeep,, drop=FALSE]
	}
	return(x)
}

reorder.comparative.data <- function(x, order = "cladewise", ...){
	
	# Uses ape reorder code
	order <- match.arg(order, c("cladewise", "pruningwise"))

	# test for existing order to avoid duplicate calls
    if (!is.null(attr(x, "order"))) 
        if (attr(x, "order") == order) 
            return(x)
	
	# exclude 2 taxon trees
    nb.node <- x$phy$Nnode
    if (nb.node == 1) 
        return(x)
	
	# otherwise
    nb.tip <- length(x$phy$tip.label)
    nb.edge <- dim(x$phy$edge)[1]
    neworder <- if (order == "cladewise") 
        .C("neworder_cladewise", as.integer(nb.tip), as.integer(x$phy$edge[, 
            1]), as.integer(x$phy$edge[, 2]), as.integer(nb.edge), 
            integer(nb.edge), PACKAGE = "ape")[[5]]
    else .C("neworder_pruningwise", as.integer(nb.tip), as.integer(nb.node), 
        as.integer(x$phy$edge[, 1]), as.integer(x$phy$edge[, 2]), as.integer(nb.edge), 
        integer(nb.edge), PACKAGE = "ape")[[6]]
	
	# apply new order to elements
    x$phy$edge <- x$phy$edge[neworder, ]
    attr(x$phy, "order") <- order

    if (!is.null(x$phy$edge.length)) 
        x$phy$edge.length <- x$phy$edge.length[neworder]

	if(! is.null(x$vcv)) x$vcv <- x$vcv[neworder, neworder]
	
	x$dat <- x$dat[neworder,]
    attr(x, "order") <- order
    
	return(x)
}

## x <- comparative.data(shorebird.tree, shorebird.data, 'Species')
## x[]
## x[,]
## x[2:3]
## x[, 2:3]
## x[1:15, ]
## x[1:15, 2:3]

## ## $ method: Don't think this is possible with S3 - want $ to be able
## ## to give back a column from data - without breaking the use
## ## of $ for subsetting. Could use name matching to figure out which
## ## but then duplicate names in data and in the class list are a huge
## ## programming gotcha.

## "$.comparative.data" <- function(x, name) {
## 	
## 	# careful to avoid using $ in here otherwise infinite recursion kicks off.
## 	return(x[['data']][, name])
## 
## }
