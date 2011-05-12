# code to create a benchmark dataset for use in testing code against
# existing implementations

library(CAIC)

# pure birth model with three continuous traits, one of which depends on the other two
# a binary trait and an ordered three level (can't go from A to C directly)

contTraits <- c(5,9,12)
names(contTraits) <- c("contResp","contExp1","contExp2")
contVar <- matrix(c(1,0.5,0.25,0.5,1,0,0.25,0,1), ncol=3)
discTraits <- list(binFact=matrix(c(0,0.05,0.05,0), ncol=2, dimnames=list(c("A","B"), c("A","B"))),
                   triFact=matrix(c(0,0.05,0,0.05,0,0.05,0,0.05,0), ncol=3, dimnames=list(c("A","B","C"), c("A","B","C"))))

# get a simulation which avoids  short branch lengths at polytomies
discard <- TRUE

while(discard){
    BENCH <- growTree(b=0.1, halt=200,
                      ct.start=contTraits, ct.change=c(0,0,0), ct.var=contVar,
                      dt.rates=discTraits, extend.proportion=1, grain=Inf)


    # make seven-ish internal polytomies...
    collapseNodes <- with( BENCH, sort(edge.length[edge[,2] > Nnode +1]))[7]
    BENCHPoly <- align.tips(di2multi(BENCH, collapseNodes))

    # check that none of the polytomies have descendents on a branch less than 1
    nChild <- table(BENCHPoly$edge[,1])
    polyParent <- as.numeric(names(nChild)[nChild > 2])
    polyDescBL <- BENCHPoly$edge.length[which(BENCHPoly$edge[,1] %in% polyParent)]
    discard <- any(polyDescBL <= 1)
    print(polyDescBL)
    cat("looping\n")
}


# basic data types 
# - continuous response
# - two continuous explanatory
# - integer response for species richness contrasts using a broken stick
# - binary factor
# - three-level factor

BenchData <- with(BENCH, cbind(ct.data, dt.data[,-1])) # drop the node labels from dt.data
BenchData <- BenchData[BenchData$node <= 200,] # drop the internal node data

# - variable that has no variance at some polytomies of the polytomous tree
#   simply by making all descendents of a polytomy share the same value
#   which is the mean at the node of the values for all descendants (_not_ phylogenetic mean!)
#   Avoid the polytomies that are too basal...

BenchData$contExp1NoVar <- BenchData$contExp1
BenchData <- BenchData[order(BenchData$node),]

rownames(BenchData) <- BenchData$node
nodeVals <- clade.stats(BenchData[,"contExp1", drop=FALSE], BENCHPoly, mean)
Desc <- clade.members.list(BENCHPoly)
nDesc <- sapply(Desc, length)
nDesc <- nDesc[as.character(polyParent)]
nDesc <- nDesc[nDesc < 20]

for(nd in seq(along=nDesc)){
    currDesc <- Desc[[names(nDesc)[nd]]]
    BenchData$contExp1NoVar[currDesc] <- nodeVals[names(nDesc)[nd], "contExp1"]    
}

# Testing handling of missing data for each variable - random 5% in each
BenchData$contRespNA <- ifelse(runif(200) < 0.05, NA, BenchData$contResp)
BenchData$contExp1NA <- ifelse(runif(200) < 0.05, NA, BenchData$contExp1)
BenchData$contExp2NA <- ifelse(runif(200) < 0.05, NA, BenchData$contExp2)
BenchData$binFactNA <- BenchData$binFact
BenchData$binFactNA[runif(200) < 0.05] <- NA
BenchData$triFactNA <- BenchData$triFact
BenchData$triFactNA[runif(200) < 0.05] <- NA

# CAIC wants categoricals as the last columns in the file.
BenchData <- BenchData[,c(1:4, 7:10, 5:6, 11:12)]

# CAIC also needs alphabetic order..
CAICData <- BenchData[order(as.character(BenchData$node)),]

# and  factors as numerics ..
CAICData$binFact <- unclass(CAICData$binFact)
CAICData$triFact <- unclass(CAICData$triFact)
CAICData$binFactNA <- unclass(CAICData$binFactNA)
CAICData$triFactNA <- unclass(CAICData$triFactNA)

write.table(CAICData, file="BenchCAIC.dat", quote=FALSE, sep="\t", row.names=FALSE, na="-9")

# make sure the three level factor is ordered
BenchData$triFact <- ordered(BenchData$triFact, levels=c("A","B","C"))
BenchData$triFactNA <- ordered(BenchData$triFactNA, levels=c("A","B","C"))


# output the trees

write.tree(BENCH, file="BenchTreeDi.tre")
write.tree(BENCHPoly, file="BenchTreePoly.tre")

write.caic(BENCH, filebase="BenchTreeDi")
write.caic(BENCHPoly, filebase="BenchTreePoly")


# now create a file of species richnesses at tips for testing MacroCAIC
# same format and order as the .phyl but with nSpp in place of CAIC code.
BenchRich <- scan("BenchTreeDi.phyl", what="character")

# no pattern with a variable, but at least follows the expected broken stick!
brks <- sample(1:4999, 199)
nSpp <- diff(c(0, sort(brks),5000))

BenchRich[seq(1,399, by=2)] <- nSpp
cat(BenchRich, file="BenchTreeDi.rich", sep="\n")
cat(BenchRich, file="BenchTreePoly.rich", sep="\n")

# now insert into the data frame for R
BenchRich <- matrix(BenchRich, ncol=2, byrow=TRUE)
BenchRich <- BenchRich[match(BenchData$node, BenchRich[,2]),]
BenchData$SppRich <- as.numeric(BenchRich[,1])

write.table(BenchData, file="BenchData.txt", quote=FALSE, sep="\t", row.names=FALSE)

# output files for analysis in Fusco

# match species richnesses into the tree string
FuscoDi <- write.tree(BENCH)
FuscoDi <- gsub(":[0-9\\.]+", "", FuscoDi) # get rid of the branch lengths
FuscoDi <- gsub("([0-9]+)","!\\1!", FuscoDi) # put some split characters around the nodes
FuscoDi <- strsplit(FuscoDi, split="!")[[1]] # split the string into node number and topology
nodeMatch <- match(FuscoDi, BenchData$node)
FuscoDi <- ifelse(is.na(nodeMatch), FuscoDi, BenchData$SppRich[nodeMatch]) # switch richness with node number, keep topology
FuscoDi <- paste(FuscoDi, sep="", collapse="") # and glue back together again

FuscoDi <- sub(";", "*", FuscoDi) # switch ; to *
# remove commas except the one between tips
FuscoDi <- gsub("([0-9]),([0-9])", "\\1!\\2", FuscoDi)
FuscoDi <- gsub(",", "", FuscoDi)
FuscoDi <- gsub("!", ",", FuscoDi)
cat("benchdi    ",FuscoDi, "\r\n", file="BenFusco.txt")

# match species richnesses into the tree string
FuscoPoly <- write.tree(BENCHPoly)
FuscoPoly <- gsub(":[0-9\\.]+", "", FuscoPoly) # get rid of the branch lengths
FuscoPoly <- gsub("([0-9]+)","!\\1!", FuscoPoly) # put some split characters around the nodes
FuscoPoly <- strsplit(FuscoPoly, split="!")[[1]] # split the string into node number and topology
nodeMatch <- match(FuscoPoly, BenchData$node)
FuscoPoly <- ifelse(is.na(nodeMatch), FuscoPoly, BenchData$SppRich[nodeMatch]) # switch richness with node number, keep topology
FuscoPoly <- paste(FuscoPoly, sep="", collapse="") # and glue back together again

FuscoPoly <- sub(";", "*", FuscoPoly) # switch ; to *
# remove commas by brackets
FuscoPoly <- gsub("(\\)),", "\\1", FuscoPoly)
FuscoPoly <- gsub(",(\\()", "\\1", FuscoPoly)

cat("benchpoly    ",FuscoPoly, "\r\n", file="BenFusco.txt", append=TRUE)
