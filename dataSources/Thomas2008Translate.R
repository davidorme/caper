BritishBirds.tree <- read.nexus('../benchmarks/Thomas2008.tre')
BritishBirds.data <- read.delim('../benchmarks/Thomas2008Data.txt', stringsAsFactors=FALSE)
save(BritishBirds.tree, BritishBirds.data, file='BritishBirds.rda')
