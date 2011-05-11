benchTreeDicho <- read.tree('BenchTreeDi.tre')
benchTreePoly <- read.tree('BenchTreePoly.tre')
benchData <-read.delim('BenchData.txt')
benchData$triFact <- ordered(benchData$triFact)
benchData$triFactNA <- ordered(benchData$triFactNA)

testTree <- read.nexus('test.nex')
testData <- read.delim('test.dat', header=FALSE)
names(testData) <- c('tips','V1','V2')

save(benchTreeDicho, benchTreePoly, benchData, testTree, testData, file="../data/benchTestInputs.rda")
