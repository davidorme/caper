benchTreeDicho <- read.tree('BenchTreeDi.tre')
benchTreePoly <- read.tree('BenchTreePoly.tre')
benchData <-read.delim('BenchData.txt')
testTree <- read.nexus('test.nex')
testData <- read.delim('test.dat')

save(benchTreeDicho, benchTreePoly, benchData, testTree, testData, file="../data/benchTestInputs.rda")