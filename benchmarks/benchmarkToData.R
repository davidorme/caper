## Script to take the benchmark data files and store them in .rda 
## files in the data/ directory. Greatly simplifies access of the benchmark
## data from vignette building and also enforces .Rd documentation

## Not that this isn't a pain in the arse, but it fits the package structure...

## NOTE THAT THE BENCHMARK FOLDER HOLDS MUCH MORE INFORMATION
## INCLUDING THE INPUT FILES TO THE ORIGINAL PROGRAMS
## BUT ONLY THE OUTPUTS AND DATA TO RUN MATCHING 


## SAVE THE CAIC OUTPUT TABLES INTO DATAFRAMES IN AN .RDA FILE

CAIC.CrDi213 <- read.delim("CAIC_outputs/BenchCAIC.dat_CrDi213")
CAIC.CrDi657 <- read.delim("CAIC_outputs/BenchCAIC.dat_CrDi657")
CAIC.CrPl213 <- read.delim("CAIC_outputs/BenchCAIC.dat_CrPl213")
CAIC.CrPl413 <- read.delim("CAIC_outputs/BenchCAIC.dat_CrPl413")
CAIC.CrPl657 <- read.delim("CAIC_outputs/BenchCAIC.dat_CrPl657")

save(CAIC.CrDi213, CAIC.CrDi657, CAIC.CrPl213, CAIC.CrPl413,
     CAIC.CrPl657, file="../data/benchCrunchOutputs.rda")

CAIC.BrDi813  <- read.delim("CAIC_outputs/BenchCAIC.dat_BrDi813")
CAIC.BrDi913  <- read.delim("CAIC_outputs/BenchCAIC.dat_BrDi913")
CAIC.BrDi1057 <- read.delim("CAIC_outputs/BenchCAIC.dat_BrDi1057")
CAIC.BrDi1157 <- read.delim("CAIC_outputs/BenchCAIC.dat_BrDi1157")
CAIC.BrPl813  <- read.delim("CAIC_outputs/BenchCAIC.dat_BrPl813")
CAIC.BrPl913  <- read.delim("CAIC_outputs/BenchCAIC.dat_BrPl913")
CAIC.BrPl1057 <- read.delim("CAIC_outputs/BenchCAIC.dat_BrPl1057")
CAIC.BrPl1157 <- read.delim("CAIC_outputs/BenchCAIC.dat_BrPl1157")

save(CAIC.BrDi813, CAIC.BrDi913, CAIC.BrDi1057, CAIC.BrDi1157, 
     CAIC.BrPl813, CAIC.BrPl913, CAIC.BrPl1057, CAIC.BrPl1157, 
     file="../data/benchBrunchOutputs.rda")

## SAVE THE MACROCAIC OUTPUTS INTO DATAFRAMES IN AN .RDA FILE
                    
MacroCAIC.DiSpp23   <- read.delim("MacroCAIC_outputs/DiSpp23.txt")
MacroCAIC.DiSpp67   <- read.delim("MacroCAIC_outputs/DiSpp67.txt")
MacroCAIC.DiTax23   <- read.delim("MacroCAIC_outputs/DiTax23.txt")
MacroCAIC.DiTax67   <- read.delim("MacroCAIC_outputs/DiTax67.txt")
MacroCAIC.PolySpp23 <- read.delim("MacroCAIC_outputs/PolySpp23.txt")
MacroCAIC.PolySpp67 <- read.delim("MacroCAIC_outputs/PolySpp67.txt")
MacroCAIC.PolyTax23 <- read.delim("MacroCAIC_outputs/PolyTax23.txt")
MacroCAIC.PolyTax67 <- read.delim("MacroCAIC_outputs/PolyTax67.txt")

save(MacroCAIC.DiSpp23, MacroCAIC.DiSpp67, MacroCAIC.DiTax23, MacroCAIC.DiTax67, 
     MacroCAIC.PolySpp23, MacroCAIC.PolySpp67, MacroCAIC.PolyTax23, MacroCAIC.PolyTax67,
     file="../data/benchMacroCAICOutputs.rda")

## IMPORT AND PROCESS THE FUSCO OUTPUT DATA - A HISTOGRAM AND NODE BY NODE TABLE - AND
## STORE THEM INTO AN OBJECT FOR EACH ANALYSIS

readFuscoOut <- function(fname){

	distTab <- read.table(fname, skip=3, nrows=10, sep="|")
	distTab <- distTab[, -c(1,6)]
	colnames(distTab) <- c("imbalanceBin","obsFreq","corr","corrFreq")

	distSum <- as.list(scan(fname, skip=15, n=10, what=list(NULL,1), sep="=")[[2]])
	names(distSum) <- c("nSpp","nTips","nInfNodes","medianI","qdI")

	nodes <- read.table(fname, skip=24, sep="|")
	nodes <- nodes[,-c(1,7)]
	colnames(nodes) <- c("node","nodeT","nodeB","nodeS","nodeI")

	return(list(distTab=distTab, distSum=distSum, nodes=nodes))
}

FuscoDiSpp   <- readFuscoOut("FUSCO_outputs/FuscoDiSpp.txt")
FuscoDiTax   <- readFuscoOut("FUSCO_outputs/FuscoDiTax.txt")
FuscoPolySpp <- readFuscoOut("FUSCO_outputs/FuscoPolySpp.txt")
FuscoPolyTax <- readFuscoOut("FUSCO_outputs/FuscoPolyTax.txt")

save(FuscoDiSpp, FuscoDiTax, FuscoPolySpp, FuscoPolyTax,
     file="../data/benchFuscoOutputs.rda")

## IMPORT AND PROCESS THE MESA FUSCO ANALYSIS

MeSA.I <- scan("MeSA_outputs/MeSA_FuscoI.txt", sep="\n", what="character")
MeSA.I <- strsplit(MeSA.I, split="\t")
MeSA.I <- lapply(MeSA.I, function(X)X[-1:-2])
MeSA.I <- as.data.frame(MeSA.I, stringsAsFactors=FALSE)[,c(-5,-7)]
names(MeSA.I) <- c("clade","I","weight","nodeSize","Iprime")
MeSA.I$I <- as.numeric(MeSA.I$I)
MeSA.I$weight <- as.numeric(MeSA.I$weight)
MeSA.I$nodeSize <- as.numeric(MeSA.I$nodeSize)
MeSA.I$Iprime <- as.numeric(MeSA.I$Iprime)

save(MeSA.I, file="../data/benchMesaOutputs.rda")

## NOW CREATE DATAFILES OF THE INPUTS

BenchData <- read.delim("BenchData.txt")
BenchData$triFact <- ordered(BenchData$triFact, levels=c("A","B","C"))
BenchData$triFactNA <- ordered(BenchData$triFactNA, levels=c("A","B","C"))


BENCH     <- read.tree("BenchTreeDi.tre")
BENCHPoly <- read.tree("BenchTreePoly.tre")

save(BenchData, BENCH, BENCHPoly, file="../data/benchTestInputs.rda")

syrphidaeTree <- read.nexus("syrphidae.nexus")
syrphidaeRich <- read.delim("syrphidae_tabbed.txt", header=FALSE, col.names=c("genus","nSpp"))

# strip off the species numbers from the genus names
syrphidaeTree$tip.label <- gsub("^[0-9]+", "", syrphidaeTree$tip.label)
syrphidaeRich$genus <- gsub("^[0-9]+", "", syrphidaeRich$genus)

save(syrphidaeTree, syrphidaeRich, file="../data/syrphidae.rda")



