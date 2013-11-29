library(ape)

chiroptera.tree  <- read.caic("bats.Phyl")
carnivora.tree   <- read.caic("Carnivora.Phyl", brlen="Carnivora.BLen")
primates.tree    <- read.caic("Primates.Phyl", brlen="Primates.BLen")
marsupialia.tree <- read.caic("Marsupialia.Phyl")

chiroptera.data  <- read.delim("Bats.txt", na.strings="-9", stringsAsFactors=FALSE)
carnivora.data   <- read.delim("Carnivores.txt", na.strings="-9", stringsAsFactors=FALSE)
primates.data    <- read.delim("Primates.txt", na.strings="-9", stringsAsFactors=FALSE)
marsupialia.data <- read.delim("Marsupial.txt", na.strings="-9", stringsAsFactors=FALSE)


# missing variables
primates.data$lendimorph    <- NA
marsupialia.data$lendimorph <- NA
marsupialia.data$grpsize    <- NA
chiroptera.data$popdens     <- NA
chiroptera.data$ibi         <- NA

# order names
primates.data$order.name    <- 'Primates'
marsupialia.data$order.name <- 'Marsupialia'
carnivora.data$order.name   <- 'Carnivora'
chiroptera.data$order.name  <- 'Chiroptera'

# richness values
primates.data$species.rich    <- 1
marsupialia.data$species.rich <- 1
carnivora.data$species.rich   <- 1
chiroptera.data$species.rich  <- 1


# get a common set of variables
primates.data <- primates.data[, c("order.name", "Species.Name", "species.rich", "lnBodyMass", "lnASM", "lnGestation", "lnIBI", "lnLitterSize", "lnPopDens", "lnGrpSize", "Ln.Dimorph.", "lendimorph")]
carnivora.data <- carnivora.data[, c("order.name", "Species.name", "species.rich", "lnBodyMass", "lnASM.d.", "lnGestation", "lnIBI", "lnLitterSize", "lnPopDens", "lnGroupSize", "lnMassDimorph", "lnLenDimorph")]
chiroptera.data  <- chiroptera.data [, c("order.name", "SPP", "species.rich",  "LNABWT", "LNFSEXMA", "LNGEST", "ibi", "LNLITSIZ", "popdens", "LNGRPSIZE", "lnSizeDimorph", "lnLenDimorph")]
marsupialia.data <- marsupialia.data[, c("order.name", "Species.Name.WR.", "species.rich", "lnAvbwt.g.", "lnAvASM.mo.", "lnGestation.d.", "lnIBI.d.", "lnLittSize", "lnPopDens", "grpsize", "ln.Dimorphism.", "lendimorph")]

commonVarNames <- c('order.name','binomial', "species.rich", 'body.mass', 'age.sexual.maturity', 'gestation', 'interbirth.interval', 'litter.size', 'population.density', 'group.size', 'mass.dimorphism', 'length.dimorphism')
names(chiroptera.data ) <- names(carnivora.data) <- names(primates.data) <- names(marsupialia.data) <- commonVarNames

save(chiroptera.tree ,carnivora.tree  ,primates.tree   ,marsupialia.tree,chiroptera.data ,carnivora.data  ,primates.data   ,marsupialia.data,  file='IsaacEtAl.rda')