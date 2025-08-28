#' Datasets used for benchmarking caper
#'
#' These data files and model objects contain input datasets for benchmarking
#' the caper package against other implementations and the results returned by
#' those other implementations. For more details, see the benchmark vignette
#' for caper.
#'
#'
#' @name caper-benchmarks
#' @aliases caper-benchmarks benchBayesTraitsOutputs benchBrunchOutputs
#' benchCrunchOutputs benchFuscoOutputs benchMacroCAICOutputs benchMesaOutputs
#' benchTestInputs BayesTraitsMods benchData benchTreeDicho benchTreePoly
#' CAIC.BrDi1057 CAIC.BrDi1157 CAIC.BrDi813 CAIC.BrDi913 CAIC.BrPl1057
#' CAIC.BrPl1157 CAIC.BrPl813 CAIC.BrPl913 CAIC.CrDi213 CAIC.CrDi657
#' CAIC.CrPl213 CAIC.CrPl413 CAIC.CrPl657 FuscoDiSpp FuscoDiTax FuscoPolySpp
#' FuscoPolyTax MacroCAIC.DiSpp23 MacroCAIC.DiSpp67 MacroCAIC.DiTax23
#' MacroCAIC.DiTax67 MacroCAIC.PolySpp23 MacroCAIC.PolySpp67
#' MacroCAIC.PolyTax23 MacroCAIC.PolyTax67 MeSA.I testData testTree nul fix kld
#' Kld kLd klD KLd kLD KlD KLD
#' @docType data
#' @seealso vignette()
#' @keywords datasets
NULL


#' Conservation status of British birds (Thomas 2008)
#'
#' The dataset contains a molecular phylogeny of 249 species of British birds
#' and a data frame containing information on the conservation status of 181 of
#' those species. The dataset is taken from Thomas (2008) and the data frame
#' contains 26 variables:
#'
#' binomial common_name Red_list Amber_list Green_list Red_amber_list IUCN
#' Red_HD Amber_HD Red_amber_HD Red_list_BDp amb_BDp Red_amb_BDp Red_list_BDr
#' Red_list_WDp amb_BDr amb_WDMp amb_spec amb_BR amb_BL amb_WL amb_BI amb_WI
#' Red_amb_BDr pop_size range_size
#'
#' @name BritishBirds
#' @docType data
#' @references Thomas, G. H. (2008). Phylogenetic distributions of British
#' birds of conservation concern. Proceedings of the Royal Society B-Biological
#' Sciences, 275(1647):2077-2083.
#' @keywords datasets
#' @examples
#' data(BritishBirds)
#' BritishBirds <- comparative.data(
#'     BritishBirds.tree, BritishBirds.data, binomial
#' )
NULL

#' @rdname BritishBirds
"BritishBirds.tree"

#' @rdname BritishBirds
"BritishBirds.data"


#' Example dataset for Fusco imbalance calculations
#'
#' This dataset contains the phylogeny of bird families and species richness
#' originally included with the FUSCO imbalance calculation programs.
#'
#' @name fuscoBird
#' @format The dataset provides a phylogeny (\code{fuscoBirdTree}) and a data
#' frame (\code{fuscoBirdData}). The phylogeny is a 137 tip tree containing one
#' polytomy and the dataframe provides tip labels and species richness for each
#' bird family. This dataset was provided with the original DOS implementation
#' of the test and the families were unlabelled in this original file.
#' @seealso fusco.test
#' @keywords datasets
NULL

#' @rdname fuscoBird
"fuscoBirdData"

#' @rdname fuscoBird
"fuscoBirdTree"


#' Example dataset for the caper package
#'
#' This data set contains four species-level comparative datasets used in Isaac
#' et al (2005).
#'
#' @name IsaacEtAl
#' @format The datafile contains species level phylogenies and accompanying
#' data frames of nine variables for each of four mammalian orders (Primates,
#' Carnivora, Chiroptera and Marsupialia). The data were published in
#' supplementary material for Isaac et al. (2005) as CAIC format files and text
#' data files and have been converted for use in 'caper'. The data files are
#' incomplete, with some variables having little or no data for some orders.
#'
#' The variables (all saved as natural log values) are:
#' \describe{
#'  \item{species.rich}{
#'       Species richness at the tips - all are set to 1 for use in
#'       \code{macrocaic}
#'  }
#'  \item{body.mass}{The average body mass in grams}
#'  \item{age.sexual.maturity}{Age at sexual maturity in months}
#'  \item{gestation}{Gestation length in days}
#'  \item{interbirth.interval}{Interbirth interval in months}
#'  \item{litter.size}{The average number of offspring in a litter}
#'  \item{population.density}{Population density}
#'  \item{group.size}{Number of individuals in a typical group}
#'  \item{mass.dimorphism}{Male mass /female mass}
#'  \item{length.dimorphism}{Male length / female length}
#' }
#' @seealso caic, pgls
#' @references Isaac, N., Jones, K., Gittleman, J., and Purvis, A. (2005).
#' Correlates of species richness in mammals: Body size, life history, and
#' ecology. American Naturalist, 165(5):600-607.
#' @keywords datasets
#' @examples
#'
#' data(IsaacEtAl)
#' chiroptera <- comparative.data(
#'     chiroptera.tree, chiroptera.data, "binomial",
#'     na.omit = FALSE
#' )
#' carnivora <- comparative.data(
#'     carnivora.tree, carnivora.data, "binomial",
#'     na.omit = FALSE
#' )
#' primates <- comparative.data(
#'     primates.tree, primates.data, "binomial",
#'     na.omit = FALSE
#' )
#' marsupialia <- comparative.data(
#'     marsupialia.tree, marsupialia.data, "binomial",
#'     na.omit = FALSE
#' )
NULL

#' @rdname IsaacEtAl
"chiroptera.tree"

#' @rdname IsaacEtAl
"chiroptera.data"

#' @rdname IsaacEtAl
"carnivora.tree"

#' @rdname IsaacEtAl
"carnivora.data"

#' @rdname IsaacEtAl
"primates.tree"

#' @rdname IsaacEtAl
"primates.data"

#' @rdname IsaacEtAl
"marsupialia.tree"

#' @rdname IsaacEtAl
"marsupialia.data"


#' Example dataset for the CAIC package
#'
#' This is a comparative dataset on Perissodactyla taken from the examples
#' include with the original CAIC program.
#'
#' @name perissodactyla
#' @format The datafile contains a phylogeny (\code{perissodactyla.tree}) of 18
#' perissodactyl species as a 'phylo' object from the \code{ape} library. The
#' tip names are the binomial names of the species. The file also contains a
#' data frame (\code{perissodactyla.data}) of variables 5 variables for 13 of
#' those species:
#'
#' \describe{
#'      \item{Binomial}{The species binomial name.}
#'      \item{log.female.wt}{Log female weight}
#'      \item{log.gestation.length}{Log gestation length}
#'      \item{log.neonatal.wt}{Log neonatal weight}
#'      \item{Territoriality}{
#'          A factor indicating whether or not the species displays territorial
#'          behaviour.
#'      }
#' }
#' The dataset is incomplete - it does not include data for each species in the
#' phylogeny and contains missing values. See the examples for the behaviour of
#' the 'comparative.data' function in handling missing data.
#' @seealso caic, pgls
#' @references Purvis, A. and Rambaut, A. (1995). Comparative Analysis by
#' Independent Contrasts (CAIC) User's Guide.
#' @keywords datasets
#' @examples
#'
#' data(perissodactyla)
#' # default behaviour is to omit incomplete data rows
#' perisso <- comparative.data(
#'     perissodactyla.tree, perissodactyla.data, Binomial
#' )
#' # but this can be turned off
#' perisso <- comparative.data(
#'     perissodactyla.tree, perissodactyla.data, Binomial,
#'     na.omit = FALSE
#' )
#' na.omit(perisso)
NULL


#' @rdname perissodactyla
"perissodactyla.data"

#' @rdname perissodactyla
"perissodactyla.tree"



#' Example dataset for the caper package
#'
#' This is a comparative dataset on the evolution of shorebird egg size taken
#' from Lislevand and Thomas (2006).
#'
#' @name shorebird
#' @format The datafile contains a phylogeny (\code{shorebird.tree}) of 71
#' shorebird species as a 'phylo' object from the \code{ape} library. The tip
#' names are the binomial names of the species. The file also contains a data
#' frame (\code{shorebird.data}) of 71 complete cases for those species. The
#' data frame contains six variables:
#'
#' \describe{
#'      \item{Species}{The species binomial name.}
#'      \item{M.Mass}{The adult male mass body mass in grams.}
#'      \item{F.Mass}{The adult female mass body mass in grams.}
#'      \item{Egg.Mass}{The fresh egg mass in grams.}
#'      \item{Cl.size}{The mean clutch size}
#'      \item{Mat.syst}{
#'           The mating system, as a three level factor: monogamous (MO),
#'           polygynous (PO) or polyandrous (PA).
#'      }
#' }
#' @seealso crunch, pgls
#' @references Lislevand, T and Thomas, G. H. (2006) Limited male incubation
#' ability and the evolution of egg size in shorebirds. Biology Letters 2, 206
#' - 208
#' @keywords datasets
NULL


#' @rdname shorebird
"shorebird.data"

#' @rdname shorebird
"shorebird.tree"



#' The syrphidae dataset of Katzourakis et al. 2001
#'
#' A genus level phylogeny of the Syrphidae (sawflies) along with data on the
#' species richness of each genus.
#'
#' @name syrphidae
#' @format A 'phylo' object (\code{syrphidaeTree}) containing a phylogeny of
#' 204 genera of sawflies and a data frame (\code{syrphidaeRich}) of the
#' species richness of each genus.
#' @seealso fusco.test
#' @references Katzourakis, A., Purvis, A., Azmeh, S., Rotheray, G., and
#' Gilbert, F. (2001). Macroevolution of hoverflies (Diptera : Syrphidae): the
#' effect of using higher-level taxa in studies of biodiversity, and correlates
#' of species richness. Journal of Evolutionary Biology, 14:219-227.
#' @keywords datasets
NULL

#' @rdname syrphidae
"syrphidaeTree"

#' @rdname syrphidae
"syrphidaeRich"
