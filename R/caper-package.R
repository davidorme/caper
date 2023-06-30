#' Comparative analysis of phylogenetics and evolution in R (caper)
#'
#' The \pkg{caper} package provides a set of functions to conduct comparative
#' analyses using both independent contrasts and phylogenetic generalised least
#' squares methods. It also provides functions to combine phylogeny objects
#' with data frames into simple comparative datasets and some utility functions
#' for manipulating phylogenies.
#'
#' In addition to linear models correcting for phylogeny, the package also
#' provides functions to test tree imbalance, calculate various measures of
#' phylogenetic diversity and simulate phylogenies and traits.
#'
#' \tabular{ll}{ Package: \tab caper\cr Type: \tab Package\cr Version: \tab
#' 0.2\cr Date: \tab 2011-05-06\cr License: \tab GPL\cr } A package to carry
#' out phylogenetic comparative analysis.
#'
#' @name caper-package
#' @aliases caper-package caper
#' @docType package
#' @author David Orme Maintainer: David Orme <d.orme@@imperial.ac.uk>
#' @seealso \code{\link[ape:ape-package]{ape}}
#' @keywords package
NULL

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
#' those species. The dataset is taken from Thomas (2008)
#'
#' The data frame contains 26 variables:
#'
#' binomial common_name Red_list Amber_list Green_list Red_amber_list IUCN
#' Red_HD Amber_HD Red_amber_HD Red_list_BDp amb_BDp Red_amb_BDp Red_list_BDr
#' Red_list_WDp amb_BDr amb_WDMp amb_spec amb_BR amb_BL amb_WL amb_BI amb_WI
#' Red_amb_BDr pop_size range_size
#'
#' @name BritishBirds
#' @aliases BritishBirds BritishBirds.tree BritishBirds.data
#' @docType data
#' @references Thomas, G. H. (2008). Phylogenetic distributions of british
#' birds of conservation concern. Proceedings of the Royal Society B-Biological
#' Sciences, 275(1647):2077-2083.
#' @keywords datasets
#' @examples
#' data(BritishBirds)
#' BritishBirds <- comparative.data(
#'     BritishBirds.tree, BritishBirds.data, binomial
#' )
NULL


#' The 'caic' S3 object class and methods
#'
#' The functions 'crunch', 'brunch', 'macrocaic' and 'piclm' all return an
#' object containing an independent contrast model. The structure of the object
#' and the available methods are described here.
#'
#'
#' @format A 'caic' object is a list containing the following: \describe{
#' \item{contrast.data}{ A list of the following: \describe{ \item{contr}{A
#' list containing matrices of the contrasts in the response variable
#' (contr\$response) and explanatory variables (contr\$explanatory).}
#' \item{nodalVals}{A list containing matrices of the nodal values in the
#' response variable (contr\$response) and explanatory variables
#' (contr\$explanatory).} \item{contrVar}{A numeric vector of the expected
#' variance for each contrast.} \item{nChild}{A vector showing the number of
#' nodes descending from each internal node} \item{nodeDepth}{A vector showing
#' the maximum number of nodes between each internal node and the tips of the
#' phylogeny (including both the node in question and the tip and hence always
#' >=2) } \item{validNodes}{A logical vector showing which internal nodes on
#' the tree have valid contrasts, given the available data and any user
#' constraints.}}} \item{data}{A 'comparative.data' object containing the
#' phylogeny used to calculate contrasts and the original data.} \item{lm}{An
#' 'lm' object containing a regression model through the origin for the
#' calculated contrast} }
#'
#' In addition, the object may have the following attributes: \describe{
#' \item{contr.method}{One of 'crunch', 'brunch' or 'piclm'.}
#' \item{macro}{Either 'RRD' or 'PDI' if the response contrasts are calculated
#' as species richness contrasts using \code{\link{macrocaic}}}
#' \item{stand.cont}{A logical value showing whether the contrasts in the
#' object have been standardized.} }
#' @keywords class
#' @export
NULL


#' Example dataset for Fusco imbalance calculations
#'
#' This dataset contains the phylogeny of bird families and species richness
#' originally included with the FUSCO imbalance calculation programs.
#'
#'
#' @aliases fuscoBirdData fuscoBirdTree
#' @format The dataset provides a phylogeny (\code{fuscoBirdTree}) and a data
#' frame (\code{fuscoBirdData}). The phylogeny is a 137 tip tree containing one
#' polytomy and the dataframe provides tip labels and species richness for each
#' bird family. This dataset was provided with the original DOS implementation
#' of the test and the families were unlabelled in this original file.
#' @seealso fusco.test
#' @keywords datasets
#' @export
NULL

#' Example dataset for the caper package
#'
#' This data set contains four species-level comparative datasets used in Isaac
#' et al (2005)
#'
#'
#' @aliases IsaacEtAl chiroptera.tree carnivora.tree primates.tree
#' marsupialia.tree chiroptera.data carnivora.data primates.data
#' marsupialia.data
#' @format The datafile contains species level phylogenies and accompanying
#' data frames of nine variables for each of four mammalian orders (Primates,
#' Carnivora, Chiroptera and Marsupialia). The data were published in
#' supplementary material for Isaac et al. (2005) as CAIC format files and text
#' data files and have been converted for use in 'caper'. The data files are
#' incomplete, with some variables having little or no data for some orders.
#'
#' The variables (all saved as natural log values) are: \describe{
#' \item{species.rich}{Species richness at the tips - all are set to 1 for use
#' in \code{macrocaic}} \item{body.mass}{The average body mass in grams}
#' \item{age.sexual.maturity}{Age at sexual maturity in months}
#' \item{gestation}{Gestation length in days}
#' \item{interbirth.interval}{Interbirth interval in months}
#' \item{litter.size}{The average number of offspring in a litter}
#' \item{population.density}{Population density} \item{group.size}{Number of
#' individuals in a typical group} \item{mass.dimorphism}{Male mass /female
#' mass} \item{length.dimorphism}{Male length / female length} }
#' @seealso caic, pgls
#' @references Isaac, N., Jones, K., Gittleman, J., and Purvis, A. (2005).
#' Correlates of species richness in mammals: Body size, life history, and
#' ecology. American Naturalist, 165(5):600-607.
#' @keywords datasets
#' @examples
#'
#' data(IsaacEtAl)
#' chiroptera <- comparative.data(chiroptera.tree, chiroptera.data, "binomial", na.omit = FALSE)
#' carnivora <- comparative.data(carnivora.tree, carnivora.data, "binomial", na.omit = FALSE)
#' primates <- comparative.data(primates.tree, primates.data, "binomial", na.omit = FALSE)
#' marsupialia <- comparative.data(marsupialia.tree, marsupialia.data, "binomial", na.omit = FALSE)
#' @export
NULL


#' Example dataset for the CAIC package
#'
#' This is a comparative dataset on Perissodactyla taken from the examples
#' include with the original CAIC program.
#'
#'
#' @aliases perissodactyla perissodactyla.data perissodactyla.tree
#' @format The datafile contains a phylogeny (\code{perissodactyla.tree}) of 18
#' perissodactyl species as a 'phylo' object from the \code{ape} library. The
#' tip names are the binomial names of the species. The file also contains a
#' data frame (\code{perissodactyla.data}) of variables 5 variables for 13 of
#' those species: \describe{ \item{Binomial}{The species binomial name.}
#' \item{log.female.wt}{Log female weight} \item{log.gestation.length}{Log
#' gestation length} \item{log.neonatal.wt}{Log neonatal weight}
#' \item{Territoriality}{A factor indicating whether or not the species
#' displays territorial behaviour.} } The dataset is incomplete - it does not
#' include data for each species in the phylogeny and contains missing values.
#' See the examples for the behaviour of the 'comparative.data' function in
#' handling missing data.
#' @seealso caic, pgls
#' @references Purvis, A. and Rambaut, A. (1995). Comparative Analysis by
#' Independent Contrasts (CAIC) User's Guide.
#' @keywords datasets
#' @examples
#'
#' data(perissodactyla)
#' # default behaviour is to omit incomplete data rows
#' (perisso <- comparative.data(
#'     perissodactyla.tree, perissodactyla.data, Binomial
#' ))
#' # but this can be turned off
#' (perisso <- comparative.data(
#'     perissodactyla.tree, perissodactyla.data, Binomial,
#'     na.omit = FALSE
#' ))
#' na.omit(perisso)
#' @export
NULL


#' Generic model methods for 'pgls' models.
#'
#' These are simple summary methods, accessor functions and summary and print
#' methods for 'pgls' models.
#'
#' Phylogenetically corrected residuals from 'pgls' models [TODO].
#'
#' Note that the r^2 values reported by \code{summary.pgls} have a specific
#' interpretation. \code{pgls} fits the intercept-only model for the data using
#' _exactly_ the same covariance matrix (phylogeny plugged through any branch
#' length transformations) as the fitted model to get a null model. The
#' r-squared and adjusted r-squared that are reported therefore hold the
#' covariance matrix constant, so show percentage of variance explained between
#' a null model and the actual model given that precise model of trait change.
#'
#' The actual ML null model for the data (optimising the BL transformation
#' independently) might be different from this - but then the r squared values
#' confound change in explanatory power from changing the model parameters and
#' from changing the trait model.
#'
#' @aliases coef.pgls residuals.pgls fitted.pgls predict.pgls print.pgls
#' summary.pgls print.summary.pgls nobs.pgls
#' @param object An object of class 'pgls'.
#' @param x An object of class 'pgls'.
#' @param phylo Return phylogenetically corrected residuals or ordinary
#' residuals (see details).
#' @param newdata Alternative data for predicting from 'pgls' models.
#' @param digits Number of digits to show in summary methods.
#' @param ... Further arguments to methods.
#' @return The 'summary' method returns an object of class 'summary.pgls'
#' containing: \item{call}{The original function call creating the model.}
#' \item{df}{A vector of the degrees of freedom used to estimate parameters and
#' the residual degrees of freedom.} \item{sigma}{The square root of the
#' estimated variance of the random error.} \item{residuals}{The
#' phylogenetically corrected residuals.} \item{coefficients}{A table of model
#' coefficient, standard errors and t values.} \item{param}{A vector of branch
#' length parameters used in the model.} \item{mlVals}{A vector showing which
#' branch length parameters have been optimised.} \item{param.CI}{A list of
#' length three containing confidence intervals and p values on parameter
#' bounds for each parameter.} \item{fstatistic}{A vector of the F value,
#' numerator and denominator degrees of freedom for the model.}
#' \item{r.squared}{The r^2 for the model.} \item{adj.r.squared}{The adjusted
#' r^2 for the model.}
#' @author Rob Freckleton, David Orme
#' @seealso \code{\link{pgls}}
#' @keywords utils stats
#' @examples
#'
#' data(shorebird)
#' shorebird <- comparative.data(
#'     shorebird.tree, shorebird.data, Species,
#'     vcv = TRUE, vcv.dim = 3
#' )
#' mod1 <- pgls(log(Egg.Mass) ~ log(M.Mass) * log(F.Mass), shorebird)
#' print(mod1)
#'
#' mod1.sum <- summary(mod1)
#' print(mod1.sum)
#' @export
NULL


#' Example dataset for the caper package
#'
#' This is a comparative dataset on the evolution of shorebird egg size taken
#' from Lislevand and Thomas (2006).
#'
#'
#' @aliases shorebird shorebird.data shorebird.tree
#' @format The datafile contains a phylogeny (\code{shorebird.tree}) of 71
#' shorebird species as a 'phylo' object from the \code{ape} library. The tip
#' names are the binomial names of the species. The file also contains a data
#' frame (\code{shorebird.data}) of 71 complete cases for those species. The
#' data frame contains six variables: \describe{ \item{Species}{The species
#' binomial name.} \item{M.Mass}{The adult male mass body mass in grams.}
#' \item{F.Mass}{The adult female mass body mass in grams.} \item{Egg.Mass}{The
#' fresh egg mass in grams.} \item{Cl.size}{The mean clutch size}
#' \item{Mat.syst}{The mating system, as a three level factor: monogamous (MO),
#' polygynous (PO) or polyandrous (PA).} }
#' @seealso crunch, pgls
#' @references Lislevand, T and Thomas, G. H. (2006) Limited male incubation
#' ability and the evolution of egg size in shorebirds. Biology Letters 2, 206
#' - 208
#' @keywords datasets
#' @export
NULL


#' The syrphidae dataset of Katzourakis et al. 2001
#'
#' A genus level phylogeny of the Syrphidae (sawflies) along with data on the
#' species richness of each genus.
#'
#'
#' @aliases syrphidae syrphidaeTree syrphidaeRich
#' @format A 'phylo' object (\code{syrphidaeTree}) containing a phylogeny of
#' 204 genera of sawflies and a data frame (\code{syrphidaeRich}) of the
#' species richness of each genus.
#' @seealso fusco.test
#' @references Katzourakis, A., Purvis, A., Azmeh, S., Rotheray, G., and
#' Gilbert, F. (2001). Macroevolution of hoverflies (Diptera : Syrphidae): the
#' effect of using higher-level taxa in studies of biodiversity, and correlates
#' of species richness. Journal of Evolutionary Biology, 14:219-227.
#' @keywords datasets
#' @export
NULL




# Declare global variables to avoid warnings about usage of variables in
# non-standard evaluation like with() etc.

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "explanatory", "response", "factor", "values", "edge",
        "Nnode", "id", "nTip", "nLin", "tip",
        "N1", "N2", "B", "m", "M", "S.odd", "w", # fusco
        "len", "nSp", "ED", "ED.cor", # pd.calc
        "edge.length", # contrCalc
        "Territoriality"
    ))
}
