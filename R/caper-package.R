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
#' @keywords internal
"_PACKAGE"



# Declare global variables to avoid warnings about usage of variables in
# non-standard evaluation like with() etc.

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "explanatory", "response", "factor", "values", "edge",
        "Nnode", "id", "nTip", "nLin", "tip",
        "N1", "N2", "B", "m", "M", "S.odd", "w", # fusco
        "len", "nSp", "ED", "ED.cor", # pd.calc
        "edge.length" # contrCalc
    ))
}

# Need to import some S3 generics from stats to define S3 methods and use
# methods dispatch (e.g. reorder, not explicit ape::reorder.phylo)
#' @importFrom stats na.omit
stats::na.omit
#' @importFrom stats nobs
stats::nobs
#' @importFrom stats reorder
stats::reorder
