\name{comparative.data}
\alias{comparative.data}
\alias{print.comparative.data}
\title{
Comparative dataset creation
}
\description{
A simple tool to combine phylogenies with datasets and ensure consistent structure and ordering for use in functions. 
}
\usage{
comparative.data(phy, data, names.col, vcv.array = FALSE, na.action = "drop", quiet = FALSE)
\method{print}{comparative.data}(x, ...)
\method{na.omit}{comparative.data}(x, ...)
\method{subset}{comparative.data}(x, subset, select, ...)
\method{reorder}{comparative.data}(x, order, ...)
\method{"["}{comparative.data}(x, i, j)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing variables that can be attributed to the taxa at the tips of a phylogeny.
}
  \item{phy}{
A phylogeny (class 'phylo') to be matched to the data above. 
}
  \item{names.col}{
The name of a column in the provided data frame that will be used to match data rows to phylogeny tips.
}
  \item{vcv}{
A logical value indicating whether to include a variance covariance array representing the phylogeny within the comparative dataset.
}
  \item{vcv.dim}{
Either 2 (a standard VCV matrix) or 3 (an array retaining the individual branches contributing to the standard values). The array form is of use for optimising some branch length transformations.
}
  \item{na.omit}{
A logical value indicating whether to reduce the comparative dataset to those tips for which all selected variables are complete. Note that some functions cannot handle missing data and will return an error.
}
  \item{x}{
An object of class 'comparative.data'.}
}
  \item{subset}{
A logical expression indicating rows of data to keep: missing values are taken as false.
}
  \item{select}{
An expression, indicating columns to select from the data frame.
}
  \item{order}{
One of 'cladewise' or 'pruningwise'. See \code{\link{reorder.phylo}}
}
  \item{i,j}{Indices specifying tips or data columns to extract. See details.}
}


\details{
The function matches rows in a data frame to tips on a phylogeny and ensures correct ordering of the data with respect to the tips. It also can add a variance covariance representation of the phylogeny. Mismatched rows and tips are removed and the taxon labels of these are stored in the 'dropped'  slot of the 'comparative.data' object. The 'print' method displays a brief summary of the dataset contents and the names of the original 'phylo' and 'data.frame' objects. If any rows or tips were dropped, 'print' will also show a venn diagram of the data shared and dropped from each source.

The 'na.omit' and 'subset' methods provide simple ways to clean up and extract parts of the comparative dataset. In particular, 'subset' acts exclusively with the data component of the object and, like subset on a data frame, expects the subset argument to produce a logical vector of data rows to include. The 'reorder' method is use to restructure all the components with the 'comparative.data' object into either pruningwise or cladewise order. This uses code from the 'ape' library: see \code{\link{reorder.phylo}}.

The '[' method allows subsets to be taken of the data. There are no replace methods ('[<-'). If only one index is specified (e.g. x[2]), then this is interpreted as extracting data columns from the object. Otherwise (e.g. x[2,], x[1,1]), the first index will specify tips to extract and the second index will specify columns. Indices for tips are permitted to be numeric, logical or character vectors or empty (missing) or NULL. Numeric values are coerced to integer as by as.integer (and hence truncated towards zero). Character vectors will be matched to the names of the object (or for matrices/arrays, the dimnames): see ‘Character indices’ below for further details.

All of these functions are in part a substitute for the considerably more sophisticated handling of such data in the package 'phylobase', which will be integrated into later releases.
}
\value{
A list of class 'comparative.data':
\item{phy      }{An object of class 'phylo'}
\item{data     }{A data frame of matched data}
\item{data.name}{The original object name of the data}
\item{phy.name }{The original object name of the phylogeny}
\item{dropped  }{A list of taxon names dropped from the dataset:
\describe{
\item{unmatched.rows :}{Data rows that do not match to tips}
\item{tips :}{Tips that do not match to data rows}}}
And optionally:
\item{vcv      }{A variance covariance matrix of the phylogeny}
\item{vcv.dim  }{The dimension of the VCV - 2 for a standard VCV matrix and 3 for an expanded array retaining individual branch lengths}
}
\author{
David Orme
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{crunch}},\code{\link{pglm}}
}
\examples{
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, 'Species')
print(shorebird)

subset(shorebird, subset=Mat.syst == 'MO')

sandpipers <- grep('Calidris', shorebird$phy$tip.label)
shorebird[-sandpipers, ]

sandpipers <- grep('Calidris', shorebird$phy$tip.label, value=TRUE)
shorebird[sandpipers, ]

shorebird[]
shorebird[,]
shorebird[2:3]
shorebird[, 2:3]
shorebird[1:15, ]
shorebird[1:15, 2:3]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

