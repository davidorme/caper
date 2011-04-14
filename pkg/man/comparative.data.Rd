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
comparative.data(data, phy, names.col, vcv.array = FALSE, na.action = "drop", quiet = FALSE)
\method{print}{comparative.data}(x, ...)
\method{na.omit}{comparative.data}(x, ...)
\method{subset}{comparative.data}(x, subset, select, ...)
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
\details{
The function matches rows in a data frame to tips on a phylogeny and ensures correct ordering of the data with respect to the tips. It also can add a variance covariance representation of the phylogeny. This is in part a placeholder for the considerably more sophisticated handling of such data in the package 'phylobase'.}
\value{
A list of class 'comparative.data':

\item{phy      }{An object of class 'phylo'}
\item{data     }{A data frame of matched data}
\item{N        }{The number of tips in the dataset}
\item{data.name}{The original object name of the data}
\item{phy.name }{The original object name of the phylogeny}
\item{na.omit  }{A logical, indicating whether the function has removed incomplete rows}
\item{dropped  }{A list of taxon names dropped from the dataset:
\describe{
\item{incomplete.rows :}{Data rows with incomplete data}
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
phy <- read.tree(text='(((B:2,A:2):1,D:3):1,(C:1,E:1):3);')
dat <- data.frame(taxa=c("A","B","C","D","E"), n.species=c(5,9,12,1,13), mass=c(4.1,4.5,5.9,3.0,6.0))
cdat <- comparative.data(data=dat, phy=phy, names.col="taxa")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

