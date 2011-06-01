\name{shorebird}
\alias{shorebird}
\alias{shorebird.data}
\alias{shorebird.tree}

\title{Example dataset for the caper package}
\description{
This is a comparative dataset on the evolution of shorebird egg size taken from Lislevand and Thomas (2006).
}
\usage{
data(shorebird)
}
\format{
The datafile contains a phylogeny (\code{shorebird.tree}) of 71 shorebird species as a 'phylo' object from the \code{ape} library. The tip names are the binomial names of the species. The file also contains a data frame (\code{shorebird.data}) of 71 complete cases for those species. The data frame contains six variables: 
\describe{
\item{Species}{The species binomial name.}
\item{M.Mass}{The adult male mass body mass in grams.}
\item{F.Mass}{The adult female mass body mass in grams.}
\item{Egg.Mass}{The fresh egg mass in grams.}
\item{Cl.size}{The mean clutch size}
\item{Mat.syst}{The mating system, as a three level factor: monogamous (MO), polygynous (PO) or polyandrous (PA).}
}
}

\references{Lislevand, T and Thomas, G. H. (2006) Limited male incubation ability and the evolution of egg size in shorebirds. Biology Letters 2, 206 - 208}

\seealso{crunch, pgls}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}

