\name{perissodactyla}
\alias{perissodactyla}
\alias{perissodactyla.data}
\alias{perissodactyla.tree}

\title{Example dataset for the CAIC package}
\description{
This is a comparative dataset on Perissodactyla taken from the examples include with the original CAIC program.
}
\usage{
data(shorebird)
}
\format{
The datafile contains a phylogeny (\code{perissodactyla.tree}) of 18 perrisodactyl species as a 'phylo' object from the \code{ape} library. The tip names are the binomial names of the species. The file also contains a data frame (\code{perissodactyla.data}) of variables 5 variables for 13 of those species: 
\describe{
\item{Binomial}{The species binomial name.}
\item{log.female.wt} {Log female weight}
\item{log.gestation.length} {Log gestation length}
\item{log.neonatal.wt} {Log neonatal weight}
\item{Territoriality} {A factor indicating whether or not the species displays territorial behaviour.}
}
The dataset is incomplete - it does not include data for each species in the phylogeny and contains missing values. See the examples for the behaviour of the 'comparative.data' function in handling missing data. 
}

\references{}
\examples{
data(perrisodactyla)
# default behaviour is to omit incomplete data rows
(perriso <- comparative.data(perrisodactyla.tree, perrisodactyla.data, Binomial))
# but this can be turned off
(perriso <- comparative.data(perrisodactyla.tree, perrisodactyla.data, Binomial, na.omit=FALSE))
na.omit(perriso)
}
\seealso{ caic, pglm }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}

