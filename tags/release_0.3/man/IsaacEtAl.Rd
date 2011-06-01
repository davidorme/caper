\name{IsaacEtAl}
\alias{IsaacEtAl}
\alias{chiroptera.tree}
\alias{carnivora.tree}
\alias{primates.tree}
\alias{marsupialia.tree}
\alias{chiroptera.data}
\alias{carnivora.data}
\alias{primates.data}
\alias{marsupialia.data}

\title{Example dataset for the caper package}
\description{
This data set contains four species-level comparative datasets used in Isaac et al (2005)
}
\usage{
data(IsaacEtAl)
}
\format{
The datafile contains species level phylogenies and accompanying data frames of nine variables for each of four mammalian orders (Primates, Carnivora, Chiroptera and Marsupialia). The data were published in supplementary material for Isaac et al. (2005) as CAIC format files and text data files and have been converted for use in 'caper'. The data files are incomplete, with some variables having little or no data for some orders.

The variables (all saved as natural log values) are:
\describe{
\item{species.rich}{Species richness at the tips - all are set to 1 for use in \code{macrocaic}}
\item{body.mass}{The average body mass in grams}
\item{age.sexual.maturity}{Age at sexual maturity in months}
\item{gestation}{Gestation length in days}
\item{interbirth.interval}{Interbirth interval in months}
\item{litter.size}{The average number of offspring in a litter}
\item{population.density}{Population density}
\item{group.size}{Number of individuals in a typical group}
\item{mass.dimorphism}{Male mass /female mass}
\item{length.dimorphism}{Male length / female length}
}
}
\examples{
data(IsaacEtAl)
chiroptera  <- comparative.data(chiroptera.tree, chiroptera.data, 'binomial', na.omit=FALSE)
carnivora   <- comparative.data(carnivora.tree, carnivora.data, 'binomial', na.omit=FALSE)
primates    <- comparative.data(primates.tree, primates.data, 'binomial', na.omit=FALSE)
marsupialia <- comparative.data(marsupialia.tree, marsupialia.data, 'binomial', na.omit=FALSE)
	
}
\references{Isaac, N., Jones, K., Gittleman, J., and Purvis, A. (2005). Correlates of species richness in mammals: Body size, life history, and ecology. American Naturalist, 165(5):600-607.}

\seealso{ caic, pgls }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datasets}

