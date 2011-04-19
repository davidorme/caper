\name{caic-class}
\alias{caic}
\alias{caic.diagnostics}
\alias{caic.table}
\alias{summary.caic}

\title{
The 'caic' S3 object class and methods
}
\description{
The functions 'crunch', 'brunch', 'macrocaic' and 'piclm' all return an object containing an independent contrast model. The structure of the object and the available methods are described here.
}
\usage{caic.diagnostics}
\usage{caic.table}
\usage{summary.caic}

\format{
A 'caic' object is a list containing the following:
\describe{
    \item{contrast.data}{ A list of the following:
    \describe{
        \item{contr}{A list containing matrices of the contrasts in the response variable (contr\$response) and explanatory variables (contr\$explanatory).}
        \item{nodalVals}{A list containing matrices of the nodal values in the response variable (contr\$response) and explanatory variables (contr\$explanatory).}
        \item{contrVar}{A numeric vector of the calculated variances for each contrast.}
        \item{nChild}{A vector showing the number of nodes descending from each internal node}
        \item{nodeDepth}{A vector showing the maximum number of nodes between each internal node and the tips of the phylogeny (including both the node in question and the tip and hence always >=2) }
        \item{validNodes}{A logical vector showing which internal nodes on the tree have valid contrasts, given the available data and any user constraints.}
    }
    \item{data}{A 'comparative.data' object containing the phylogeny used to calculate contrasts and the original data.}
    \item{lm}{An 'lm' object containing a regression model through the origin for the calculated contrast}
}

In addition, the object may have the following attributes:
\describe{

\item{contr.method}{One of 'crunch', 'brunch' or 'piclm'.}
\item{macro}{Either 'RRD' or 'PDI' if the response contrasts are calculated as species richness contrasts using \code{\link{macrocaic}}}
\item{stand.cont}{A logical value showing whether the contrasts in the object have been standardized.}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

}
\keyword{class}
