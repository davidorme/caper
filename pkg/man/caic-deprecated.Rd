\name{caper-deprecated}
\alias{max.lik.lambda}
\alias{max.lik.kappa}
\alias{profile.lambda}
\alias{profile.kappa}
\alias{all.clades}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Deprecated functions in the caper package}
\description{These functions are deprecated and are now simply redirections to the replacement functions}
\usage{
\method{profile}{lambda}(x, dat, V)
\method{profile}{kappa}(x, V)
\method{max.lik}{lambda}(x, V, dat)
\method{max.lik}{kappa}(x, V)
\method{all}{clades}(phy, tips=FALSE, tip.labels=FALSE)
%profile.lambda(x, dat, V)
%profile.kappa(x, V)
%max.lik.lambda(x, V, dat)
%max.lik.kappa(x, V)
%all.clades(phy, tips=FALSE, tip.labels=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of data values matched to the tips of the variance covariance matrix.}
  \item{dat}{A data frame, containing a variable to use for likelihood calculation.}
  \item{V}{A variance-covariance matrix representing a phylogeny.}
  \item{phy}{An object of class 'phylo'. }
  \item{tips}{ A logical indicating whether to include external node membership in the list.}
  \item{tip.labels}{A logical flag indicating whether to return the node numbers of the tips or their tip labels.}
}

\details{
Although odd to deprecate a function on it's first incarnation in a package, the profile and max.lik functions come form the pglm code which saw extensive use as a stand alone file of functions. These functions were renamed in order to avoid clashes in the naming with the \code{max} and \code{profile} generic S3 functions and the deprecated versions are included in order to ease the transition of users from the stand alone file to the packaged versions.

The function \code{all.clades} has been deprecated for the same reasons.
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ documentation }
