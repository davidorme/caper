\name{growTree}
\alias{growTree}
\alias{as.comparative.data.growTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Tree simulation with traits.  }
\description{
This function provides a very general environment in which to simulate trees. The basic philosophy is that the user provides a series of expressions that define speciation rates, extinction rates and trait evolution. These expressions can make use of information about the internal state of the tree, allowing for very flexible definitions of rules for tree growth.}

\usage{
growTree(b = 1, d = 0, halt = 20, grain = 0.1, linObj = NULL, ct.start = NULL,
         ct.change = NULL, ct.var = NULL, dt.rates = NULL, inheritance = NULL, 
         trace.events = FALSE, trace.cladesize = FALSE, output.lineages = FALSE,
         neg.rates = "abort", inf.rates = "abort", stall.time = 10, 
         extend.proportion=0)
\method{as.comparative.data}{growTree}(x, ...)
}

\arguments{
  \item{b}{A speciation rate. This can be a numeric constant, as in the default, which specifies a single speciation rate for the simulation. Alternatively, this can be an expression, or a list of expressions which define speciation rate in terms of the properties of the tree. See details for discussion of those properties.}
  \item{d}{An extinction rate, described as above.}
  \item{halt}{A rule use to halt the simulation. The default is the number of tips in the simulation, specified as a single integer, but this can also be an expression or list of expressions on the properties of the tree. The simulation is halted when any of these expressions becomes true.}
  \item{grain}{Where rates depend on time or trait values, it becomes necessary to allow time to pass discretely in order to re-evaluate waiting times under the changing values. This sets the amount of time that is allowed to pass before re-evaluation. If rates do not depend on such changing parameters, it is sensible to set this to infinity - this will ensure that the flow of the simulation is not slowed by checking.}
  \item{linObj}{This can be used to supply an existing simulation object, which will then continue to grow under the provided rules. This allows the user to simulate trees with different sets of rules operating in different epochs. The function \code{linToApe} will convert such an object to a 'phylo' object, retaining additional trait data as extra components of the 'phylo' object list. }
  \item{ct.start}{A numeric vector specifying the starting values for continuous traits. If unnamed these will be sequentially named as 'ct1', 'ct2' etc. The names of traits may be used in expressions governing tree growth rules.}
  \item{ct.change}{A numeric vector describing the mean change per unit time in continuous trait values, used to simulate a directional bias in character evolution. If ct.change is NULL, then this is assumed to be zero for each species.}
  \item{ct.var}{Either a vector of variances for each trait or a square matrix describing the variances and covariances amongst the continuous traits. If this is NULL, then uncorrelated traits with a variance of 1 are assumed.}
  \item{dt.rates}{A list of matrices describing the rate of transition between discrete character traits. Each matrix defines a trait and, as with ct.start, the list names are used to identify the traits in the simulation and default to 'dt1', 'dt2', etc. The dimnames of the matrix are used to identify the states of the trait and default to 'st1', 'st2', etc. The matrix need not be symmetrical: the rates are defined from the states in the columns to the states in the rows, hence the diagonal should probably be zero. Each trait is assumed to start the simulation in the first state in the matrix.}
  \item{inheritance}{A list of rules that are applied after a speciation and can be used to modify trait values for the descendent lineages. The names of the list specify which traits are to be modified and, for each trait specified, should return a vector of length two which replaces the existing values. }
  \item{trace.events}{A logical flag, indicating whether or not report speciation, extinction and discrete character evolution events.}
  \item{trace.cladesize}{A positive integer giving an increment size for the simulation to report clade size if required.}
  \item{output.lineages}{A logical flag indicating whether to return the internal lineages object.}
  \item{neg.rates}{One of 'abort', 'warn' or 'quiet', defining the behaviour when a rate calculation produces a negative number. With 'warn'  and 'quiet', negative rates are set to zero and the simulation continues.}
  \item{inf.rates}{One of 'abort', 'warn' or 'quiet', as for\code{neg.rates}. With 'warn' or 'quiet', infinite rates are left in place, resulting in events happening instantly. This may, in some cases, be desirable!}
  \item{stall.time}{If the all rates within the simulation are zero then only this length of time is allowed to pass before the simulation exits with a 'stalled' status. If \code{grain} is infinite, then the simulation stalls immediately when all rates are zero.}
  \item{extend.proportion}{This option allows the simulation to continue running for a given proportion of the time to the next speciation. This makes sense when growing a clade to a given number of extant taxa; with the default setting of zero, the resulting tree ends at a bifurcation with zero branch lengths and this option allows the tree to grow (and taxa to go extinct and traits to evolve).}
  \item{x}{A lineage table output from growTree}
  \item{...}{Further arguments to as.}
}
\details{
The main idea behind this function (which is still in development) is to provide a flexible framework for simulating tree growth and trait evolution. The user provides expressions for the main arguments (\code{b}, \code{d} and \code{halt}) which act as rules defining speciation and extinction and the ending of the simulation. These can be simple constants, but can also make use of the properties of the environment of the evolving tree. This includes both lineage specific properties (as described in the \code{lineages} section of the returned value) or properties of the clade as a whole (as described in the \code{clade} section of the returned value). For example, a extinction rate might increase with lineage age (\code{d=0.01*lin.age}) or a speciation rate might decrease according to a density dependent process (\code{b=1 - (nExtantTip/500)}). Halt expressions will typically use clade properties (\code{halt=clade.age >= 5} or \code{nTips >= 50}) but could use lineage properties, for example stopping when a trait value hits a certain value (\code{halt=any(ct1 >= 10)}). It is not permitted to use '==' in a halt function of clade.age because it will allow the simulation to run away if the actual value steps over the test value.

Discrete traits are defined using a matrices of rates for transitions between states for each trait. At present, these are fixed for the duration of a simulation epoch and cannot be set as expressions of tree variables.

Continuous trait evolution currently employs a simple Brownian model, given a starting value and variance per unit time. The traits can have defined co-variance (the simulation uses mvrnorm at present) and can also have a defined mean change, allowing for a directional walk in the trait values. At present, it is not possible for the trait variance to vary according to the internal state of the tree; continuous characters retain the same variance and covariance for the whole of the simulation epoch.

Whilst none of the \code{halt} rules are TRUE, then the function evaluates the birth, death and discrete trait rates and converts these to waiting times using random variates from a exponential distribution with the calculated rates. These competing waiting times are compared both to each other and the \code{grain} of the simulation, the shortest waiting time is found and the relevant event is then triggered. The winning event is identified in the character vector \code{winnerName} in order to allow inheritance rules to differentiate events.
}
\value{
Depending on the value of \code{output.phylo}, either an object of class 'phylo' or an object of class 'growTree' with the following structure:
  \item{lineages}{A data frame with a row for each lineage in the tree. Each row identifies the \code{parent.id} and \code{id} of the row along with the total age of the lineage (\code{lin.age}) and the times at which the lineage was born (\code{birth.time}). If the species went extinct (or speciated) then the \code{death.time} is recorded and \code{extinct} is set to TRUE. Speciating lineages have  \code{tip} set to FALSE. Each row also records the \code{caic.code} of the lineage - this is used as a sorting code for conversion to a 'phylo' object and is a kludge. If traits are defined in the simulation then the values or states are recorded in this table. These are the current values for extant tips and the values at extinction for extinct tips and internal nodes.}
  \item{clade}{A list containing:\code{clade.age}, the total age of the simulation; \code{nLin}, the total number of lineages; \code{nTip}, the total number of tips, differentated into \code{nExtantTip} and \code{nExtinctTip}.}
  \item{rules}{A list reporting the birth (\code{b}), death (\code{d}) and stopping (\code{halt}) rules and any inheritance rules.}
  \item{ct.set}{If continuous characters were simulated, a list of the \code{ct.start}, \code{ct.change} and \code{ct.var} details provided.}
  \item{dt.rates}{If discrete characters were simulated, a list containing the \code{dt.rates} details provided.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{David Orme, drawing heavily on discussions with Paul-Michael Agapow.}
\examples{

## see the package vignette for a much fuller discussion of examples.

# A basic 200 tip tree, output as a 'comparative.data' object
    tree <- growTree(halt=200, grain=Inf)
    plot(tree$phy)

# A basic tree of age 4 time units, output as a 'comparative.data' object
    tree <- growTree(halt=expression(clade.age >= 4), grain=Inf)
    plot(tree$phy)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
