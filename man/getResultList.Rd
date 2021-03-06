% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Design-methods.R
\name{getResultList}
\alias{getResultList}
\title{Helper function to obtain simulation results list}
\usage{
getResultList(fun, nsim, vars, parallel = NULL)
}
\arguments{
\item{fun}{the simulation function for a single iteration, which takes as
single parameter the iteration index}

\item{nsim}{number of simulations to be conducted}

\item{vars}{names of the variables}

\item{parallel}{shall the iterations be parallelized across the cores?
if NULL, then no parallelization will be done. If scalar positive number, 
then so many cores will be used.}
}
\value{
the list with all simulation results (one iteration corresponds
to one list element)
}
\description{
The function \code{fun} can use variables that are visible to itself. The
names of these variables have to given in the vector \code{vars}.
}
\author{
Daniel Sabanes Bove \email{sabanesd@roche.com}
}
\keyword{internal}
\keyword{programming}
