% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data-methods.R
\name{update,Data-method}
\alias{update,Data-method}
\title{Update method for the "Data" class}
\usage{
\S4method{update}{Data}(
  object,
  x,
  y,
  ID = (if (length(object@ID)) max(object@ID) else 0L) + seq_along(y),
  newCohort = TRUE,
  ...
)
}
\arguments{
\item{object}{the old \code{\linkS4class{Data}} object}

\item{x}{the dose level (one level only!)}

\item{y}{the DLT vector (0/1 vector), for all patients in this cohort}

\item{ID}{the patient IDs}

\item{newCohort}{logical: if TRUE (default) the new data are assigned
to a new cohort}

\item{\dots}{not used}
}
\value{
the new \code{\linkS4class{Data}} object
}
\description{
Add new data to the \code{\linkS4class{Data}} object
}
\examples{

# Create some data of class 'Data'
myData <- Data(x=c(0.1,0.5,1.5,3,6,10,10,10),
               y=c(0,0,0,0,0,0,1,0),
               doseGrid=c(0.1,0.5,1.5,3,6,
                          seq(from=10,to=80,by=2)))

## update the data with a new cohort
myData <- update(myData,
                 x=20,
                 y=c(0,1,1))


}
\keyword{methods}
