% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Design-methods.R
\name{simulate,DualResponsesDesign-method}
\alias{simulate,DualResponsesDesign-method}
\title{This is a methods to simulate dose escalation procedure using both DLE and efficacy responses.
This is a method based on the \code{\linkS4class{DualResponsesDesign}} where DLEmodel used are of
\code{\linkS4class{ModelTox}} class object and efficacy model used are of \code{\linkS4class{ModelEff}}
class object. In addition, no DLE and efficacy samples are involved or generated in the simulation 
process}
\usage{
\S4method{simulate}{DualResponsesDesign}(
  object,
  nsim = 1L,
  seed = NULL,
  trueDLE,
  trueEff,
  trueNu,
  args = NULL,
  firstSeparate = FALSE,
  parallel = FALSE,
  nCores = min(parallel::detectCores(), 5),
  ...
)
}
\arguments{
\item{object}{the \code{\linkS4class{DualResponsesDesign}} object we want to simulate the data from}

\item{nsim}{the number of simulations (default :1)}

\item{seed}{see \code{\link{setSeed}}}

\item{trueDLE}{a function which takes as input a dose (vector) and returns the true probability 
(vector) of the occurrence of a DLE. Additional arguments can be supplied in \code{args}.}

\item{trueEff}{a function which takes as input a dose (vector) and returns the expected efficacy
responses (vector). Additional arguments can be supplied in \code{args}.}

\item{trueNu}{the precision, the inverse of the variance of the efficacy responses}

\item{args}{data frame with arguments for the \code{trueDLE} and
\code{trueEff} function. The column names correspond to the argument
names, the rows to the values of the arguments. The rows are appropriately
recycled in the \code{nsim} simulations.}

\item{firstSeparate}{enroll the first patient separately from the rest of
the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
in this patient.}

\item{parallel}{should the simulation runs be parallelized across the
clusters of the computer? (not default)}

\item{nCores}{how many cores should be used for parallel computing?
Defaults to the number of cores on the machine, maximum 5.}

\item{\dots}{not used}
}
\value{
an object of class \code{\linkS4class{PseudoDualSimulations}}
}
\description{
This is a methods to simulate dose escalation procedure using both DLE and efficacy responses.
This is a method based on the \code{\linkS4class{DualResponsesDesign}} where DLEmodel used are of
\code{\linkS4class{ModelTox}} class object and efficacy model used are of \code{\linkS4class{ModelEff}}
class object. In addition, no DLE and efficacy samples are involved or generated in the simulation 
process
}
\examples{
##Simulate dose-escalation procedure based on DLE and efficacy responses where no DLE 
## and efficacy samples are used
## we need a data object with doses >= 1:
data <- DataDual(doseGrid=seq(25,300,25),placebo=FALSE)
##First for the DLE model 
##The DLE model must be of 'ModelTox' (e.g 'LogisticIndepBeta') class 
DLEmodel <- LogisticIndepBeta(binDLE=c(1.05,1.8),
                              DLEweights=c(3,3),
                              DLEdose=c(25,300),
                              data=data)

##The efficacy model of 'ModelEff' (e.g 'Effloglog') class 
Effmodel<-Effloglog(Eff=c(1.223,2.513),Effdose=c(25,300),
                    nu=c(a=1,b=0.025),data=data,c=0)

##The escalation rule using the 'NextBestMaxGain' class
mynextbest<-NextBestMaxGain(DLEDuringTrialtarget=0.35,
                            DLEEndOfTrialtarget=0.3)


##The increments (see Increments class examples) 
## 200\% allowable increase for dose below 300 and 200\% increase for dose above 300
myIncrements<-IncrementsRelative(intervals=c(25,300),
                                 increments=c(2,2))
##cohort size of 3
mySize<-CohortSizeConst(size=3)
##Stop only when 36 subjects are treated
myStopping <- StoppingMinPatients(nPatients=36)
##Now specified the design with all the above information and starting with a dose of 25

##Specified the design(for details please refer to the 'DualResponsesDesign' example)
design <- DualResponsesDesign(nextBest=mynextbest,
                              model=DLEmodel,
                              Effmodel=Effmodel,
                              stopping=myStopping,
                              increments=myIncrements,
                              cohortSize=mySize,
                              data=data,startingDose=25)
##Specify the true DLE and efficacy curves
myTruthDLE<- function(dose)
{ DLEmodel@prob(dose, phi1=-53.66584, phi2=10.50499)
}

myTruthEff<- function(dose)
{Effmodel@ExpEff(dose,theta1=-4.818429,theta2=3.653058)
}

##The true gain curve can also be seen
myTruthGain <- function(dose)
{return((myTruthEff(dose))/(1+(myTruthDLE(dose)/(1-myTruthDLE(dose)))))}


## Then specified the simulations and generate the trial 
##For illustration purpose only 1 simulation is produced (nsim=1). 
options<-McmcOptions(burnin=100,step=2,samples=200)
mySim <-simulate(object=design,
                 args=NULL,
                 trueDLE=myTruthDLE,
                 trueEff=myTruthEff,
                 trueNu=1/0.025,
                 nsim=1,
                 seed=819,
                 parallel=FALSE)
}
\keyword{methods}
