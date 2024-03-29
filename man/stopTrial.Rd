% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rules-methods.R
\name{stopTrial}
\alias{stopTrial}
\alias{stopTrial,StoppingList,ANY,ANY,ANY,ANY-method}
\alias{stopTrial,StoppingAll,ANY,ANY,ANY,ANY-method}
\alias{stopTrial,StoppingAny,ANY,ANY,ANY,ANY-method}
\alias{stopTrial,StoppingCohortsNearDose,numeric,ANY,ANY,Data-method}
\alias{stopTrial,StoppingPatientsNearDose,numeric,ANY,ANY,Data-method}
\alias{stopTrial,StoppingMinCohorts,ANY,ANY,ANY,Data-method}
\alias{stopTrial,StoppingMinPatients,ANY,ANY,ANY,Data-method}
\alias{stopTrial,StoppingTargetProb,numeric,Samples,Model,ANY-method}
\alias{stopTrial,StoppingMTDdistribution,numeric,Samples,Model,ANY-method}
\alias{stopTrial,StoppingTargetBiomarker,numeric,Samples,DualEndpoint,ANY-method}
\alias{stopTrial,StoppingHighestDose,numeric,ANY,ANY,Data-method}
\alias{stopTrial,StoppingTDCIRatio,ANY,Samples,ModelTox,ANY-method}
\alias{stopTrial,StoppingTDCIRatio,ANY,missing,ModelTox,ANY-method}
\alias{stopTrial,StoppingGstarCIRatio,ANY,Samples,ModelTox,DataDual-method}
\alias{stopTrial,StoppingGstarCIRatio,ANY,missing,ModelTox,DataDual-method}
\title{Stop the trial?}
\usage{
stopTrial(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingList,ANY,ANY,ANY,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingAll,ANY,ANY,ANY,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingAny,ANY,ANY,ANY,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingCohortsNearDose,numeric,ANY,ANY,Data}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingPatientsNearDose,numeric,ANY,ANY,Data}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingMinCohorts,ANY,ANY,ANY,Data}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingMinPatients,ANY,ANY,ANY,Data}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingTargetProb,numeric,Samples,Model,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingMTDdistribution,numeric,Samples,Model,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingTargetBiomarker,numeric,Samples,DualEndpoint,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingHighestDose,numeric,ANY,ANY,Data}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingTDCIRatio,ANY,Samples,ModelTox,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingTDCIRatio,ANY,missing,ModelTox,ANY}(stopping, dose, samples, model, data, ...)

\S4method{stopTrial}{StoppingGstarCIRatio,ANY,Samples,ModelTox,DataDual}(
  stopping,
  dose,
  samples,
  model,
  data,
  TDderive,
  Effmodel,
  Effsamples,
  Gstarderive,
  ...
)

\S4method{stopTrial}{StoppingGstarCIRatio,ANY,missing,ModelTox,DataDual}(stopping, dose, model, data, Effmodel, ...)
}
\arguments{
\item{stopping}{The rule, an object of class
\code{\linkS4class{Stopping}}}

\item{dose}{the recommended next best dose}

\item{samples}{the \code{\linkS4class{Samples}} object}

\item{model}{The model input, an object of class \code{\linkS4class{Model}}}

\item{data}{The data input, an object of class \code{\linkS4class{Data}}}

\item{\dots}{additional arguments}

\item{TDderive}{the function which derives from the input, a vector of the posterior samples called 
\code{TDsamples} of the dose
which has the probability of the occurrence of DLE equals to either the targetDuringTrial or
targetEndOfTrial, the final next best TDtargetDuringTrial (the dose with probability of the 
occurrence of DLE equals to the targetDuringTrial)and TDtargetEndOfTrial estimate.}

\item{Effmodel}{the efficacy model of \code{\linkS4class{ModelEff}} class object}

\item{Effsamples}{the efficacy samples of \code{\linkS4class{Samples}} class object}

\item{Gstarderive}{the function which derives from the input, a vector of the posterior Gstar (the dose
which gives the maximum gain value) samples 
called \code{Gstarsamples}, the final next best Gstar estimate.}
}
\value{
logical value: \code{TRUE} if the trial can be stopped, \code{FALSE}
otherwise. It should have an attribute \code{message} which gives the reason
for the decision.
}
\description{
This function returns whether to stop the trial.
}
\section{Functions}{
\itemize{
\item \code{stopTrial(
  stopping = StoppingList,
  dose = ANY,
  samples = ANY,
  model = ANY,
  data = ANY
)}: Stop based on multiple stopping rules

\item \code{stopTrial(
  stopping = StoppingAll,
  dose = ANY,
  samples = ANY,
  model = ANY,
  data = ANY
)}: Stop based on fulfillment of all multiple stopping
rules

\item \code{stopTrial(
  stopping = StoppingAny,
  dose = ANY,
  samples = ANY,
  model = ANY,
  data = ANY
)}: Stop based on fulfillment of any stopping rule

\item \code{stopTrial(
  stopping = StoppingCohortsNearDose,
  dose = numeric,
  samples = ANY,
  model = ANY,
  data = Data
)}: Stop based on number of cohorts near to next best dose

\item \code{stopTrial(
  stopping = StoppingPatientsNearDose,
  dose = numeric,
  samples = ANY,
  model = ANY,
  data = Data
)}: Stop based on number of patients near to next best
dose

\item \code{stopTrial(
  stopping = StoppingMinCohorts,
  dose = ANY,
  samples = ANY,
  model = ANY,
  data = Data
)}: Stop based on minimum number of cohorts

\item \code{stopTrial(
  stopping = StoppingMinPatients,
  dose = ANY,
  samples = ANY,
  model = ANY,
  data = Data
)}: Stop based on minimum number of patients

\item \code{stopTrial(
  stopping = StoppingTargetProb,
  dose = numeric,
  samples = Samples,
  model = Model,
  data = ANY
)}: Stop based on probability of target tox interval

\item \code{stopTrial(
  stopping = StoppingMTDdistribution,
  dose = numeric,
  samples = Samples,
  model = Model,
  data = ANY
)}: Stop based on MTD distribution

\item \code{stopTrial(
  stopping = StoppingTargetBiomarker,
  dose = numeric,
  samples = Samples,
  model = DualEndpoint,
  data = ANY
)}: Stop based on probability of targeting biomarker

\item \code{stopTrial(
  stopping = StoppingHighestDose,
  dose = numeric,
  samples = ANY,
  model = ANY,
  data = Data
)}: Stop when the highest dose is reached

\item \code{stopTrial(
  stopping = StoppingTDCIRatio,
  dose = ANY,
  samples = Samples,
  model = ModelTox,
  data = ANY
)}: Stop based on 'StoppingTDCIRatio' class when 
reaching the target ratio of the upper to the lower 95% credibility 
interval of the estimate (TDtargetEndOfTrial). This is a stopping rule which incorporate only 
DLE responses and DLE samples are given

\item \code{stopTrial(
  stopping = StoppingTDCIRatio,
  dose = ANY,
  samples = missing,
  model = ModelTox,
  data = ANY
)}: Stop based on 'StoppingTDCIRatio' class
when reaching the target ratio of the upper to the lower 95% credibility 
interval of the estimate (TDtargetEndOfTrial). This is a stopping rule which incorporate only 
DLE responses and no DLE samples are involved

\item \code{stopTrial(
  stopping = StoppingGstarCIRatio,
  dose = ANY,
  samples = Samples,
  model = ModelTox,
  data = DataDual
)}: Stop based on reaching the target ratio of the upper to the lower 95% credibility 
interval of the estimate (the minimum of Gstar and TDtargetEndOfTrial). This is a stopping rule which 
incorporate DLE and efficacy responses and DLE and efficacy samples are also used.

\item \code{stopTrial(
  stopping = StoppingGstarCIRatio,
  dose = ANY,
  samples = missing,
  model = ModelTox,
  data = DataDual
)}: Stop based on reaching the target ratio of the upper to the lower 95% credibility 
interval of the estimate (the minimum of Gstar and TDtargetEndOfTrial). This is a stopping rule which 
incorporate DLE and efficacy responses without DLE and efficacy samples involved.

}}
\examples{

## Example of combining stopping rules with '&' and/or '|' operators

myStopping1 <- StoppingMinCohorts(nCohorts=3)
myStopping2 <- StoppingTargetProb(target=c(0.2, 0.35),
                                  prob=0.5)
myStopping3 <- StoppingMinPatients(nPatients=20)

myStopping <- (myStopping1 & myStopping2) | myStopping3





# Create some data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rules
myStopping1 <- StoppingMinCohorts(nCohorts=3)
myStopping2 <- StoppingTargetProb(target=c(0.2, 0.35),
                                  prob=0.5)
myStopping3 <- StoppingMinPatients(nPatients=20)

# Create a list of stopping rules (of class 'StoppingList') which will then be
# summarized (in this specific example) with the 'any' function, meaning that the study 
# would be stopped if 'any' of the single stopping rules is TRUE.
mystopping <- StoppingList(stopList=c(myStopping1,myStopping2,myStopping3),
                           summary=any)

# Evaluate if to stop the Trial
stopTrial(stopping=myStopping, dose=doseRecommendation$value,
          samples=samples, model=model, data=data)




# Create some data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rules
myStopping1 <- StoppingMinCohorts(nCohorts=3)
myStopping2 <- StoppingTargetProb(target=c(0.2, 0.35),
                                  prob=0.5)
myStopping3 <- StoppingMinPatients(nPatients=20)

# Combine the stopping rules, obtaining (in this specific example) a list of stopping 
# rules of class 'StoppingAll'
myStopping <- (myStopping1 | myStopping2) & myStopping3

# Evaluate if to stop the Trial
stopTrial(stopping=myStopping, dose=doseRecommendation$value,
          samples=samples, model=model, data=data)




# Create some data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rules
myStopping1 <- StoppingMinCohorts(nCohorts=3)
myStopping2 <- StoppingTargetProb(target=c(0.2, 0.35),
                                  prob=0.5)
myStopping3 <- StoppingMinPatients(nPatients=20)

# Combine the stopping rules, obtaining (in this specific example) a list of stopping 
# rules of class 'StoppingAny'
myStopping <- (myStopping1 | myStopping2) | myStopping3

# Evaluate if to stop the Trial
stopTrial(stopping=myStopping, dose=doseRecommendation$value,
          samples=samples, model=model, data=data)




# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if at least 3 
# cohorts were already dosed within 1 +/- 0.2 of the next best dose
myStopping <- StoppingCohortsNearDose(nCohorts = 3,
                                      percentage = 0.2)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          data=data)


# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if at least 9 
# patients were already dosed within 1 +/- 0.2 of the next best dose
myStopping <- StoppingPatientsNearDose(nPatients = 9,
                                       percentage = 0.2)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          data=data)


# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if at least 6 
# cohorts were already dosed 
myStopping <- StoppingMinCohorts(nCohorts = 6)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          data=data)


# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if at least 20 
# patients were already dosed 
myStopping <- StoppingMinPatients(nPatients = 20)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          data=data)


# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if there is at least 
# 0.5 posterior probability that [0.2 =< Prob(DLT | next-best-dose) <= 0.35] 
myStopping <- StoppingTargetProb(target=c(0.2, 0.35),
                                 prob=0.5)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          samples=samples,
          model=model,
          data=data)


# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
             y=c(0, 0, 0, 0, 0, 0, 1, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if there is at least 
# 0.9 probability that MTD > 0.5*next_best_dose. Here MTD is defined as the dose for 
# which prob(DLE)=0.33
myStopping <- StoppingMTDdistribution(target = 0.33,
                                      thresh = 0.5,
                                      prob = 0.9)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          samples=samples,
          model=model,
          data=data)


# Create the data
data <- DataDual(
  x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10,
      20, 20, 20, 40, 40, 40, 50, 50, 50),
  y=c(0, 0, 0, 0, 0, 0, 1, 0,
      0, 1, 1, 0, 0, 1, 0, 1, 1),
  w=c(0.31, 0.42, 0.59, 0.45, 0.6, 0.7, 0.55, 0.6,
      0.52, 0.54, 0.56, 0.43, 0.41, 0.39, 0.34, 0.38, 0.21),
  doseGrid=c(0.1, 0.5, 1.5, 3, 6,
             seq(from=10, to=80, by=2)))

# Initialize the Dual-Endpoint model (in this case RW1)
model <- DualEndpointRW(mu = c(0, 1),
                        Sigma = matrix(c(1, 0, 0, 1), nrow=2),
                        sigma2betaW = 0.01,
                        sigma2W = c(a=0.1, b=0.1),
                        rho = c(a=1, b=1),
                        smooth = "RW1")

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=500)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# In this case target a dose achieving at least 0.9 of maximum biomarker level (efficacy)
# and with a probability below 0.25 that prob(DLT)>0.35 (safety)
myNextBest <- NextBestDualEndpoint(target=c(0.9, 1),
                                   overdose=c(0.35, 1),
                                   maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples,
                               model=model,
                               data=data)

# Define the stopping rule such that the study would be stopped if if there is at 
# least 0.5 posterior probability that the biomarker (efficacy) is within the 
# biomarker target range of [0.9, 1.0] (relative to the maximum for the biomarker). 
myStopping <- StoppingTargetBiomarker(target = c(0.9, 1),
                                      prob = 0.5)

# Evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          samples=samples,
          model=model,
          data=data)

# Create the data
data <- Data(x=c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10, 20, 20, 20, 40, 40, 40,
                 80, 80, 80),
             y=c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             cohort=c(0, 1, 2, 3, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8),
             doseGrid=
               c(0.1, 0.5, 1.5, 3, 6,
                 seq(from=10, to=80, by=2)))

# Initialize the CRM model used to model the data
model <- LogisticLogNormal(mean=c(-0.85, 1),
                           cov=
                             matrix(c(1, -0.5, -0.5, 1),
                                    nrow=2),
                           refDose=56)

# Set-up some MCMC parameters and generate samples from the posterior
options <- McmcOptions(burnin=100,
                       step=2,
                       samples=2000)
set.seed(94)
samples <- mcmc(data, model, options)

# Define the rule for dose increments and calculate the maximum dose allowed
myIncrements <- IncrementsRelative(intervals=c(0, 20),
                                   increments=c(1, 0.33))
nextMaxDose <- maxDose(myIncrements,
                       data=data)

# Define the rule which will be used to select the next best dose
# based on the class 'NextBestNCRM'
myNextBest <- NextBestNCRM(target=c(0.2, 0.35),
                           overdose=c(0.35, 1),
                           maxOverdoseProb=0.25)

# Calculate the next best dose
doseRecommendation <- nextBest(myNextBest,
                               doselimit=nextMaxDose,
                               samples=samples, model=model, data=data)

# Define the stopping rule such that the study would be stopped if there is at least 
# 0.5 posterior probability that [0.2 =< Prob(DLT | next-best-dose) <= 0.35] 
stopTarget <- StoppingTargetProb(target=c(0.2, 0.35),
                                 prob=0.5)

## now use the StoppingHighestDose rule:
stopHigh <-
  StoppingHighestDose() &
  StoppingPatientsNearDose(nPatients=3, percentage=0) &
  StoppingTargetProb(target=c(0, 0.2),
                     prob=0.5)

## and combine everything:
myStopping <- stopTarget | stopHigh

# Then evaluate if to stop the trial
stopTrial(stopping=myStopping, 
          dose=doseRecommendation$value,
          samples=samples,
          model=model,
          data=data)



##define the stopping rules based on the 'StoppingTDCIRatio' class
##Using only DLE responses with samples
## we need a data object with doses >= 1:
data<-Data(x=c(25,50,50,75,150,200,225,300),
           y=c(0,0,0,0,1,1,1,1),
           doseGrid=seq(from=25,to=300,by=25))

##model can be specified of 'Model' or 'ModelTox' class
##For example, the 'logisticIndepBeta' class model
model<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)
##define MCMC options
##for illustration purpose we use 10 burn-in and generate 50 samples
options<-McmcOptions(burnin=10,step=2,samples=50)
##samples of 'Samples' class
samples<-mcmc(data,model,options)
##define the 'StoppingTDCIRatio' class
myStopping <- StoppingTDCIRatio(targetRatio=5,
                                       targetEndOfTrial=0.3)
##Find the next Recommend dose using the nextBest method (plesae refer to nextbest examples)
tdNextBest<-NextBestTDsamples(targetDuringTrial=0.35,targetEndOfTrial=0.3,
                              derive=function(TDsamples){quantile(TDsamples,probs=0.3)})


RecommendDose<-nextBest(tdNextBest,doselimit=max(data@doseGrid),samples=samples,
                        model=model,data=data)
##use 'stopTrial' to determine if the rule has been fulfilled
##use 0.3 as the target proability of DLE at the end of the trial

stopTrial(stopping=myStopping,dose=RecommendDose$nextdose,
          samples=samples,model=model,data=data)
## RecommendDose$nextdose refers to the next dose obtained in RecommendDose
##define the stopping rules based on the 'StoppingTDCIRatio' class
##Using only DLE responses 
## we need a data object with doses >= 1:
data<-Data(x=c(25,50,50,75,150,200,225,300),
           y=c(0,0,0,0,1,1,1,1),
           doseGrid=seq(from=25,to=300,by=25))

##model must be of 'ModelTox' class
##For example, the 'logisticIndepBeta' class model
model<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)
##define the 'StoppingTDCIRatio' class
myStopping <- StoppingTDCIRatio(targetRatio=5,
                                targetEndOfTrial=0.3)
##Find the next Recommend dose using the nextBest method (plesae refer to nextbest examples)
tdNextBest<-NextBestTD(targetDuringTrial=0.35,targetEndOfTrial=0.3)


RecommendDose<-nextBest(tdNextBest,doselimit=max(data@doseGrid),model=model,data=data)
##use 'stopTrial' to determine if the rule has been fulfilled
##use 0.3 as the target proability of DLE at the end of the trial

stopTrial(stopping=myStopping,dose=RecommendDose$nextdose,
          model=model,data=data)
## RecommendDose$nextdose refers to the next dose obtained in RecommendDose
##define the stopping rules based on the 'StoppingGstarCIRatio' class
##Using both DLE and efficacy responses
## we need a data object with doses >= 1:
data <-DataDual(x=c(25,50,25,50,75,300,250,150),
                y=c(0,0,0,0,0,1,1,0),
                w=c(0.31,0.42,0.59,0.45,0.6,0.7,0.6,0.52),
                doseGrid=seq(25,300,25),
                placebo=FALSE)

##DLEmodel must be of 'ModelTox' class
##For example, the 'logisticIndepBeta' class model
DLEmodel<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)

##Effmodel must be  of 'ModelEff' class
##For example, the 'Effloglog' class model
Effmodel<-Effloglog(Eff=c(1.223,2.513),Effdose=c(25,300),nu=c(a=1,b=0.025),data=data,c=0)
##for illustration purpose we use 10 burn-in and generate 50 samples
options<-McmcOptions(burnin=10,step=2,samples=50)
##DLE and efficacy samples must be of 'Samples' class
DLEsamples<-mcmc(data,DLEmodel,options)
Effsamples<-mcmc(data,Effmodel,options)

##define the 'StoppingGstarCIRatio' class
myStopping <- StoppingGstarCIRatio(targetRatio=5,
                                          targetEndOfTrial=0.3)
##Find the next Recommend dose using the nextBest method (plesae refer to nextbest examples)
mynextbest<-NextBestMaxGainSamples(DLEDuringTrialtarget=0.35,
                                   DLEEndOfTrialtarget=0.3,
                                   TDderive=function(TDsamples){
                                     quantile(TDsamples,prob=0.3)},
                                   Gstarderive=function(Gstarsamples){
                                     quantile(Gstarsamples,prob=0.5)})

RecommendDose<-nextBest(mynextbest,doselimit=max(data@doseGrid),samples=DLEsamples,model=DLEmodel,
                        data=data,Effmodel=Effmodel,Effsamples=Effsamples)
##use 'stopTrial' to determine if the rule has been fulfilled
##use 0.3 as the target proability of DLE at the end of the trial 

stopTrial(stopping=myStopping,
          dose=RecommendDose$nextdose,
          samples=DLEsamples,
          model=DLEmodel,
          data=data,
          TDderive=function(TDsamples){
            quantile(TDsamples,prob=0.3)},
          Effmodel=Effmodel,
          Effsamples=Effsamples,
          Gstarderive=function(Gstarsamples){
            quantile(Gstarsamples,prob=0.5)})

## RecommendDose$nextdose refers to the next dose obtained in RecommendDose
##define the stopping rules based on the 'StoppingGstarCIRatio' class
##Using both DLE and efficacy responses
## we need a data object with doses >= 1:
data <-DataDual(x=c(25,50,25,50,75,300,250,150),
                y=c(0,0,0,0,0,1,1,0),
                w=c(0.31,0.42,0.59,0.45,0.6,0.7,0.6,0.52),
                doseGrid=seq(25,300,25),
                placebo=FALSE)

##DLEmodel must be of 'ModelTox' class
##For example, the 'logisticIndepBeta' class model
DLEmodel<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)

##Effmodel must be  of 'ModelEff' class
##For example, the 'Effloglog' class model
Effmodel<-Effloglog(Eff=c(1.223,2.513),Effdose=c(25,300),nu=c(a=1,b=0.025),data=data,c=0)


##define the 'StoppingGstarCIRatio' class
myStopping <- StoppingGstarCIRatio(targetRatio=5,
                                   targetEndOfTrial=0.3)
##Find the next Recommend dose using the nextBest method (plesae refer to nextbest examples)
mynextbest<-NextBestMaxGain(DLEDuringTrialtarget=0.35,DLEEndOfTrialtarget=0.3)

RecommendDose<-nextBest(mynextbest,doselimit=max(data@doseGrid),model=DLEmodel,
                        Effmodel=Effmodel,data=data)

##use 'stopTrial' to determine if the rule has been fulfilled
##use 0.3 as the target proability of DLE at the end of the trial


stopTrial(stopping=myStopping,dose=RecommendDose$nextdose,model=DLEmodel,
          data=data, Effmodel=Effmodel)


## RecommendDose$nextdose refers to the next dose obtained in RecommendDose
}
\keyword{methods}
