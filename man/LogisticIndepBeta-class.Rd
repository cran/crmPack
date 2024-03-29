% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model-class.R
\docType{class}
\name{LogisticIndepBeta-class}
\alias{LogisticIndepBeta-class}
\alias{.LogisticIndepBeta}
\title{No initialization function
Standard logistic model with prior in form of pseudo data}
\description{
This is a class for the two-parameter logistic regression DLE model with prior expressed
in form of pseudo data. This model describe the relationhship of the binary DLE (dose-limiting
events) responses and the dose levels. More specifically, this DLE model reprsents the relationship 
of the probabilities of the occurrence of a DLE with their corresponding dose levels in log scale.
This model is specified as 
\deqn{p(d_{(j)})= \frac{exp(\phi_1+\phi_2 log(d_{(j)}))}{1+exp(\phi_1+\phi_2 log(d_{(j)}))}}
for any dose j where \eqn{p(d_{(j)})} is the probability of the occurrence of a DLE at dose j.
The two parameters of this model is the intercept \eqn{\phi_1} and the slope \eqn{\phi_2}
It inherits all slots from \code{\linkS4class{ModelTox}} class.
}
\details{
The pseudo data can be interpreted as as if we obtain some observations before the trial starts.
These pseudo data can be used to express our prior, the initial beliefs for the model parameter(s).
The pseudo data are expressed in the following way. First, fix at least two dose levels which are
Then ask for experts' opinion how many subjects are to be treated at each of these dose levels and
the number of subjects observed with DLE are observed. At each dose level, the number of subjects 
observed with a DLE divided by the total number of subjects treated is the probability of the 
occurrence of a DLE at that particular dose level. The probabilities of the occurrence of a DLE
based on these pseudo data are independent Beta distributions. Therefore, the joint prior probability 
density function of all these probabilities can be obtained. Hence, by a change of variable, the 
joint prior probability density function of the two parameters in this model can also be obtained.
In addition, a conjugate joint perior density function of the two paramaters in the model is used.
For details about the form of all these joint prior and posterior probability density function, please 
refers to Whitehead and Willamson (1998). 


When expressing the pseudo data, \code{binDLE},\code{DLEdose} and \code{DLEweights} are used.
The \code{binDLE} represents the number of subjects observed with DLE. Note that, since the imaginary 
nature of the pseudo data, the number of subjects observed wtih DLE is not necesssary to be integer(s)
but any scalar value.
The \code{DLEdose} represents the dose levels at which the pseudo DLE responses (\code{binDLE}) are 
observed.
The \code{DLEweights} represents the total number of subjects treated.
Since at least two DLE pseudo responses are needed to obtain prior modal estimates (same as the maximum 
likelihood estimates) for the model parameters. \code{binDLE}, \code{DLEdose} and \code{DLEweights} must
all be vectors of at least length 2. Since given one pseudo DLE responses, the number of subjects observed 
with a DLE relates to at which dose level they are treated and the total number of of subjects treated at
this dose level. Therefore, each of the elements in any of the vectors of \code{binDLE}, \code{DLEdose} and 
\code{DLEweights} must have a corresponding elements in the other two vectors. A set of three values with
one of each in the vectors of \code{binDLE}, \code{DLEdose} and \code{DLEweights}. In this model, each of 
these three values must be specified in the same position as in each of the vector of \code{binDLE}, 
\code{DLEdose} and \code{DLEweights}. The order of the values or elements in one of the vector \code{binDLE}, 
\code{DLEdose} and \code{DLEweights} must corresponds to the values or elements specified in the other two 
vectors.
}
\section{Slots}{

\describe{
\item{\code{binDLE}}{represents the vector of pseudo DLE responses. This must be at least f length 2 and the 
order of its elements must corresponds to values specified in \code{DLEdose} and \code{DLEweights}. 
(see details from above)}

\item{\code{DLEdose}}{represents the vector of the corresponding dose levels observed at each of the 
pseudo DLE responses (\code{binDLE}). This mus be at least of length 2 and the order of its elements
must corresponds to values specified in \code{binDLE} and \code{DLEweights}.
(see details from above)}

\item{\code{DLEweights}}{refers to the total number of subjects treated at each of the pseudo dose level 
(\code{DLEdose}). This must be of length of at least 2 and the oreder of its elements must corresponds
to values specified in \code{binDLE} and \code{DLEdose}. (see details from above)}

\item{\code{phi1}}{refers the intercept of the model. This slot is used in output to display the resulting prior 
or posterior modal estimate of the intercept obtained based on the pseudo data and (if any) 
observed data/responses.}

\item{\code{phi2}}{refers to slope of the model. This slot is used in output to display the resulting prior or 
posterior modal estimate of the slope obtained based on the pseudo data and (if any) the observed data/responses.}

\item{\code{Pcov}}{refers to the covariance matrix of the intercept (phi1) and the slope parameters (phi2) of the 
model. This is used in output to display the resulting prior and posterior covariance matrix of phi1 and 
phi2 obtained, based on the pseudo data and (if any) the observed data and responses. This slot is needed for 
internal purposes.}
}}

\examples{
##Obtain prior modal estimates given the pseudo data.
##First we used an empty data set such that only the dose levels under investigations are given.
##In total, 12 dose levels are under investigation ranging from 25 to 300 mg with increments of 25
##(i.e 25, 50, 75, ..., 300).
emptydata<- Data(doseGrid=seq(25,300,25))

##specified our data set is the empty data
data<-emptydata
## Given the pseudo data such that
## Fix two dose level 25 and 300 mg and specified in (DLEdose slot).
## Total number of subjects treated in each of these levels is 3, specified in (DLEweights slot).
## The number of subjects observed with a DLE is 1.05 at dose 25 mg and 1.8 at dose 300 mg, 
## and specified in (binDLE slot).
## the data set we used in the emptydata set, and specified in (data slot).
## Then to modal estimates of the model parameters.
model<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)

##using a data set (see data -class example specification) with observed DLE responses 
##to obtain posterior modal estimates. 
##for the model given the pseudo data

data<-Data(x=c(25,50,50,75,100,100,225,300),y=c(0,0,0,0,1,1,1,1),
           doseGrid=seq(25,300,25))

model<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)

}
\keyword{classes}
