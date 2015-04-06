\name{Insecticides}
\alias{Insecticides}
\docType{data}
\title{Insecticides}
\description{
 Grain yield of rice resulting from use of different 
 foliar and granular insecticides for the control of 
 brown planthoppers and stem borers from a Completely 
 Randomized Design Experiment with Four Replicates and
 7 treatments. 
}
\usage{data(Insecticides)}
\format{
  A data frame with 28 observations on the following 3 variables.
  \describe{
    \item{\code{Insecticide}}{a factor with levels \code{Azodin} \code{Control} \code{DDT + phi-BHC} \code{Dimercon-Boom} \code{Dimercon-Knap} \code{Dol-Mix (1kg)} \code{Dol-Mix (2kg)}}
    \item{\code{Rep}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{GrainYield}}{a numeric vector}
  }
}
\source{
 Gomez, Kwanchai A. and Arturo A. Gomez. (1984) 
 Statistical Procedures for Agricultural Research 2nd 
 Edition. John Wiley & Sons, Inc. Page 14.  
}
\examples{
data(Insecticides)
}
\keyword{datasets}
