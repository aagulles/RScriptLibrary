\name{Inoculation}
\alias{Inoculation}
\docType{data}
\title{Inoculation}
\description{
An experiment was conducted on the effect of 6 times of inoculation of S. linicola 
on the oil content of Redwing Flaxseed. The experiment was conducted using an RCB 
design with 4 replications.
}
\usage{data(Inoculation)}
\format{
  A data frame with 24 observations on the following 3 variables.
  \describe{
    \item{\code{Rep}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{InoculationTime}}{a factor with levels \code{Early bloom} \code{Full bloom} \code{Full bloom (1/100)} \code{Ripening} \code{Seedling} \code{Uninoculated}}
    \item{\code{OilContent}}{a numeric vector}
  }
}
\source{
Steel, R.G.D., J.H. Torrie, and D.A. Dickey. Principles and Procedures of Statistics: A Biometrical Approach.
}
\examples{
data(Inoculation)
## maybe str(Inoculation) ; plot(Inoculation) ...
}
\keyword{datasets}
