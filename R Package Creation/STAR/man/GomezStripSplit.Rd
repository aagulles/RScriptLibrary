\name{GomezStripSplit}
\alias{GomezStripSplit}
\docType{data}
\title{GomezStripSplit}
\description{Data for Strip Split Plot Design}
\usage{data(GomezStripSplit)}
\format{
  A data frame with 108 observations on the following 5 variables.
  \describe{
    \item{\code{REP}}{a factor with levels \code{I} \code{II} \code{III}}
    \item{\code{PLANTING}}{a factor with levels \code{Broadcast} \code{Transplanted}}
    \item{\code{NITROGEN}}{a factor with levels \code{0} \code{60} \code{120}}
    \item{\code{VARIETY}}{a factor with levels \code{V1} \code{V2} \code{V3} \code{V4} \code{V5} \code{V6}}
    \item{\code{YIELD}}{a numeric vector}
  }
}
\examples{
data(GomezStripSplit)
## maybe str(GomezStripSplit) ; plot(GomezStripSplit) ...
}
\keyword{datasets}
