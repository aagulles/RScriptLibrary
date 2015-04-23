\name{recodeDiallelData}
\alias{recodeDiallelData}
\alias{recodeDiallelData.default}
\title{
recode diallel data
}
\description{
This function recodes the input data to be suitable for diallel analysis.
}
\usage{
recodeDiallelData(design = c("diallel1", "diallel2", "diallel3", "diallel4"), data, p1, p2, rep, block=NULL, row=NULL, column=NULL, respvar)
}
\arguments{
  \item{design}{experiment design}
  \item{data}{name of data frame}
  \item{p1}{string, name of p1 factor}
  \item{p2}{string, name of p2 factor}
  \item{rep}{string, name of rep factor}
  \item{block}{string, name of block factor}
  \item{row}{string, name of row factor}
  \item{column}{string, name of column factor}
  \item{respvar}{vector of response variables}
}
\author{
Nellwyn L. Sales
}
\keyword{recodeDiallel}
