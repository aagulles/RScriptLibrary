\name{nc3Test}
\alias{nc3Test}
\alias{nc3Test.default}
\title{nc3Test}
\description{
 performs analysis for North Carolina Design III experiments
}
\usage{
nc3Test(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, tester, f2lines, rep = NULL, block = NULL, row = NULL, column = NULL, individual = NULL, environment = NULL)
}
\arguments{
  \item{design}{experiment design}
  \item{data}{name of data frame}
  \item{respvar}{vector of response variables}
  \item{tester}{string, name of tester factor}
  \item{f2lines}{string, name of F2 lines factor}
  \item{rep}{string, name of rep factor}
  \item{block}{string, name of block factor}
  \item{row}{string, name of row factor}
  \item{column}{string, name of column factor}
  \item{individual}{string, name of individual factor}
  \item{environment}{string, name of environment factor}
}
\author{Nellwyn L. Sales, Violeta Bartolome, Guoyou Ye}
\keyword{nc3}
