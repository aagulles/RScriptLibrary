\name{nc2TestME}
\alias{nc2TestME}
\alias{nc2TestME.default}
\title{nc2TestME}
\description{
 performs multi-environment analysis for North Carolina Design II experiments
}
\usage{
nc2TestME(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, female, male, rep = NULL, block = NULL, row = NULL, column = NULL, inbred = TRUE, individual = NULL, environment)
}
\arguments{
  \item{design}{experiment design}
  \item{data}{name of data frame}
  \item{respvar}{vector of response variables}
  \item{female}{string, name of female factor}
  \item{male}{string, name of male factor}
  \item{rep}{string, name of rep factor}
  \item{block}{string, name of block factor}
  \item{row}{string, name of row factor}
  \item{column}{string, name of column factor}
  \item{inbred}{logical}
  \item{individual}{string, name of individual factor}
  \item{environment}{string, name of environment factor}
}
\author{Nellwyn L. Sales, Violeta Bartolome, Guoyou Ye}
\keyword{nc2me}
