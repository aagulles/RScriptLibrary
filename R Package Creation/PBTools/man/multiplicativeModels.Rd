\name{multiplicativeModels}
\alias{multiplicativeModels}
\alias{multiplicativeModels.default}
\title{multiplicativeModels}
\description{
This function performs ammi analysis, gge analysis and writes the results to a text file
}
\usage{
multiplicativeModels(data, outFileName, respvar, environment, genotype, numberOfReps, residualVariances, responsePlot=FALSE, doAMMI=TRUE, biplotPC12 = FALSE, biplotPC13 = FALSE, biplotPC23 = FALSE, ammi1 = FALSE, adaptMap = FALSE, doGGE=FALSE, graphSym=FALSE, graphEnv=FALSE, graphGeno=FALSE)
}
\arguments{
  \item{data}{data frame}
  \item{outFileName}{path and filename of the text file to be created}
  \item{respvar}{vector of the names of response variable}
  \item{environment}{string; name of environment factor}
  \item{genotype}{string; name of genotype factor}
  \item{numberOfReps}{vector of the column names corresponding to the number of rep; or vector of numeric values of number of rep}
  \item{residualVariances}{vector of the column names corresponding to the residual variance; or vector of numeric values of residual variance}
  \item{responsePlot}{TRUE if response plots will be created}
  \item{doAMMI}{TRUE if AMMI analysis will be performed}
  \item{biplotPC12}{TRUE if biplot (PC1 vs PC2) will be created}
  \item{biplotPC13}{TRUE if biplot (PC1 vs PC3) will be created}
  \item{biplotPC23}{TRUE if biplot (PC2 vs PC3) will be created}
  \item{ammi1}{TRUE if ammi1 biplot will be created}
  \item{adaptMap}{TRUE if adaptation map will be created}
  \item{doGGE}{TRUE if GGE analysis will be performed}
  \item{graphSym}{TRUE if GGE biplot in symmetric view will be created}
  \item{graphEnv}{TRUE if GGE biplot in environment view will be created}
  \item{graphGeno}{TRUE if GGE biplot in genotype view will be created}
}
\author{Nellwyn L. Sales}
\keyword{multiplicativeModels}
