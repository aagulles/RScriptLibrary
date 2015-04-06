\name{gge.analysis}
\alias{gge.analysis}
\alias{gge.analysis.default}
\title{gge.analysis}
\description{
performs gge analysis
}
\usage{
gge.analysis(ENV, GEN, REP, Y, MSE = 0, number = TRUE, yVar = "Y1", f=0.5, graphSym=FALSE, graphEnv=FALSE, graphGeno=FALSE, ...)
}
\arguments{
  \item{ENV}{vector containing values of environment factor}
  \item{GEN}{vector containing values of genotype factor}
  \item{REP}{vector containing values of rep factor; or a single value of the number of replicates}
  \item{Y}{vector containing the values of the response variable}
  \item{MSE}{vector containing values of residual variance; or a single value of the residual variance}
  \item{number}{TRUE if levels of genotype and environment factors are numeric}
  \item{yVar}{name of response variable}
  \item{f}{0.5 if symmetric view; 0 if environment view; 1 if genotype view}
  \item{graphSym}{TRUE if GGE biplot in symmetric view will be created}
  \item{graphEnv}{TRUE if GGE biplot in environment view will be created}
  \item{graphGeno}{TRUE if GGE biplot in genotype view will be created}
  \item{\dots}{other graph parameters}
}
\author{Felipe de Mendiburu, Nellwyn L. Sales, Violeta Bartolome}
\keyword{gge}
