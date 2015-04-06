\name{gge.analysis2}
\alias{gge.analysis2}
\alias{gge.analysis2.default}
\title{gge.analysis2}
\description{
This function performs gge analysis for multi-env two-stage. It computes several values first before calling gge.analysis function.
}
\usage{
gge.analysis2(ENV, GEN, REP, Y, SIGMA2, number = TRUE, graphSym=FALSE, graphEnv=FALSE, graphGeno=FALSE, respVar, f=0.5)
}
\arguments{
  \item{ENV}{vector containing values of environment factor}
  \item{GEN}{vector containing values of genotype factor}
  \item{REP}{vector containing values of rep factor}
  \item{Y}{vector containing the values of the response variable}
  \item{SIGMA2}{vector containing values of residual variance}
  \item{number}{TRUE if levels of genotype and environment factors are numeric}
  \item{graphSym}{TRUE if GGE biplot in symmetric view will be created}
  \item{graphEnv}{TRUE if GGE biplot in environment view will be created}
  \item{graphGeno}{TRUE if GGE biplot in genotype view will be created}
  \item{respVar}{name of response variable}
  \item{f}{0.5 if symmetric view; 0 if environment view; 1 if genotype view}
}
\author{Alaine A. Gulles, Nellwyn L. Sales}
\keyword{gge2}
