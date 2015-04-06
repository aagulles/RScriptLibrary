\name{GEOneStage_resid}
\alias{GEOneStage_resid}
\alias{GEOneStage_resid.default}
\title{GEOneStage_resid}
\description{
This function consolidates the residuals generated after running GEOneStage.test
}
\usage{
GEOneStage_resid(geOneResult, respvar, is.genoRandom = FALSE)
}
\arguments{
  \item{geOneResult}{object generated after running GEOneStage.test}
  \item{respvar}{a vector of strings; variable names of the response variables}
  \item{is.genoRandom}{logical; indicating whether genotype/treatment is random or not; default value is FALSE (FIXED factor)}
}
\author{Alaine Gulles, Rose Imee Zhella Morantte, Nellwyn Sales}
\keyword{residuals}
