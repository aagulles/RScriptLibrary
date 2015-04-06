\name{factorialRegression}
\alias{factorialRegression}
\alias{factorialRegression}
\title{
 factorialRegression
}
\description{
 This function performs factorial regression and saves results in a text file.
}
\usage{
factorialRegression(data, outFileName, respvar, genotype, environment, covariateData, covariateEnvironment, covariateGenotype)
}
\arguments{
  \item{data}{data frame of GxE means}
  \item{outFileName}{path and filename of the text file to be created}
  \item{respvar}{vector of the names of response variable}
  \item{genotype}{string; name of genotype factor}
  \item{environment}{string; name of environment factor}
  \item{covariateData}{data frame of covariate data}
  \item{covariateEnvironment}{string; name of environment factor in the covariate data; NULL if genotype characteristics are used as covariate}
  \item{covariateGenotype}{string; name of genotype factor in the covariate data; NULL if environment characteristics are used as covariate}
}
\author{
Nellwyn Sales
}
\keyword{factorialRegression}
