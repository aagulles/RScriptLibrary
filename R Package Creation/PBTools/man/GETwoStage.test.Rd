\name{GETwoStage.test}
\alias{GETwoStage.test}
\alias{GETwoStage.test.default}
\title{GETwoStage.test}
\description{
This function performs multi-environment two stage test
}
\usage{
GETwoStage.test(data, respvar, stderr = NULL, sigma2, numrep, geno, env, weight = c("none", "stderr", "stdmean"), is.genoRandom = FALSE)
}
\arguments{
  \item{data}{data frame}
  \item{respvar}{vector of strings; variable names of the reponse variables}
  \item{stderr}{string; column name corresponding to the standard error of the means}
  \item{sigma2}{string; column name corresponding to the residual variances}
  \item{numrep}{string; column name corresponding to the number of replicates}
  \item{geno}{string; name of genotype factor}
  \item{env}{string; name of environment factor}
  \item{weight}{weight options are: "none", "stderr" and "stdmean"}
  \item{is.genoRandom}{logical; TRUE if genotype factor is random}
}
\author{
Violeta Bartolome, Alaine A. Gulles, Rose Imee Zhella Morantte, Nellwyn Sales
}
\keyword{ twostage }
