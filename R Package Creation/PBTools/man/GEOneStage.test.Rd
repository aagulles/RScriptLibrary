\name{GEOneStage.test}
\alias{GEOneStage.test}
\alias{GEOneStage.test.default}
\title{GEOneStage.test}
\description{
This function performs multi-environment one stage test
}
\usage{
GEOneStage.test(exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"), data, respvar, geno, row, column = NULL, rep = NULL, env, is.genoRandom = FALSE)
}
\arguments{
  \item{exptl.design}{RCB, AugRCB, AugLS, Alpha, RowCol, LatinAlpha, LatinRowCol}
  \item{data}{a string; name of the dataframe}
  \item{respvar}{a vector of strings; variable names of the response variables}
  \item{geno}{a string; variable name of the treatment/genotype variable}
  \item{row}{a string; variable name of the blocking variable or row variable}
  \item{column}{a string; variable name of the column variable; NULL, if design is RCB, Alpha, LatinAlpha}
  \item{rep}{a string; variable name of the replication variable; NULL, if design is RCB}
  \item{env}{a string; variable name of the environment variable}
  \item{is.genoRandom}{logical; indicating whether genotype/treatment is random or not; default value is FALSE (FIXED factor)}
}
\author{
Violeta Bartolome, Alaine A. Gulles, Rose Imee Zhella A. Morantte, Nellwyn L. Sales
}
\keyword{multi-site}
