\name{ssaTestPrep}
\alias{ssaTestPrep}
\alias{ssaTestPrep.default}
\alias{ssaTestPrepMain}
\alias{ssaModeling}
\alias{ssaSpatialSED}
\alias{ssaEstimateLSMeans}
\alias{ssaTestPrepResid}
\alias{ssaTestPrepSummary}
\title{Single Site Analysis for p-rep design}
\description{Single Site Analysis for p-rep design}
\usage{
ssaTestPrep(data, respvar, geno, row, column, env = NULL, is.random = FALSE, excludeCheck = FALSE, checkList = NULL, moransTest = FALSE, spatialStruc = c("none", "CompSymm", "Gaus", "Exp", "Spher"), descriptive = FALSE, varCorr = FALSE, heatmap = FALSE, diagplot = FALSE, histogram = FALSE, boxplot = FALSE, outputPath = NULL)

ssaTestPrepMain(data, respvar, geno, row, column, is.random = FALSE, excludeCheck = FALSE, checkList = NULL, moransTest = FALSE, spatialStruc = c("none", "CompSymm", "Gaus", "Exp", "Spher"), varCorr = FALSE)

ssaModeling(data, row, column, formula, spatialStruc)
ssaSpatialSED(model)
ssaEstimateLSMeans(model)

ssaTestPrepResid(object)
ssaTestPrepSummary(object)
}
\arguments{
  \item{data}{a}
  \item{respvar}{a}
  \item{geno}{a}
  \item{row}{a}
  \item{column}{a}
  \item{env}{a}
  \item{is.random}{a}
  \item{excludeCheck}{a}
  \item{checkList}{a}
  \item{moransTest}{a}
  \item{spatialStruc}{a}
  \item{descriptive}{a}
  \item{varCorr}{a}
  \item{heatmap}{a}
  \item{diagplot}{a}
  \item{histogram}{a}
  \item{boxplot}{a}
  \item{outputPath}{a}
  \item{formula}{a}
  \item{model}{a}
  \item{object}{a}
}
\author{Violeta Bartolome, Alaine A. Gulles}
\keyword{analysis}