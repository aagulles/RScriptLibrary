\name{OutlierDetection}
\alias{OutlierDetection}
\alias{OutlierDetection.default}
\alias{OutlierDetectionMethod1}
\alias{OutlierDetectionMethod2}
\alias{OutlierGraph2}
\alias{OutlierGraph2Sub}
\title{Outlier Detection}
\description{Detecting Outliers}
\usage{
OutlierDetection(data, var, grp = NULL, rep = NULL, path = NULL, method = c("method1", "method2"))
OutlierDetectionMethod1(data, var, grp = NULL)
OutlierDetectionMethod2(data, var, grp)
OutlierGraph2(data, var, grp, rep, path)
OutlierGraph2Sub(data, var, grp, ylabel, rep, min, max, mean)
}
\arguments{
  \item{data}{a data frame}
  \item{var}{a character vector}
  \item{grp}{NULL or a character vector}
  \item{rep}{NULL or character vector}
  \item{path}{a}
  \item{method}{a}
  \item{ylabel}{a}
  \item{min}{a}
  \item{max}{a}
  \item{mean}{a}
}
\author{Violeta I. Bartolome, Alaine Gulles}
\keyword{outlier}
\keyword{analysis}
