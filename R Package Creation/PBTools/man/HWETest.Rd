\name{HWETest}
\alias{HWETest}
\alias{HWETest.default}
\title{
HWETest
}
\description{
This function performs Test on Hardy-Weinberg Equilibrium. This uses the HWE.test.genind function of adegenet package.
}
\usage{
HWETest(genindObject, display=c("full", "matrix"))
}
\arguments{
  \item{genindObject}{genind object}
  \item{display}{value is "full" if the list of complete tests will be displayed; value is "matrix" if matrix of p-values will be displayed}
}
\author{
Nellwyn Sales
}
\keyword{HWE}
