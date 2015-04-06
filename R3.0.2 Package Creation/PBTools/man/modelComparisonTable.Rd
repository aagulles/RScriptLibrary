\name{modelComparisonTable}
\alias{modelComparisonTable}
\alias{modelComparisonTable.default}
\title{Model Comparison Table}
\description{
This function constructs model comparison table when model1 and model2 are fitted using REML=TRUE and the term to be tested is random.
}
\usage{
modelComparisonTable(model1, model2)
}
\arguments{
  \item{model1}{resulting object when lmer function is called using REML=TRUE (full model)}
  \item{model2}{resulting object when lmer function is called using REML=TRUE (reduced model)}
}
\author{
Nellwyn L.Sales
}
\keyword{modelComparison}
