\name{designAugmented}
\alias{designAugmented}
\alias{designAugmented.default}
\title{Randomization for Augmented Design}
\description{
  Generate randomization of augmented design in
  randomized complete block design and latin square.
}
\usage{
designAugmented(checkTrmt, newTrmt, r = NULL, trial = 1, design = c("rcbd", "lsd"), factorName = NULL, display = TRUE, file = NULL)
}
\arguments{
  \item{checkTrmt}{levels of replicate treatment or number of levels} 
  \item{newTrmt}{levels of unreplicate treament or the number of levels}
  \item{r}{number of replicate or NULL}
  \item{trial}{number of trial}
  \item{design}{experimental design used}
  \item{factorName}{string or NULL}
  \item{display}{logical}
  \item{file}{path}
}
\author{AAG}
\keyword{ design }
