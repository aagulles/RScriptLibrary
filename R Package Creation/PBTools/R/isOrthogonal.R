# ------------------------------------------------------------------------------
# File and Script Created by: Alaine A. Gulles 08.05.2014
#                             for International Rice Research Institute
# ------------------------------------------------------------------------------
# Description: determine if the linear contrast are pairwise orthogonal
# Arguments: contrast - a contrast matrix
# Returned Value: logical
# ------------------------------------------------------------------------------

isOrthogonal <- function(contrast) UseMethod("isOrthogonal")

isOrthogonal.default <- function(contrast) {
     numContrast <- nrow(contrast)
     if (numContrast == 1) {
          stop("The number of contrast should be greater than or equal to 2.")
     }
     result <- FALSE
     theContrastName <- rownames(contrast)
     pairedContrast <- combn(theContrastName, 2)
     for (i in ncol(pairedContrast)) {
          index <- match(pairedContrast[,i], theContrastName)
          if (sum(apply(contrast[index,],2,prod)) != 0) {
               result <- FALSE
               break
          } else { result <- TRUE }
     }
     return(result)
}