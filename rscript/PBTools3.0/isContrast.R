# ------------------------------------------------------------------------------
# File and Script Created by: Alaine A. Gulles 08.05.2014
#                             for International Rice Research Institute
# ------------------------------------------------------------------------------
# Description: determine if the linear function is a contrast or not
# Arguments: contrast - a contrast matrix
# Returned Value: logical
# ------------------------------------------------------------------------------

isContrast <- function(contrast) UseMethod("isContrast")

isContrast.default <- function(contrast) {
     result <- FALSE
     theSum <- rowSums(contrast)
     if (all(theSum == 0)) { result <- TRUE }
     return(result)
}

