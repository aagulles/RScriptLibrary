# ------------------------------------------------------------------------------------------
# Description: These functions return all factors of an integer
# Source: http://stackoverflow.com/questions/6424856/r-function-for-returning-all-factors
# ------------------------------------------------------------------------------------------

allFactors <- function(x) UseMethod("allFactors")

allFactors.default <- function(x) {
     x <- as.integer(x)
     div <- seq_len(abs(x))
     factors <- div[x %% div == 0L]
     return(factors)
}