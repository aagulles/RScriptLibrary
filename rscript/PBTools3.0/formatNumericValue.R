#-------------------------------------------------
# This function formats a number to a specific format
# Author: Nellwyn Sales
#-------------------------------------------------

formatNumericValue <- function(number) UseMethod("formatNumericValue")

formatNumericValue.default <- function(number) {
  if (!is.nan(number)) {
    if (number == 0) {
      newNumber<-formatC(number, digits=6, format="f")
    } else {
      if (number < 0.000001) { newNumber<-formatC(number, digits=2, format="e")
      } else {newNumber<-formatC(number, digits=6, format="f") }
    }
  } else {newNumber<-number }
  return(newNumber)
}