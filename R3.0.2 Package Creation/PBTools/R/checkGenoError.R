# checkGenoError - function for identifying genotyping errors
#
# REQUIRED input: 
# crossData - R object of class cross
# lodCutOff - cutoff for error LOD scores
# errorProb - genotyping error rate

#
# OUTPUT: a list with the following elements:
#  crossObj - an R object of class "cross"

##########################################################

checkGenoError <- function(crossData, lodCutOff, errorProb) UseMethod("checkGenoError")

checkGenoError.default <- function(crossData, lodCutOff, errorProb) {
  crossData2 <- crossData
  newmap <- est.map(crossData2, errorProb)
  crossData2 <- replace.map(crossData2, newmap)
  crossData2 <- calc.errorlod(crossData2)
  top <- top.errorlod(crossData2, lodCutOff)
  
#   return(topOut = top)
  return(list(crossObj = crossData2, topOut = top))
}