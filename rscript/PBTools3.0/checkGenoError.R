##############################################################
#  checkGenoError
#' Function for identifying genotyping errors
#'
#  Parameters:
#' @param crossData R object of class cross
#' @param lodCutOff cutoff for error LOD scores
#' @param errorProb genotyping error rate
#'
#' @return list with the elements crossObj and topOut
#   where:
#   crossObj - an R object of class "cross"
#   topOut - table of genotypes with errorlods above cut-off
#
# Package required: qtl
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