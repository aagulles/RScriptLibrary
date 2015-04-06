####################################################################################
#  checkMarkerOrder 
#' Function for checking placement of markers
#
#  Parameters:
#' @param crossData R object of class cross
#' @param lodThreshold minimum increase in maximum 2-point LOD that will be flagged
#'
#' @return a list with the elements crossObj and chkOut
#  where:
#  crossObj - an R object of class "cross"
#  chkOut - a data frame containing the flagged markers
#
# Package required: qtl
####################################################################################

checkMarkerOrder <- function(crossData, lodThreshold) UseMethod("checkMarkerOrder")

checkMarkerOrder.default <- function(crossData, lodThreshold) {
  crossData2 <- crossData
  crossData2 <- est.rf(crossData)
  chkAllOut <- checkAlleles(crossData2, threshold = lodThreshold)

  return(list(crossObj = crossData2, chkOut = chkAllOut))
}