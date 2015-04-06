##########################################################
#  testQTLSegregation 
#' Function for performing test for segregation distortion
#'
#  Parameters:
#' @param crossData R object of class cross
#' @param pvalCutOff borderline probability value for Test of Segregation Distortion, below which the markers will be deleted
#'
#' @return list with the elements crossObj and distTable
#  where:
#  crossObj - an R object of class "cross"
#  distTable - table of markers and corresponding columns deleted from the data))#
##########################################################

testQTLsegregation <- function(crossData, pvalCutOff) UseMethod("testQTLsegregation")

testQTLsegregation.default <- function(crossData, pvalCutOff) {
  
  gt <- geno.table(crossData)
  pvals <- gt$P.value
  distOut <- gt[!is.na(pvals) & pvals < pvalCutOff,] 
  if (dim(distOut)[1] > 0)
    crossDataNoDist <- drop.markers(crossData, markers = rownames(distOut))
  else crossDataNoDist = crossData

  return(list(crossObj = crossDataNoDist, distTable = distOut))
}