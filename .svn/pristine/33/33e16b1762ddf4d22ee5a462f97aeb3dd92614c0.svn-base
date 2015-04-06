####################################################################################
#  checkQTLdata 
#' Function for data quality checking
#'
#  Parameters:
#' @param outputPath folder where graph(s) will be saved
#' @param crossData cross object
#' @param crossType type of cross used, whether "f2", "bc", "risib", "riself", "bcsft"
#' @param bcNum bc generation, if crosstype is "bcsft"
#' @param fNum filial generation, if crosstype is "bcsft"
#' @param doMissing logical; whether to manage missing genotypic data or not
#' @param deleteMiss logical; whether to delete missing data or perform imputation
#' @param cutOff percentage cutoff for deleting rows/columns with missing data
#' @param doDistortionTest logical; whether to perform Segregation Distortion Test or not
#' @param pvalCutOff borderline probability value for Test of Segregation Distortion, below which the markers will be deleted
#' @param doCompareGeno logical; whether to compare genotypes or not
#' @param cutoffP cutoff probability for comparing genotypes #?
#' @param doCheckMarkerOrder logical; whether to check marker orders or not
#' @param lodThreshold minimum increase in maximum 2-point LOD that will be flagged
#' @param doCheckGenoErrors logical; whether to check for genotyping errors or not
#' @param lodCutOff cutoff for error LOD scores
#' @param errorProb genotyping error rate  
#' 
#' @return list with the elements crossObj, cross, bcGen, fGen
#  where:
#  crossObj - an R object of class "cross"
#  cross - type of cross
#  bcGen - bc generation
#  fGen - filial generation
#
#####################################################################################

checkQTLdata <- function(outputPath, crossData, crossType, bcNum, fNum,
                         doMissing = FALSE, deleteMiss = FALSE, cutOff = NULL,
                         doDistortionTest = FALSE, pvalCutOff, doCompareGeno, cutoffP, 
                         doCheckMarkerOrder, lodThreshold, doCheckGenoErrors, lodCutOff, errorProb
                         ) UseMethod("checkQTLdata")

checkQTLdata.default <- function(outputPath, crossData, crossType, bcNum, fNum,
                         doMissing = FALSE, deleteMiss = FALSE, cutOff = NULL,
                         doDistortionTest = FALSE, pvalCutOff, doCompareGeno, cutoffP, 
                         doCheckMarkerOrder, lodThreshold, doCheckGenoErrors, lodCutOff, errorProb
                         ) {

  
  cat("QTL ANALYSIS\n\nDATA QUALITY CHECK\n\n")
  
  # call function to manage missing data
  if (doMissing) {
    dataDoManMiss <- manageQTLmissing(crossData, deleteMiss, cutOff)
    # - result: genotype data with imputed values/fewer missing values
    crossData <- dataDoManMiss$crossObj
    #print info on missing data
    if (deleteMiss) {
      cat("CHECKING MISSING DATA\n\n")
      if (!is.null(dataDoManMiss$missingRow)) {
        missRows <- data.frame(dataDoManMiss$missingRow)
        colnames(missRows) <- "Row"
        cat("The following rows have too many missing genotypes:\n\n")
        printDataFrame(missRows)
        cat("\n")
      }
      if (!is.null(dataDoManMiss$missingCol)) {
        cat("The following markers have too many missing genotypes:\n\n")
        missCols <- data.frame(dataDoManMiss$missingCol)
        colnames(missCols) <- "Markers"
        printDataFrame(missCols)
        cat("\n")
      }
    }
    cat("\n\n\n")
  }
  
  # call function to perform segregation distortion test
  if (doDistortionTest) {
    dataDoSegTest <- testQTLsegregation(crossData, pvalCutOff)
    crossData <- dataDoSegTest$crossObj
    cat("CHECKING SEGREGATION DISTORTION\n\n")
    if (dim(dataDoSegTest$distTable)[1] > 0) {
      cat("Genotype frequencies at each marker:\n\n")
      printDataFrame(as.data.frame(dataDoSegTest$distTable))
    } else cat(" No markers exhibiting segregation distortion at", pvalCutOff, "level of significance.")
    cat("\n\n\n")
  }
  
  # call function for comparing genotypes
  if (doCompareGeno) {
    dataCompareGeno <- compareGenotypes(crossData, outputPath, cutoffP)
    crossData <- dataCompareGeno$crossObj
    cat("CHECKING FOR SIMILAR GENOTYPES\n\n")
    if (dim(dataCompareGeno$simGenoOut)[1] > 0 && dim(dataCompareGeno$simGenoOut)[2] > 0) {
      cat("Pairs of individuals with similar genotypes:\n\n")
      printDataFrame(as.data.frame(dataCompareGeno$simGenoOut))
    } else cat(" No individuals with more than ", cutoffP*100, "% similar genotype data.", sep ="")
    cat("\n\n\n")
  }
  
  # call function for checking marker order
  if (doCheckMarkerOrder) {
    cat("CHECKING MARKER ORDER\n\n ")
    checkMO <- checkMarkerOrder(crossData, lodThreshold)   
    if (!is.null(checkMO$chkOut)) {
      cat("\nAlleles potentially switched at markers:\n\n")
      printDataFrame(as.data.frame(checkMO$chkOut))
    } else # cat("No flagged markers.")
    crossData <- checkMO$crossObj
    cat("\n\n\n")
  }
  
  # call function for identifying genotyping errors
  if (doCheckGenoErrors) {
    cat("CHECKING GENOTYPING ERRORS\n\n")
    checkGE <- checkGenoError(crossData, lodCutOff, errorProb)
    if (!is.null(checkGE$topOut)) {
      cat("Genotypes with error LOD scores >", lodCutOff)
      printDataFrame(as.data.frame(checkGE$topOut))
    } # else cat("No errorlods above cutoff.")
    crossData <- checkGE$crossObj
    cat("\n\n")
  }
    
  return(invisible(list(crossObj = crossData, cross = crossType, bcGen = bcNum, fGen = fNum)))
}