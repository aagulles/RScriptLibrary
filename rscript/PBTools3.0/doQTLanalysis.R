#########################################################################################
#  doQTLanalysis
#' Function for QTL analysis
#
#  Parameters:
#' @param outputPath folder where graph(s) will be saved
#' @param crossData cross object
#' @param traitType type of traits to be analyzed, whether "Continuous", "Binary", or "Ordinal"
#' @param yVars traits to be analyzed
#' @param mMethod mapping method to be used
#' @param stepCalc step size or the maximum distance in cM between positions at which genotype probabilities are calculated
#' @param errCalc genotyping error rate used in the calculation of the penetrance Pr(observed|true genotype)
#' @param mapCalc map function used when converting genetic distances into recombination fractions, whether "haldane","kosambi","c-f","morgan"
#' @param lodCutoffM cutoff for error LOD scores
#' @param phenoModel phenotype model, whether "normal","binary","2part","np"
#' @param alMethod indicates whether to use "em","imp","hk","ehk"
#' @param nPermutations number of permutation replicates
#' @param numCovar number of marker covariates to use 
#' @param winSize window size in cM
#' @param genoName name of genotype variable
#' @param thresholdWUR theshold for determining LOD cutoff, whether numeric or based on "Li&Ji"
#' @param minDist minimum distance
#' @param stepSize step size used in interval mapping
#' @param addModel logical; whether "additive" (not "dominance" model) is used
#' @param numCofac number of cofactors to set
#' @param mlAlgo logical; whether to use maximum likelihood instead of the default REML as statistical method to use to estimate variance components
#' @param setupModel logical; whether to set-up QTL model or not
#' @param includeEpistasis logical; whether epistasis is included in the model or not
#' @param useDepPrior logical; whether to use dependent prior for indicating variables of epistatic effects or not
#' @param priorMain prior expected number of main effect QTLs 
#' @param priorAll prior expected number for all QTLs on all chromosomes including QTLs with main effects, epistatic 
#'               effects and gene-environment interactions; default is priorMain + 3
#' @param maxQTLs maximum number of QTLs allowed in the model
#' @param priorProb prior inclusion probabilities for epistatic effects: 0.5, 0.1 or 0.05 when both (default), one or none
#'               of the main effects of the two interacting QTL are included in the model           
#
# Packages required: qtl, qtlbim, lattice
#########################################################################################

doQTLanalysis <- function(outputPath, crossData, traitType = c("Continuous", "Binary", "Ordinal"), yVars, mMethod = c("IM", "CIM", "MQM", "BIM"), 
                          stepCalc = 0, errCalc = 0.01, mapCalc = c("haldane","kosambi","c-f","morgan"), lodCutoffM = 3,
                          phenoModel = c("normal","binary","2part","np"), alMethod = c("em","imp","hk","ehk"), nPermutations = 100, 
                          numCovar = 1, winSize = 10, genoName, thresholdWUR =  "Li&Ji", minDist = 10, stepSize = 5.0,
                          addModel = TRUE, numCofac = 1, mlAlgo = TRUE, 
                          setupModel = TRUE, includeEpistasis = FALSE, useDepPrior = FALSE,
                          priorMain = 3, priorAll = priorMain + 3, maxQTLs = NULL, priorProb = c(0.5, 0.1, 0.05)) UseMethod("doQTLanalysis") 
                          


doQTLanalysis.default <- function(outputPath, crossData, traitType = c("Continuous", "Binary", "Ordinal"), yVars, mMethod = c("IM", "CIM", "MQM", "BIM"), 
                          stepCalc = 0, errCalc = 0.01, mapCalc = c("haldane","kosambi","c-f","morgan"), lodCutoffM = 3,
                          phenoModel = c("normal","binary","2part","np"), alMethod = c("em","imp","hk","ehk"), nPermutations = 100, 
                          numCovar = 1, winSize = 10, genoName, thresholdWUR =  "Li&Ji", minDist = 10, stepSize = 5.0,
                          addModel = TRUE, numCofac = 1, mlAlgo = TRUE, 
                          setupModel = TRUE, includeEpistasis = FALSE, useDepPrior = FALSE,
                          priorMain = 3, priorAll = priorMain + 3, maxQTLs = NULL, priorProb = c(0.5, 0.1, 0.05)) { 
                          
  
  cat("QTL ANALYSIS\n\n")
  
  crossData2 <- crossData

  if (mMethod == "IM") {
    
    doQTL <- doQTL_IM(outputPath, crossData2, yVars, stepCalc, errCalc, mapCalc,
         lodCutoffM, phenoModel, alMethod, nPermutations)
      
  } else if (mMethod == "CIM") {
    
    doQTL <- doQTL_CIM(outputPath, crossData2, yVars, stepCalc, errCalc, mapCalc, lodCutoffM, alMethod, 
                       nPermutations, numCovar, winSize)

  } else if (mMethod == "CIM2") {
    
    cat("Method: Composite Interval Mapping\n")
    
    for (i in 1:length(yVars)) { 
      
      cat("\nTrait: ", yVars[i], "\n\n")
      
      QTLselected <- doQTL_CIM2(crossobj = crossData2, QTL.path = outputPath, geno = genoName, env.label = NULL, env = NULL, trait = yVars[i], 
            step = stepCalc, method = "SIM", threshold = thresholdWUR, distance = minDist, cofactors = NULL, window.size = winSize)
    
      if (!is.null(QTLselected$selected$marker)) {
        QTLOutCIMWUR <- doQTL_CIM2(crossobj = crossData2, QTL.path = outputPath, geno = genoName, env.label = NULL, env = NULL, trait = yVars[i], 
           step = stepCalc, method = "CIM", threshold = thresholdWUR, distance = minDist, cofactors = QTLselected$selected$marker, window.size = winSize)  
      }
      
      if (!is.null(QTLOutCIMWUR$selected)) {
        cat("Markers above the threshold:\n\n") 
        printDataFrame(QTLOutCIMWUR$selected)
      }
      
    }
    
  } else if (mMethod == "MQM") {
    
    doQTL <- doQTL_MQM(outputPath, crossData2, yVars, lodCutoffM, nPermutations, winSize, stepSize,
                       addModel, numCofac, mlAlgo)
    
  } else if (mMethod == "BIM") {

    doQTL <- doQTL_BIM(outputPath, crossData2, traitType, yVars, lodCutoffM, nPermutations, setupModel, includeEpistasis, 
		       useDepPrior, priorMain, priorAll, maxQTLs, priorProb)
      
  }
  
}