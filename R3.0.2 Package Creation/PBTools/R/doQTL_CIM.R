######################################################################################
#  doQTL_CIM
#' Function for QTL analysis using Composite Interval Mapping of qtl package
#'
#  Parameters: 
#' @param outputPath folder where graph(s) will be saved
#' @param crossData2 cross object
#' @param yVars traits to be analyzed
#' @param stepCalc step size or the maximum distance in cM between positions at which genotype probabilities are calculated
#' @param errCalc genotyping error rate used in the calculation of the penetrance Pr(observed|true genotype)
#' @param mapCalc map function used when converting genetic distances into recombination fractions, whether "haldane","kosambi","c-f","morgan"
#' @param lodCutoffM cutoff for error LOD scores
#' @param alMethod indicates whether to use "em","imp","hk","ehk"
#' @param nPermutations number of permutation replicates
#' @param numCovar number of marker covariates to use 
#' @param winSize window size in cM
#
# Packages required: qtl
#######################################################################################

doQTL_CIM <- function(outputPath, crossData2, yVars, stepCalc = 0, errCalc = 0.01, mapCalc = c("haldane","kosambi","c-f","morgan"), 
                      lodCutoffM = 3, alMethod = c("em","imp","hk","ehk"), nPermutations = 100, numCovar = 1, winSize = 10)
		UseMethod("doQTL_CIM")
  
doQTL_CIM.default <- function(outputPath, crossData2, yVars, stepCalc = 0, errCalc = 0.01, mapCalc = c("haldane","kosambi","c-f","morgan"), 
                      lodCutoffM = 3, alMethod = c("em","imp","hk","ehk"), nPermutations = 100, numCovar = 1, winSize = 10) {
  
  crossData2 <- calc.genoprob(crossData2, step = stepCalc, error.prob = errCalc, map.function = mapCalc)
  
  cat("Method: Composite Interval Mapping\n")
  
  for (i in 1:length(yVars)) {
    
    #find column location for phenotypes
    colno <- match(yVars[i], names(crossData2$pheno))
    
    cat("\nTrait: ", yVars[i], "\n\n")
    
    if (is.na(lodCutoffM)) {
      
      cimOut <- cim(crossData2, pheno.col = colno, n.marcovar = numCovar, window = winSize,
                    method = alMethod, n.perm = nPermutations)
      pVal = summary(cimOut, alpha = c(0.01, 0.05))
      cat("\n")
      print(pVal)
      lodCutoffM = pVal[2]
    } 
    
    cimOut2 <- cim(crossData2, pheno.col = colno, n.marcovar = numCovar, window = winSize, method = alMethod)
    cimOutDF <- cbind(Marker = rownames(cimOut2),as.data.frame(cimOut2))
    cimOutDFSel <- cimOutDF[which(cimOutDF$lod > lodCutoffM),]
    cat("Markers with lod scores above the cut-off:\n\n") 
    printDataFrame(cimOutDFSel)
    
    #save to file
    graphFilename = paste(outputPath,"/QTLmap_cim_", yVars[i], ".png", sep = "")
    png(graphFilename)
    plot(cimOut2, col = "red")
    #       if (nPermutations == 0) 
    abline(h = lodCutoffM,lty=2, col = "black")
    dev.off()
    
    ##genetic map?  
  }
}
