#############################################################################
#  doQTL_IM 
#' Function for doing QTL analysis using Simple Interval Mapping
#
#  Parameters: 
#' @param outputPath folder where graph(s) will be saved
#' @param crossData2 cross object
#' @param yVars traits to be analyzed
#' @param stepCalc step size or the maximum distance in cM between positions at which genotype probabilities are calculated
#' @param errCalc genotyping error rate used in the calculation of the penetrance Pr(observed|true genotype)
#' @param mapCalc map function used when converting genetic distances into recombination fractions, whether "haldane","kosambi","c-f","morgan"
#' @param lodCutoffM cutoff for error LOD scores
#' @param phenoModel phenotype model, whether "normal","binary","2part","np"
#' @param alMethod indicates whether to use "em","imp","hk","ehk"
#' @param nPermutations number of permutation replicates
#
# Packages required: qtl
#############################################################################

doQTL_IM <- function(outputPath, crossData2, yVars, stepCalc = 0, errCalc = 0.01, mapCalc = c("haldane","kosambi","c-f","morgan"),
                 lodCutoffM = 3, phenoModel = c("normal","binary","2part","np"), alMethod = c("em","imp","hk","ehk"),
                 nPermutations = 100) UseMethod("doQTL_IM")

doQTL_IM.default <- function(outputPath, crossData2, yVars, stepCalc = 0, errCalc = 0.01, mapCalc = c("haldane","kosambi","c-f","morgan"),
                 lodCutoffM = 3, phenoModel = c("normal","binary","2part","np"), alMethod = c("em","imp","hk","ehk"),
                 nPermutations = 100) {
  
  crossData2 <- calc.genoprob(crossData2, step = stepCalc, error.prob = errCalc, map.function = mapCalc)
  
  cat("Method: Interval Mapping\n")
  
  for (i in 1:length(yVars)) {
    
    #find column location for phenotypes
    colno <- match(yVars[i], names(crossData2$pheno))
    
    cat("\nTrait: ", yVars[i], "\n\n")
    
    if (is.na(lodCutoffM)) {
      
      imOut <- scanone(crossData2, pheno.col = colno, model = phenoModel, method = alMethod, n.perm = nPermutations)
      pVal = summary(imOut, alpha = c(0.01, 0.05))
      cat("\n")
      print(pVal)
      lodCutoffM = pVal[2]
      
    } 
    
    imOut2 <- scanone(crossData2, pheno.col = colno, model = phenoModel, method = alMethod)
    imOutDF <- cbind(Marker = rownames(imOut2),as.data.frame(imOut2))
    imOutDFSel <- imOutDF[which(imOutDF$lod > lodCutoffM),]
    
    cat("Markers with lod scores above the cut-off:\n\n") 
    printDataFrame(imOutDFSel)
    
    #save to file
    graphFilename = paste(outputPath,"/QTLmap_im_", yVars[i], ".png", sep = "")
    png(graphFilename)
    plot(imOut2, col = "red")
    #       if (nPermutations == 0) 
    abline(h = lodCutoffM, lty=2, col = "black")
    dev.off()
    
    ##genetic map?
    #       imOut <- scanone(crossData, pheno.col=colno, model=c("normal","binary","2part","np"),
    #                        method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
    #                        addcovar=NULL, intcovar=NULL, weights=NULL,
    #                        use=c("all.obs", "complete.obs"), upper=FALSE,
    #                        ties.random=FALSE, start=NULL, maxit=4000,
    #                        tol=1e-4, n.perm, perm.Xsp=FALSE, perm.strata=NULL, verbose,
    #                        batchsize=250, n.cluster=1, ind.noqtl)

  }
  
}
