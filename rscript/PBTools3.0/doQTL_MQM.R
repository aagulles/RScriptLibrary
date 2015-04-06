#########################################################################
#  doQTL_MQM
#' Function for QTL analysis using MQM
#
#  Parameters:
#' @param outputPath folder where graph(s) will be saved
#' @param crossData2 cross object
#' @param yVars traits to be analyzed
#' @param lodCutoffM cutoff for error LOD scores
#' @param nPermutations number of permutation replicates
#' @param winSize window size in cM
#' @param stepSize step size used in interval mapping
#' @param addModel logical; whether "additive" (not "dominance" model) is used
#' @param numCofac number of cofactors to set
#' @param mlAlgo logical; whether to use maximum likelihood instead of the default REML as statistical method to use to estimate variance components
#
#  Packages required: qtl
#########################################################################

doQTL_MQM <- function(outputPath, crossData2, yVars, lodCutoffM = 3, nPermutations = 100, winSize = 10, stepSize = 5.0,
                      addModel = TRUE, numCofac = 1, mlAlgo = TRUE) 
		UseMethod("doQTL_MQM")
  
doQTL_MQM.default <- function(outputPath, crossData2, yVars, lodCutoffM = 3, nPermutations = 100, winSize = 10, stepSize = 5.0,
                      addModel = TRUE, numCofac = 1, mlAlgo = TRUE) {
  
  cat("Method: Multiple QTL Mapping\n")
  
  for (i in 1:length(yVars)) {
    
    #find column location for phenotypes
    colno <- match(yVars[i], names(crossData2$pheno))
    
    
    cat("\nTrait: ", yVars[i], "\n\n")
    crossDataAug <- mqmaugment(crossData2)
    if (addModel) { modelMQM = "additive"} else modelMQM = "dominance"
    cofacList <- mqmautocofactors(crossDataAug, numCofac)
    
    if (is.na(lodCutoffM)) {
      
      mqmOut <- mqmpermutation(crossDataAug, scanfunction = mqmscan, pheno.col = colno, multicore = FALSE, 
                               n.perm = nPermutations, cofactors = cofacList)
      mqmOutRes  <- mqmprocesspermutation(mqmOut)
      pVal = summary(mqmOutRes, alpha = c(0.01, 0.05))
      cat("\n")
      print(pVal)
      lodCutoffM = pVal[2]
    } 
    mqmOut2 <- mqmscan(crossDataAug, cofactors = cofacList, pheno.col = colno, 
                       model = modelMQM, forceML = mlAlgo, window.size = winSize, step.size = stepSize)
    mqmOutMarkers <- mqmextractmarkers(mqmOut2)
    mqmOutDF <- cbind(Marker = rownames(mqmOutMarkers),as.data.frame(mqmOutMarkers))
    varname = paste("LOD", yVars[i])
    mqmOutDFSel <- mqmOutDF[which(mqmOutDF[,varname] > lodCutoffM),]
    
    cat("Markers with lod scores above the cut-off:\n\n") 
    printDataFrame(mqmOutDFSel)
    
    #save to file
    graphFilename = paste(outputPath,"/QTLmap_mqm_", yVars[i], ".png", sep = "")
    png(graphFilename)
    plot(mqmOut2, col = "red")
    #       if (nPermutations == 0) 
    abline(h = lodCutoffM,lty=2, col = "black")
    dev.off()
    
    #plot genetic map
    graph2Filename = paste(outputPath,"/mqm_geneticmap_", yVars[i], ".png", sep = "")
    png(graph2Filename)
    plot(mqmgetmodel(mqmOut2))
    dev.off()
    
    #       #plot directed qtl
    #       graph3Filename = paste(outputPath,"/mqm_directed_", yVars[i], ".png", sep = "")
    #       png(graph3Filename)
    #       plot(mqmplot.directedqtl(crossData2, mqmOut, pheno.col = colno))
    #       dev.off()
  }   
}