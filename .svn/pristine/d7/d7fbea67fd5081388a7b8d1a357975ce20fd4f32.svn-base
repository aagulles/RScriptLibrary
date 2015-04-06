##################################################################################################
#  doQTL_BIM 
#' Function for QTL analysis using BIM
#'
#  Parameters: 
#' @param outputPath - folder where graph(s) will be saved
#' @param crossData2 - cross object
#' @param traitType - type of traits to be analyzed, whether "Continuous", "Binary", or "Ordinal"
#' @param yVars - traits to be analyzed
#' @param lodCutoffM - cutoff for error LOD scores
#' @param nPermutations - number of permutation replicates
#' @param setupModel - logical; whether to set-up QTL model or not
#' @param includeEpistasis - logical; whether epistasis is included in the model or not
#' @param useDepPrior - logical; whether to use dependent prior for indicating variables of epistatic effects or not
#' @param priorMain - prior expected number of main effect QTLs 
#' @param priorAll - prior expected number for all QTLs on all chromosomes including QTLs with main effects, epistatic 
#'               effects and gene-environment interactions; default is priorMain + 3
#' @param maxQTLs - maximum number of QTLs allowed in the model
#' @param priorProb - prior inclusion probabilities for epistatic effects: 0.5, 0.1 or 0.05 when both (default), one or none
#'               of the main effects of the two interacting QTL are included in the model           
#
# Packages required: qtl, qtlbim
##################################################################################################

doQTL_BIM <- function(outputPath, crossData2, traitType = c("Continuous", "Binary", "Ordinal"), yVars, lodCutoffM = 3,
                      nPermutations = 100, setupModel = TRUE, includeEpistasis = FALSE, useDepPrior = FALSE,
                      priorMain = 3, priorAll = priorMain + 3, maxQTLs = NULL, priorProb = c(0.5, 0.1, 0.05)) UseMethod("doQTL_BIM")


doQTL_BIM.default <- function(outputPath, crossData2, traitType = c("Continuous", "Binary", "Ordinal"), yVars, lodCutoffM = 3,
                      nPermutations = 100, setupModel = TRUE, includeEpistasis = FALSE, useDepPrior = FALSE,
                      priorMain = 3, priorAll = priorMain + 3, maxQTLs = NULL, priorProb = c(0.5, 0.1, 0.05)) {
  
  cat("Method: Bayesian Interval Mapping\n")
  
  for (i in 1:length(yVars)) {
    
    #find column location for phenotypes
    colno <- match(yVars[i], names(crossData2$pheno))
    
    cat("\nTrait: ", yVars[i], "\n\n")
    
    if (traitType == "Continuous") { traitType2 = "normal"
    } else traitType2 = tolower(traitType) 
    
    
    if (setupModel) { 
      capture.output(bimMcmc <- qb.mcmc(crossData, pheno.col = colno, trait = traitType2,
                                        epistasis = includeEpistasis, depend = useDepPrior, main.nqtl = priorMain, mean.nqtl = priorAll, max.nqtl = maxQTLs, prop = priorProb))
      
    } else 
      capture.output(bimMcmc <- qb.mcmc(crossData, pheno.col = colno, trait = traitType2))
    
    sumBimOut <- list()
    bimOutput <- list()
    scanTypes <- c("heritability", "LPD", "LR", "deviance", "detection", "variance", "estimate", "cellmean", "count", "logposterior", "BF", "nqtl")
    
    graphFilename = paste(outputPath,"/QTLmap_bim_", yVars[i], ".png", sep = "")
    png(graphFilename, width = 1920, height = 1440)
    par(mfrow=c(4,3), pty =  "s")
    
    for (j in 1:12) {
      bimOutput[[j]] <- qb.scanone(bimMcmc, epistasis = includeEpistasis, type.scan = scanTypes[j], sum.scan = "no")
      sumBimOut[[j]] <- summary(bimOutput[[j]])
      if (j == 8) { names(sumBimOut[[j]])[5:7] <- c(paste(scanTypes[j],"A",sep="_"), paste(scanTypes[j],"H",sep="_"), paste(scanTypes[j],"B",sep="_"))
      } else { names(sumBimOut[[j]])[5] <- scanTypes[j] }
      #save to file
      #         graphFilename = paste(outputPath,"/QTLmap_bim_", yVars[i], "_", scanTypes[j], ".png", sep = "")
      #         png(graphFilename)
      plot(bimOutput[[j]], col = "red")
      #         dev.off()
    }
    
    dev.off()
    par(mfrow = c(1,1))
    sumBimScan <- sumBimOut[[1]]
    for (k in 2:12 )
      sumBimScan <- merge(sumBimScan, sumBimOut[[k]], by = c("chr", "n.qtl", "pos", "m.pos")) 
    printDataFrame(sumBimScan)
    
    
    #       #save to file
    #       graphFilename = paste(outputPath,"/bim_", yVars[i], ".png", sep = "")
    #       png(graphFilename)
    #       plot(bimOut)
    #       if (nPermutations == 0) 
    #??? cutoff? abline(h = lodCutoffM,lty=2, col = "blue")
    dev.off()
  }
}