###########
#  doGBLUP
#' Function for performing GBLUP
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param pheno_file name of phenotypic file
#' @param geno_file name of genotypic file
#' @param markerFormat marker format type
#' @param importRel logical; whether a relationship matrix file is to be imported or not (calculated from genotpic file)
#' @param rel_file name of file containing relationship matrix (if to be imported)
#' @param rMatType type of relationship matrix to be calculated, if importRel is FALSE, whether "t1" = "sm-smin", "t2" = "realized", "t3" = "realizedAB", or "t4" = "sm"
#' @param map_file name of map file
#' @param traitNames names of traits
#' @param covariates name(s) of covariate(s), if any
#' @param doCV logical;
#@param varCompEst = c("BL", "BRR"), 
#' @param samplingStrat sampling strategy to employ, whether "random" or "within popStruc"
#' @param popStruc_file name of file describing the population structure, if samplingStrat is "within popStruc"
#' @param nfolds number of folds, default is 2 
#' @param nrep number of replications
#
###########

doGBLUP <- function(outputPath, pheno_file, geno_file = NULL, markerFormat = c(1, 2, 3, 4, 5), 
                    importRel = FALSE, rel_file = NULL, rMatType = c("t1", "t2", "t3", "t4"), 
                    map_file = NULL, # ped_file = NULL, #peFormat = NULL, #data quality check options, ...,
                    traitNames, covariates = NULL, doCV = FALSE, 
#                     varCompEst = c("BL", "BRR"), 
                    samplingStrat = c("random","within popStruc"), popStruc_file = NULL, nfolds = 2, nrep = 1) {
  
  phenoDataAll <- read.table(pheno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",")    
  
  phenoData <- data.frame(phenoDataAll[, traitNames])
  rownames(phenoData) <- rownames(phenoDataAll)
  colnames(phenoData) <- traitNames
  
  if (!is.null(geno_file)) {
    genoData <- read.table(geno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",")    
  } else geno_file = NULL
  
  if (!is.null(map_file)) {
    mapData <- read.table(map_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",")    
  } else mapData = NULL

  #create rel file if not available
  if (!importRel) {
    setRelMatType <- function(rType) {
      switch(rType, 
             t1 = "realized", 
             t2 = "realizedAB",
             t3 = "sm",
             t4 = "sm-smin")
    }
    
    relMatType = setRelMatType(rMatType)
    rel_file = paste("markerRelation_", relMatType, sep = "") ##need?

    gen_RM <- marker_relationship(outputPath, geno_file, fFormat = "csv", type = relMatType, outputName = rel_file)
  }
  
  if (!is.null(covariates)) { covData = phenoData[, covariates] 
  } else covData = NULL

    gsObject <- create.gpData(pheno = phenoData, geno = genoData, map = mapData, covar = covData) #, pedData)
  
  #receive "output" input files from data preparation/quality check
  #format: files saved in csv or combined to create gpData? ###SECOND OPTION

  #generate genetic value (and selection accuracy- not in new script)
  #gsObject <- gsData$gsObject

  blupOut <- NULL
  blupOutL <- list()
  blupOutCV <- NULL
  blupOutCVL <- list()
  for (i in 1:length(traitNames)) {
    blupOutL[[i]] <- BLUP_synbreed_gv(gsObject, traitNames[i], gen_RM$relMat)
    if (i == 1) { blupOut <- blupOutL[[1]]
    } else blupOut <- merge(blupOut, blupOutL[[i]], by = "row.names")

    #do cv if asked (in separate function?)
    if (doCV) 
    blupOutCVL[[i]] <- BLUP_synbreed_cv(gsObject, traitNames[i], gen_RM$relMat, samplingStrat, popStruc_file,
                       nfolds, nrep)
#     if (i == 1) { blupOutCV <- blupOutCVL[[1]]
#     } else blupOutCV <- merge(blupOutCV, blupOutCVL[[i]], by = "row.names")
  }

    if (length(traitNames) > 1) { write.csv(blupOut[,-1], paste0(outputPath, "blupOut.csv"))
    } else write.csv(blupOut, paste0(outputPath, "blupOut.csv"))

    if (doCV && (!is.null(blupOutCVL))) {
      cat("\nGENOMIC SELECTION\n")
      cat("\nGenome-based Best Linear Unbiased Prediction\n")
      
      cat("\nCross Validation\n")
      
      for (i in 1:length(traitNames)) {
        cat("\nTrait: ", traitNames[i], "\n\n")
        
        print(summary(blupOutCVL[[i]]))
        
        cat("\nSpearman's Rank Correlation\n")
        print(blupOutCVL[[i]]$rankCor)
        
        cat("\nPredictive Ability\n")
        print(blupOutCVL[[i]]$PredAbi)
        
        cat("\nMean squared error\n")
        print(blupOutCVL[[i]]$mse)
        
        cat("\nMean of observed values for the 10% best predicted\n")
        print(blupOutCVL[[i]]$m10)
      }
    }

  #use separate function for generating graphs
  createGSPlots(outputPath, gsObject, gen_RM$relMat, traitNames) #relData) #, rel_file, rFormat, importRel)

}
