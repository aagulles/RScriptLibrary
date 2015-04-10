doGBLUP <- function(pheno_File, pFormat, pSep, geno_file, gFormat, gSep, rel_file = NULL, rFormat = NULL, rSep = NULL,
                    traitName, doCrossVal = FALSE, varCompEst = c("BL", "BRR"), samplingStrat = c("random","within popStruc"), nfolds = 2, nrep = 1) {
  
  ###TO DO:
  #pass name of pheno_file XXXor only get selected trait from pheno_file?
  #det. trait index from list of traits
  
  #recode genotypic file
  code<-recode(ori_file = geno_file, type = markerFormat) #*** wait for revised Data Quality/Preprocessing script from FW  
  
  #pass on name of rel_file or calculate from genotypic file and save to a file for passing on
  if (calculate) {
    gen_RM <- marker_relationship(outputPath, geno_file = code$synbreed, fFormat = fileFormat, type = relType, outputName = "relMatFile")
    rel_file <- paste(outputPath, relMathFile, ".m")
  } else 
  BLUP_synbreed(pheno_file, pFormat, traitName, rel_file, doCV = doCrossVal, varCompEst, sampling = samplingStrat, popStruc_file = NULL, k = nfolds, Rep = nrep)
    
  
  
}
