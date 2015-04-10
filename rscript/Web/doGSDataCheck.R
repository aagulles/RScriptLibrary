###########
#  doQualityCheck
#' Function for checking GS data quality
#
#  Parameters:
#' @param outputPath
#' @param genoData
#' @param markerFormat, 
#' @param importRel
#' @param maxCorr
#' @param maxMissingP
#' @param minMAF = 0.05
#' 
############

# @param rel_file
# @param rFormat
# @param rMatType
# @param map_file
# @param mFormat

doGSDataCheck <- function(outputPath, geno_data, markerFormat = c(1, 2, 3, 4, 5), 
                          importRel = FALSE, rel_file = NULL, rFormat = NULL, 
                          rMatType = c("t1", "t2", "t3", "t4"), # ped_file = NULL, #peFormat = NULL, #data quality check options, ...,
                          maxCorr = 0.90, maxMissingP = 0.10, minMAF = 0.05, packageFormat                      
) {
  
# doGSDataCheck <- function(outputPath, pheno_file, pFormat, geno_file = NULL, markerFormat = c(1, 2, 3, 4, 5), 
#                       gFormat = NULL, importRel = FALSE, rel_file = NULL, rFormat = NULL, 
#                       rMatType = c("t1", "t2", "t3", "t4"), map_file = NULL, mFormat = NULL,  # ped_file = NULL, #peFormat = NULL, #data quality check options, ...,
#                       maxCorr = 0.90, maxMissingP = 0.10, minMAF = 0.05                      
#   ) {
 
#   #input files into data preparation/quality check function
#   gsData <- GSDataPrep(outputPath, pheno_file, geno_file,  map_file, rel_file, 
#                        # ped_file = NULL, 
#                        pFormat, gFormat, mFormat, rFormat
#                        # , peFormat = NULL
#   ) 
  
  # TODO: perform Data Quality Check using recode_synbreed
#   geno_rec <- recode(ori_data = gsData$genoData, type = markerFormat, nmiss = maxMissingP, maf = minMAF, cor_threshold = maxCorr)
  geno_rec <- recode(outputPath, ori_data = geno_data, type = markerFormat, nmiss = maxMissingP, maf = minMAF, cor_threshold = maxCorr, packageFormat = packageFormat)
  geno_data <- read.table(geno_rec, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",")    

#   gsData$genoData = genoData

#   #create rel file if not available
#   if (!importRel) {
#     setRelMatType <- function(rType) {
#       switch(rType, 
#            t1 = "sm-smin", 
#            t2 = "realized",
#            t3 = "realizedAB",
#            t4 = "sm")
#   }
#   
#   relMatType = setRelMatType(rMatType)
# #   genSimOut <- doGenSim(outputPath, fileName = geno_file, fileFormat = gFormat, doPedigree = FALSE, relType = relMatType, outFileName = rel_file, markerFormat = markerFormat)
#   rFormat = "ctxt"
#   rel_file = paste(outputPath, rel_file, ".m", sep = "") ##need?
#   
#   gen_RM <- marker_relationship(outputPath, geno_file = geno_rec$synbreed, fFormat = "csv", type = relMatType, outputName = rel_file)
# 
# }
# 
# #   gsObject <- create.gpData(pheno = phenoData, geno = genoData, map = mapData) #, pedData)
#   gsObject <- create.gpData(pheno = gsData$phenoData, geno = gsData$genoData, map = gsData$mapData) #, pedData)
#   gsObject$relData <- relData #, pedData)
  
#   return(gsObject = gsObject)
  return(geno_data)
}