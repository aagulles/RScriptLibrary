###########
#  doGenSim
#' Function for computing relationhsip matrix
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param fileName name of input file which can be pedigree or genotypic file 
#' @param fileFormat format of input file, whether comma-separated file ("csv"), or text file using comma ("ctxt"), space ("stxt"), tab ("ttxt"), or semi-colon ("sctxt") as delimiter
#' @param doPedigree logical; whether to compute pedigree-based relationship matrix or not (marker-based)
#' @param relType type of relationship matrix to compute, whether "add" or "dom" for pedigree-based or "sm-smin","realized","realizedAB","sm", or "additive" for marker-based
#' @param outFileName name of file to save the relationship matrix into
#' @param markerFormat format type
#' 
#
#  Output list containing
#   relFileName - gen_RM$relFileName, 
#   relMat - relationship matrix
###########

doGenSim <- function(outputPath, fileName, fileFormat = c("csv", "ctxt", "stxt", "ttxt", "sctxt"), 
                     doPedigree = TRUE, relType = c("dom", "add", "sm-smin","realized","realizedAB","sm","additive"), outFileName, markerFormat = 1) {
#inputs: fileName, fileFormat, relType

  if (doPedigree) {
  
  #   library(synbreed)
  #   library(pedigreemm)
  
    gen_RM <- pedigree_relationship(outputPath, ped_file = fileName, fFormat = fileFormat, type = relType, outputName = outFileName)
  
  } else {
  
  #   library(synbreed)
  #   library(rrBLUP)
    
    #allow different file formats for geno_file:
    formatsep <- function(fformat) {
        switch(fformat, 
               csv = ",",
               ctxt = ",",
               stxt = " ",
               ttxt = "\t",
               sctxt = ";")
      }
      
    gSep = formatsep(fileFormat)
    
    #Z#  data=as.matrix(read.csv(ori_file,sep=",",row.names=1,stringsAsFactors=FALSE,check.names=FALSE))
###     data <- as.matrix(read.table(fileName, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = gSep, check.names=FALSE))
#   ori_data <- as.matrix(read.table(fileName, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = gSep))
###     code<-recode(ori_file = fileName, type = markerFormat)  #*** wait for revised Data Quality/Preprocessing script from FW
##     code <- recode(outputPath, ori_data, type = markerFormat, nmiss = 0.1, maf = 0.05, cor_threshold = 0.90, packageFormat = "synbreed")

#     code <- GSDataCheck(outputPath, ori_data, type = markerFormat, nmiss = 0.1, maf = 0.05, cor_threshold = 0.90)
    code <- GSDataCheck(outputPath, fileName, type = markerFormat, nmiss = 0.1, maf = 0.05, cor_threshold = 0.90)
#     gen_RM <- marker_relationship(outputPath, geno_file = code$outFileName, fFormat = fileFormat, type = relType, outputName = outFileName)
    gen_RM <- marker_relationship(outputPath, geno_file = code, fFormat = fileFormat, type = relType, outputName = outFileName)
  
  }


  #verify output to return
  return(list(relFileName = gen_RM$relFileName, relMat = gen_RM$relMat))
}

