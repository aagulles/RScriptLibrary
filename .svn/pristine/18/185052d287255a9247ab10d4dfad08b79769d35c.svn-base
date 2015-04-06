###########
#  marker_relationship
#' Function for computing marker-based relationship matrix
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param geno_file name of genotypic input file
#' @param fFormat format of genotypic input file, whether comma-separated file ("csv"), or text file using comma ("ctxt"), space ("stxt"), tab ("ttxt"), or semi-colon ("sctxt") as delimiter
#' @param type type of relationship matrix to compute, whether "sm-smin","realized","realizedAB","sm", or "additive"
#' @param outputName name of file to save the relationship matrix into
#' 
#  Output list containing
#   relFileName - gen_RM$relFileName, 
#   relMat - relationship matrix
###########

# compute marker-based relationship matrix

# input
## geno_file: recoded SNP markers (synbreed or rrBLUP format, limited to method type)
## type: method for computing relationship matrix
### "realized","realizedAB", "sm", and "sm-smin": methods provided by package synbreed
### "additive": method provided by package rrBLUP

### genotype file (csv file *** or tab-delimited tab file)
# rows: individual
# rownames: individual ID
# columns: marker
# columnnames: marker ID

# output
## marker-based relationship matrix


# library(synbreed)
# library(rrBLUP)

marker_relationship<-function(outputPath, geno_file, fFormat = c("csv", "ctxt", "stxt", "ttxt", "sctxt"), type=c("sm-smin","realized","realizedAB","sm","additive"), outputName){
  type=match.arg(type) 
#   file_prefix=substr(geno_file,1,nchar(geno_file)-4)

  if (is.character(geno_file)) {
    formatsep <- function(fformat) {
      switch(fformat, 
             csv = ",",
             ctxt = ",",
             stxt = " ",
             ttxt = "\t",
             sctxt = ";")
    }
    fSep = formatsep(fFormat)
    
    #   #   pedigree=as.matrix(read.table(ped_file,sep=fSep,stringsAsFactors=FALSE))
    #   pedigree <- as.matrix(read.table(ped_file, header = TRUE, na.strings = c("NA","."), blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = sep3))
    
    #   geno=read.csv(geno_file,sep=",",row.names=1,stringsAsFactors=FALSE)
    
    geno <- read.table(geno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = fSep)
  } else geno <- geno_file #if geno_file is data frame containing genotype data

  if(type=="additive") {
    rel=A.mat(geno)
  } else {
    mygpC<-create.gpData(geno=geno)
    mygpC$info$codeGeno=TRUE
    rel=kin(mygpC,ret=type)
  }

  write.csv(rel, paste0(outputPath, outputName, ".m"))
  return (list(relFileName=paste0(outputPath, outputName, ".m"), relMat = rel))
}
  

