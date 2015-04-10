###########
#  pedigree_relationship
#' Function for computing pedigree-based relationship matrix
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param ped_file name of genotypic input file
#' @param fFormat format of genotypic input file, whether comma-separated file ("csv"), or text file using comma ("ctxt"), space ("stxt"), tab ("ttxt"), or semi-colon ("sctxt") as delimiter
#' @param type type of relationship matrix to compute, whether "add", "kin","dom", or "gam", for now, use either "add" or "dom"
#' @param outputName name of file to save the relationship matrix into
#' 
#  Output list containing
#   relMat - relationship matrix
#   relFileName - gen_RM$relFileName, 
###########
# compute pedigree-based relationship matrix

# input
## ped_file: pedigree information
## type: method for computing relationship matrix
### "add": methods provided by package pedigreemm
### "kin", "dom", and "gam": methods provided by package synbreed
#### NOTE:
###### methods provided by package pedigreemm and GeneticsPed are same as "add" in synbreed
###### package pedigreemm and GeneticsPed are much faster than "add" in synbreed
###### "kin" is identical with "gam", they are half of "add"

## ped_file requirement
### rows: individuals
### columns: id, par1,par2
### without ancestor: use 0 

# output
## pedigree-based relationship matrix


# library(synbreed)
# library(pedigreemm)

pedigree_relationship<-function(outputPath, ped_file, fFormat = c("csv", "ctxt", "stxt", "ttxt", "sctxt"), type=c("add","kin","dom","gam"), outputName){
  
  type=match.arg(type)
#   file_prefix=substr(ped_file,1,nchar(ped_file)-4)

  formatsep <- function(fformat) {
    switch(fformat, 
           csv = ",",
           ctxt = ",",
           stxt = " ",
           ttxt = "\t",
           sctxt = ";")
  }
  fSep = formatsep(fFormat)
  
#   pedigree=as.matrix(read.table(ped_file,sep=fSep,stringsAsFactors=FALSE))
  pedigree <- as.matrix(read.table(ped_file, header = TRUE, na.strings = c("NA","."), blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = fSep))
    
#   head(pedigree) #
  id=pedigree[,1]
  par1=pedigree[,2]
  par2=pedigree[,3]

  if(type=="dom"){
    ped = create.pedigree(id,par1,par2)
    gp = create.gpData(pedigree=ped)
    rel = kin(gp,ret=type)
  } else{    
    ped = pedigree(sire = par1, dam  = par2, label= id)
    rel<-as.matrix(getA(ped))
    if(type!="add") rel=rel/2
  }
  
  write.table(rel,paste0(outputPath, outputName,".p"),col.names=NA,sep=",")  
  return (list(relMat = rel,
               relFileName=paste0(outputPath, outputName,".p")))
}