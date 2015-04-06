###########
#  GSDataImputation
#' Function for imputing data
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param geno_file name of genotypic file
#' @param impType type of imputation to use, whether random or family
#' @param pheno_file name of phenotypic file
#' @param familyTrait name of variable to describe family structure, if impType is family
#' @param packageFormat package format used in saving genotype file, whether "synbreed", "rrBLUP" or "BGLR"
  
#  Output
#   genoData - R data containing genotype information
###########

GSDataImputation <- function(outputPath, geno_file, impType = c("random", "family"), pheno_file = NULL, familyTrait = NULL, packageFormat = c("synbreed", "rrBLUP", "BGLR")) {
    
  genoData <- read.table(geno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",")    

  phenoData = NULL
  if (impType == "family") {
    phenoData <- read.table(family_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",")    
  }
  
  genoGp <- create.gpData(geno = genoData, family = phenoData[, familyTrait])

#Z#   capture.output(syn_code<-codeGeno(mygp,impute=TRUE,impute.type="random",nmiss=nmiss,maf=maf,verbose=TRUE),file=paste0(file_prefix,"-QCSummary-",nmiss,"-",maf,".txt"))
  capture.output(syn_code<-codeGeno(genoGp,impute = TRUE,impute.type = impType, nmiss = NULL, maf = NULL,verbose=TRUE),file=paste0(outputPath,"QCSummary_imp.txt"))

  #save using synbreed format
  if (packageFormat == "synbreed") {
  #     syn_code=discard.markers(syn_code, snps[dis])  # remove high-correlated markers #Z# regardless of package format?
      syn_code$geno[syn_code$geno==1]=0        # set heterozygotes to AA
      #Z#   write.csv(syn_code$geno,paste0(file_prefix,"-synbreed-",nmiss,"-",maf,"-",cor_threshold,".csv"))  
      outputName <- paste0(outputPath,"synbreed.csv")
      write.csv(syn_code$geno, outputName)  
    } else if (packageFormat == "rrBLUP") {
      # recode into rrBLUP package
      rrBLUP_code=syn_code$geno
      rrBLUP_code[which(syn_code$geno==0,arr.ind=T)]=1
      rrBLUP_code[which(syn_code$genod==1,arr.ind=T)]=0
      rrBLUP_code[which(syn_code$geno==2,arr.ind=T)]=-1
      #Z#   write.csv(rrBLUP_code,paste0(file_prefix,"-rrBLUP-",nmiss,"-",maf,"-",cor_threshold,".csv"))
      outputName <- paste0(outputPath,"rrBLUP.csv")
      write.csv(rrBLUP_code, outputName)  
    } else if (packageFormat == "BGLR") {
      # recode into BGLR package
      BGLR_code=syn_code$geno
      BGLR_code[which(syn_code$geno==2,arr.ind=T)]=3
      BGLR_code[which(is.na(syn_code$geno))]=2
      #Z#   write.csv(BGLR_code,paste0(file_prefix,"-BGLR-",nmiss,"-",maf,"-",cor_threshold,".csv"))
      outputName <- paste0(outputPath,"BGLR.csv") 
      write.csv(BGLR_code, outputName)
    }
  
  #Z# return(list(synbreed=paste0(file_prefix,"-synbreed-",nmiss,"-",maf,"-",cor_threshold,".csv"),BGLR=paste0(file_prefix,"-BGLR-",nmiss,"-",maf,"-",cor_threshold,".csv"),rrBLUP=paste0(file_prefix,"-rrBLUP-",nmiss,"-",maf,"-",cor_threshold,".csv")))
  #Z# return(list(synbreed=paste0(outputPath,"synbreed.csv"),BGLR=paste0(outputPath,"BGLR.csv"),rrBLUP=paste0(outputPath,"rrBLUP.csv")))
  return(geno_data = syn_code$geno)
}