###########
#  BLUP_synbreed_gv
#' Function for calculating genetic value (GBLUP)
#
#  Parameters:
#' @param gsObject R gpData object
#' @param traitName names of trait
#' @param relData R object of class relationshipMatrix
#
###########

# genetic value calculation: BLUP by synbreed

## input
# pheno_file: phenotype file (for multiple traits and one replicate for each trait)
# trait: trait index (or name) used for GS
# rel_file: marker-based relationship matrix (GBLUP) or pedigree-based relationship matrix (PBLUP)


### phenotype file (csv file)
# rows: individual
# rownames: individual ID  
# columns: trait
# columnnames: trait name
# NAs are acceptable
# all individuals should show up in relationship matrix


## output
# predicted genetic values for all individuals

# library(synbreed) 
# library(Matrix)

#Z# BLUP_synbreed_gv<-function(pheno_file,trait=1,rel_file){
# BLUP_synbreed_gv <- function(outputPath, gsObject, traitName
#                              #for cv
#                              , doCV=FALSE, varCompEst = c("BL", "BRR"), sampling=c("random","within popStruc"), popStruc_file = NULL,
#                              nfolds, nrep
#                              ){

BLUP_synbreed_gv <- function(gsObject, traitName, relData){
  
#   file_prefix=substr(rel_file,1,nchar(rel_file)-4)              # get directory of phenotype file, results are saved at the same directory

  # read phenotype data
#Z#   pheno=as.matrix(read.csv(pheno_file,sep=',',row.names=1,stringsAsFactors=F))
  #not matrix??
  pheno = gsObject$pheno
  traitNames = colnames(pheno)
  trait = match(traitName, colnames(pheno))
  
  if(trait<1 | trait>ncol(pheno))
    stop("Trait ID out of range")
 
  # select trait for GS and move NAs
#Z#   phenotype=as.matrix(na.omit(pheno[,trait]))
  phenotype=as.matrix(na.omit(pheno[, trait, 1]))
    
  # read relationship matrix
#Z#   U=as.matrix(read.csv(rel_file,sep=',',row.names=1,stringsAsFactors=F))         #Diagonal elements approach 2
  U=as.matrix(relData)

  colnames(U)=rownames(U)
  U=U[order(rownames(U)),]
  U=U[,order(colnames(U))]
  
  if(max(diag(U))<=1) U=U*2
  
  phenotype=as.matrix(phenotype[rownames(phenotype)%in%rownames(U),])
  phenotype=as.matrix(phenotype[order(rownames(phenotype)),])  

  # generate NAs for unphenotyped individuals
  NAs=rep(NA,nrow(U)-nrow(phenotype))
  names(NAs)=setdiff(rownames(U),rownames(phenotype))
  phenotype=as.matrix(c(phenotype[,1],NAs))
  phenotype=as.matrix(phenotype[order(rownames(phenotype)),])
    
  mygpC<-create.gpData(pheno=as.data.frame(phenotype))
  mygpC$info$codeGeno=TRUE
  

  myBLUP<-tryCatch({gpMod(mygpC,model="BLUP",kin=U)},error=function(err){
    # adjust relationship matrix to be positive definite
    U_adjusted<-U
    U_adjusted<-matrix(2*nearPD(U/2,corr=TRUE,keepDiag=TRUE,maxit=3000)$mat@x,ncol(U))
    dimnames(U_adjusted)<-dimnames(U)
    U<-U_adjusted
    return(res=gpMod(mygpC,model="BLUP",kin=U))
    })   

  pred1=data.frame(gv=predict(myBLUP,newdata=rownames(phenotype)[which(is.na(phenotype))]))
  pred2=data.frame(gv=myBLUP$fit$predicted)
  gv=as.matrix(rbind(pred1,pred2))
  gv=gv[order(rownames(gv)),]
  
  result=cbind(phenotype[,1],gv)
  
  colnames(result)=c(colnames(pheno)[trait], paste0(colnames(pheno)[trait],"_gv"))

  # save the result
#Z#   write.csv(result,paste0(file_prefix,"-",colnames(pheno)[trait],"-synbreed-GV.csv"))    
  
  return(gv = result)
}
