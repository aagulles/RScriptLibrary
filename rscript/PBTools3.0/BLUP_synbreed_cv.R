###########
#  BLUP_synbreed_cv
#' Function for performing cross validation for GBLUP
#
#  Parameters:
#' @param gsObject R gpData object
#' @param traitName names of trait
#' @param relData R object of class relationshipMatrix
#' @param sampling sampling strategy to employ, whether "random" or "within popStruc"
#' @param popStruc_file name of file describing the population structure, if samplingStrat is "within popStruc"
#' @param k number of folds, default is 2 
#' @param Rep number of replications
#
###########

# cross validation: BLUP by synbreed

## input
# pheno_file: phenotype file (for multiple traits and one replicate for each trait)
# trait: trait index (or name) used for GS
# rel_file: marker-based relationship matrix (GBLUP) or pedigree-based relationship matrix (PBLUP)
# doCV: get predictive ability of BLUP by cross valiation or not

### phenotype file (csv file)
# rows: individual
# rownames: individual ID  
# columns: trait
# columnnames: trait name
# NAs are acceptable
# all individuals should show up in relationship matrix


## output
# predicted genetic values and corresponding selection accuracy for all individuals
# cross validation information (if doCV=TRUE)

# library(synbreed) 
# library(Matrix)

#(outputPath, gsObject, traitName, relData)
# blupOutCVL[[i]] <- BLUP_synbreed_cv(gsObject, traitNames[i], relData, sampling, popStruc_file, nfolds, nrep)
# BLUP_synbreed_cv<-function(pheno_file,trait=1,rel_file,sampling=c("random","within popStruc"),k=5,Rep=20,popStruc_file=NULL){
BLUP_synbreed_cv <- function(gsObject, traitName, relData, sampling = c("random","within popStruc"), popStruc_file = NULL, k = 5, Rep = 20){

  sampling=match.arg(sampling)
  if ( sampling=="within popStruc" & is.null(popStruc_file) ) 
    stop("no popStruc was given")
  
  popStruc=NULL
  # read phenotype data
#   pheno=as.matrix(read.csv(pheno_file,sep=',',row.names=1,stringsAsFactors=F))
  pheno <- gsObject$pheno
  traitNames = colnames(pheno)
  trait = match(traitName, colnames(pheno))

  if(trait<1 | trait>ncol(pheno))
    stop("Trait ID out of range")
 
  # select trait for GS and move NAs
  phenotype = as.matrix(na.omit(pheno[,trait,1]))
  ##CONT FRI  
  # read relationship matrix
#   U=as.matrix(read.csv(rel_file,sep=',',row.names=1,stringsAsFactors=F))         #Diagonal elements approach 2
  U=as.matrix(relData)
  colnames(U)=rownames(U)
  
  phenotype=as.matrix(phenotype[rownames(phenotype)%in%rownames(U),])
  phenotype=as.matrix(phenotype[order(rownames(phenotype)),])  
  colnames(phenotype)=colnames(pheno)[trait]
  U=U[order(rownames(U)),]
  U=U[,order(colnames(U))]
  
  if(max(diag(U))<=1) U=U*2
  
  genotype=matrix(sample(c(0,2),nrow(U)*2,replace=TRUE,prob=c(0.8,0.2)),nrow(U),2)
  genotype=as.data.frame(genotype)
  rownames(genotype)=rownames(U)
  
  mygpC<-create.gpData(pheno=as.data.frame(phenotype),geno=genotype)
  mygpC$info$codeGeno=TRUE
  
  
  myBLUP<-tryCatch({gpMod(mygpC,model="BLUP",kin=U)},error=function(err){
    # adjust relationship matrix to be positive definite
    U_adjusted<-U
    U_adjusted<-matrix(2*nearPD(U/2,corr=TRUE,keepDiag=TRUE,maxit=3000)$mat@x,ncol(U))
    dimnames(U_adjusted)<-dimnames(U)
    U<-U_adjusted
    return(list(U=U,res=gpMod(mygpC,model="BLUP",kin=U)))
    }) 
  if(!is.null(myBLUP$U)){
    U=myBLUP$U
    sig=myBLUP$res$fit$sigma
  }
  else 
    sig=myBLUP$fit$sigma
  
  if (sampling=="within popStruc") {
    popStruc=as.matrix(na.omit(read.csv(popStruc_file,sep=',',row.names=1,stringsAsFactors=F)))
    
    common=intersect(intersect(rownames(popStruc),rownames(genotype)),rownames(phenotype))
    popStruc=as.matrix(popStruc[rownames(popStruc)%in%common,])
    popStruc=popStruc[order(rownames(popStruc)),1]
    
    genotype=genotype[rownames(genotype)%in%common,]
    genotype=genotype[order(rownames(genotype)),] 
    
    phenotype=as.matrix(phenotype[rownames(phenotype)%in%common,])
    phenotype=as.matrix(phenotype[order(rownames(phenotype)),]) 
    
    
    mygpC<-create.gpData(pheno=as.data.frame(phenotype),geno=genotype)
    mygpC$info$codeGeno=TRUE
    
    myBLUP<-tryCatch({gpMod(mygpC,model="BLUP",kin=U)},error=function(err){
      # adjust relationship matrix to be positive definite
      U_adjusted<-U
      U_adjusted<-matrix(2*nearPD(U/2,corr=TRUE,keepDiag=TRUE,maxit=3000)$mat@x,ncol(U))
      dimnames(U_adjusted)<-dimnames(U)
      U<-U_adjusted
      return(list(U=U,res=gpMod(mygpC,model="BLUP",kin=U)))
    }) 
    if(!is.null(myBLUP$U)){
      U=myBLUP$U
      sig=myBLUP$res$fit$sigma
    }
    else 
      sig=myBLUP$fit$sigma
  }  
  
  cv<-crossVal(mygpC,cov.matrix=list(U),k=k,Rep=Rep,sampling=sampling,popStruc=popStruc,varComp=sig,VC.est="commit")

  return(cv)

}
