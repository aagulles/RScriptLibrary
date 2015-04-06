###########
#  GSDataCheck
#' Function for preparing/matching data for genomic selection
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param geno_file name of genotypic file
#' @param type marker format type
#' @param nmiss maximum proportion of missing values per marker
#' @param maf minimum minor allele frequency (MAF)
#' @param cor_threshold maximum correlation of markers

#  Output
#   genoData - R data containing genotype information
###########


# filter and recode SNP marker into synbreed, rrBLUP and BGLR format

# correlation filtering: my method (fast)

## input
# ori_file:SNP markers
# type: marker format type (1, 2, 3 or 4)
# nmiss=0.1: proportion of missing marker threshold
# maf=0.05: minor allele frequency threshold
# cor_threshold=0.90: correlation threshold
# NOTE:
# for looser filtering, set: nmiss=0.2,maf=0.01,cor_threshold=0.95 

### SNP markers (csv file)
# rows: individual
# rownames: individual ID
# columns: marker
# columnnames: marker ID

## marker format type
# type 1: AA,AB,BB,NA
# type 2: AG, CT,...,NA
# type 3: A/G, C/T,...,NA
# type 4: 0(AA),1(AB),2(BB),NA (already in synbreed code schema)
# type 5: do not contain heterozygous

## QC processes
# filter by monomorphic
# filter by nmiss threshold
# randomly imputation
# filter by maf threshold
# remove high_correlated markers 

## output
# XXX-QCSummary.txt: summary information of QC process
# XXX-QCreport.csv: QC results for all snps

# synbreed: 0(AA),1(AB),2(BB),NA
# rrBLUP: -1(BB),0(AB),1(AA),NA
# BGLR: 0(AA),1(AB),2(NA),3(BB)

# library(synbreed)

# recode<-function(ori_file,type,nmiss=0.1,maf=0.05,cor_threshold=0.90){
GSDataCheck <- function(outputPath, geno_file, type = c(1, 2, 3, 4, 5), nmiss = 0.1, maf = 0.05, cor_threshold = 0.90) {
    
  if(nmiss<0 | nmiss>1)
    stop("'nmiss' have to be in [0,1]")
  
  if(maf<0 | maf>1)
    stop("'maf' have to be in [0,1]")
  
  if(! type %in% c(1,2,3,4,5))
    stop("'type' should be 1, 2, 3, 4, or 5")
  
#   file_prefix=substr(ori_file,1,nchar(ori_file)-4)
    
  # remove monomorphic markers
#Z#   data=as.matrix(read.csv(ori_file,sep=",",row.names=1,stringsAsFactors=FALSE,check.names=FALSE))
  ori_data <- read.table(geno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = ",", check.names = FALSE)
  data = as.matrix(ori_data)

  counts=apply(data,2,table)
  
  # deal with vector or array result for apply
  if(is.array(counts)) counts = as.data.frame(counts)
  if(!is.list(counts)){
    counts = as.list(counts)
  }
  
  lens=unlist(lapply(counts,length))
  idx=which(lens==1)
  monos=NULL
  if(length(idx)){
    SNPs=NULL
    for(i in 1:length(idx))
      SNPs=c(SNPs,names(counts[[idx[i]]]))
    monos=data.frame(SNPname=names(lens)[idx],major=SNPs,minor="/",MAF=0,missing_proportion=NA,heter_proportion=NA,filtered_nmiss_or_maf=NA,filtered_highCorrelation=NA,filtered_monomorphic=TRUE,stringsAsFactors=FALSE)
    data=data[,!colnames(data)%in%names(lens)[idx]]
  }
  
  
#   if(type==1 & length(setdiff(unique(as.vector(data)),c("AA","AB","BB",NA))))
#     stop("input genotype data are in wrong code schema!")
#   
#   if(type==4 & length(setdiff(unique(as.vector(data)),c(0,1,2,NA))))
#     stop("input genotype data are in wrong code schema!")
  
  mygp<-create.gpData(geno=data)
  rm(data)
  gc()

  is.heter<-function(x){substr(x,1,1)!=substr(x,2,2)}

# {
  if(type==1)    
    syn_code<-codeGeno(mygp,label.heter="AB",print.report=TRUE)
  else if(type==2)
    syn_code<-codeGeno(mygp,label.heter=is.heter,print.report=TRUE)
  else if(type==3)
    syn_code<-codeGeno(mygp,label.heter="alleleCoding",print.report=TRUE)
  else if(type==4)
    syn_code<-codeGeno(mygp,label.heter="1",print.report=TRUE)
  else
    syn_code<-codeGeno(mygp,print.report=TRUE) #no label.heter = NULL?
# }

  report=read.table("SNPreport.txt",header=TRUE,stringsAsFactors=FALSE)[,1:4]
  heters=apply(syn_code$geno,2,function(x) sum(x==1,na.rm=TRUE))
  NAs=apply(syn_code$geno,2,function(x) sum(is.na(x)))
  report=data.frame(report,missing_proportion=NAs/nrow(syn_code$geno),heter_proportion=heters/nrow(syn_code$geno),stringsAsFactors=FALSE)
   
  syn_code$geno[syn_code$geno==1]=NA        # detect heterozygotes and set them to NAs
  syn_code$geno[syn_code$geno==0]="AA"
  syn_code$geno[syn_code$geno==2]="BB"


  mygp<-create.gpData(geno=syn_code$geno)
#   capture.output(syn_code<-codeGeno(mygp,impute=TRUE,impute.type="random",nmiss=nmiss,maf=maf,verbose=TRUE),file=paste0(file_prefix,"-QCSummary-",nmiss,"-",maf,".txt"))
  capture.output(syn_code<-codeGeno(mygp,impute=FALSE,nmiss=nmiss,maf=maf,verbose=TRUE),file=paste0(outputPath,"QCSummary.txt"))
  filtered=!report[,1]%in%colnames(syn_code$geno)
  report=data.frame(report,filtered_nmiss_or_maf=filtered,stringsAsFactors=FALSE) # filter markers by nmiss and maf threshold

  snps=colnames(syn_code$geno)

  # compute correlationship matrix and create a relationship link table
  if(ncol(syn_code$geno)<2000){    # for small dataset
#Z#     corMat=cor(syn_code$geno)
    corMat=cor(syn_code$geno, use = "pairwise.complete.obs")
    
    diag(corMat)=0
    iden=apply(corMat, 1, function(x) which(abs(x)>cor_threshold))
    rm(corMat)
    gc()
  }
  else{                           # for large dataset
    corMat<-bigcor(syn_code$geno,verbose=FALSE)
    diag(corMat)=0
    iden1=ffapply(X=corMat, MARGIN=1, AFUN=function(x) which(abs(x)>cor_threshold), FF_RETURN=FALSE, RETURN=TRUE,CFUN="list")
    close(corMat)
    
    lens=unlist(lapply(iden1,length))
    iden=NULL
    blocksize=max(lens)
    for(i in 1:length(iden1)){
      if(length(iden1[[i]]))
        iden=c(iden,iden1[[i]])
      else if(i!=length(iden1))
        iden=c(iden,vector("list",blocksize))
      else
        iden=c(iden,vector("list",ncol(corMat)-blocksize*(length(iden1)-1)))
    }
  }
  
  # deal with vector or array result for apply
  if(is.array(iden)) iden = as.data.frame(iden)
  if(!is.list(iden)){
    iden = as.list(iden)
    iden = lapply(iden, function(x){names(x) <- snps[x]; return(x)})
  }
  
  lens=unlist(lapply(iden,length))
  
  # remove high_correlated markers step by step
  idx_max=which.max(lens)
  dis=NULL
  while(length(iden)){
    if(lens[idx_max]>1){
      dis=c(dis,idx_max)
      for(i in iden[[idx_max]])
        iden[[i]]=setdiff(iden[[i]],idx_max)
      iden[[idx_max]]=setdiff(iden[[idx_max]],iden[[idx_max]])
      lens=unlist(lapply(iden,length))
      idx_max=which.max(lens)
    }
    else{
      rest=which(lens==1)
      pairs=data.frame(snp1=rest,snp2=0,stringsAsFactors=FALSE)
      for(i in 1:length(rest)){
        snp2=iden[[rest[i]]]
        pairs[i,2]=ifelse(snp2%in%pairs[1:(i-1),1],NA,snp2)
      }
      pairs=na.omit(pairs)
      for(i in 1:nrow(pairs)){
        inf1=round(report[report$SNPname==snps[pairs[i,1]],4:5],digits=3)
        inf2=round(report[report$SNPname==snps[pairs[i,2]],4:5],digits=3)
        dis1=ifelse(inf1[,2]>inf2[,2],pairs[i,1],ifelse(inf1[,1]<inf2[,1],pairs[i,1],pairs[i,2]))
        dis=c(dis,dis1)
      }
      iden=NULL
    }
  }

  report=data.frame(report,filtered_highCorrelation=rownames(report)%in%snps[dis],filtered_monomorphic=FALSE,stringsAsFactors=FALSE)
  report=rbind(monos,report)
#Z#   write.csv(report,paste0(file_prefix,"-QCreport-",nmiss,"-",maf,"-",cor_threshold,".csv"),row.names=FALSE)
  write.csv(report,paste0(outputPath,"QCreport.csv"),row.names=FALSE)


  syn_code=discard.markers(syn_code, snps[dis])  # remove high-correlated markers

#   write.csv(syn_code$geno,paste0(file_prefix,"-synbreed-",nmiss,"-",maf,"-",cor_threshold,".csv"))  
  write.csv(syn_code$geno,paste0(outputPath,"genoData_qc.csv"))  

  #Z##Do this after imputation in imputation function
#   syn_code$geno[syn_code$geno==1]=0        # set heterozygotes to AA
#   write.csv(syn_code$geno,paste0(file_prefix,"-synbreed-",nmiss,"-",maf,"-",cor_threshold,".csv"))  
#   
#   # recode into rrBLUP package
#   rrBLUP_code=syn_code$geno
#   rrBLUP_code[which(syn_code$geno==0,arr.ind=T)]=1
#   rrBLUP_code[which(syn_code$genod==1,arr.ind=T)]=0
#   rrBLUP_code[which(syn_code$geno==2,arr.ind=T)]=-1
#   write.csv(rrBLUP_code,paste0(file_prefix,"-rrBLUP-",nmiss,"-",maf,"-",cor_threshold,".csv"))
#   
# # recode into BGLR package
#   BGLR_code=syn_code$geno
#   BGLR_code[which(syn_code$geno==2,arr.ind=T)]=3
#   BGLR_code[which(is.na(syn_code$geno))]=2
#   write.csv(BGLR_code,paste0(file_prefix,"-BGLR-",nmiss,"-",maf,"-",cor_threshold,".csv"))

#Z#   return(list(synbreed=paste0(file_prefix,"-synbreed-",nmiss,"-",maf,"-",cor_threshold,".csv"),BGLR=paste0(file_prefix,"-BGLR-",nmiss,"-",maf,"-",cor_threshold,".csv"),rrBLUP=paste0(file_prefix,"-rrBLUP-",nmiss,"-",maf,"-",cor_threshold,".csv")))
return(geno_data = syn_code$geno)

}


