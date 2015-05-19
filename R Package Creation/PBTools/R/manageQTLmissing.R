##########################################################
#manageQTLmissing - function for handling missing data
#
# REQUIRED input: 
# crossData - R object of class cross
# deleteMiss - logical; whether to delete missing data or not (imputation)
# cutOff - minimum proportion of missing data per row/column for which the rows/columns will be deleted
#
# OUTPUT: a list with the following elements:
#  crossObj - an R object of class "cross"
#  nMissingRow - number of rows deleted
#  missingRow - rows deleted
#  nMissingCol - number of columns deleted
#  missingCol - columns deleted
#
##########################################################

manageQTLmissing <- function(crossData, deleteMiss, cutOff = NULL) UseMethod("manageQTLmissing") 

manageQTLmissing.default <- function(crossData, deleteMiss, cutOff = NULL) {
  
  # delete missing
  if (deleteMiss) {
#     rownames(gData) <- rownames(pData)
    rownmiss <- nmissing(crossData, what = "ind")
#     rowpmiss <- nmiss/totmar(crossData)
    rowpmiss <- rownmiss/totmar(crossData)
    rowTMiss <- rowpmiss >= cutOff
    nRowMiss <- sum(rowTMiss)
    missGeno <- which(rowTMiss==TRUE)
    
    colnmiss <- nmissing(crossData, what = "mar")
#     colpmiss <- nmiss/nind(crossData)
    colpmiss <- colnmiss/nind(crossData)
    colTMiss <- colpmiss >= cutOff
    nColMiss <- sum(colTMiss)
    missMarkers <- markernames(crossData)[which(colTMiss==TRUE)]
    
    crossDataNM <- drop.markers(crossData, markers = missMarkers)
#     crossDataNM <- crossDataNM(,)
#     crossDataNM <- crossDataNM[,which(rowTooManyMiss==FALSE)]
    crossDataNM <- crossDataNM[,which(rowTMiss==FALSE)]

    #     colNMiss <- NULL
    #     rowNMiss <- NULL
    #     colPctMiss <- NULL
    #     rowPctMiss <- NULL
    #     colTooManyMiss <- NULL
    #     rowTooManyMiss <- NULL
    
    #     for (i in 1:ncol(gData)) {
    #       colNMiss[j] <- sum(is.na(gData[,i]))
    #       colPctMiss[i] <- sum(is.na(gData[,i]))/nrow(gData)  
    #       colTooManyMiss[i] <- colPctMiss[i] >= cutOff
    #     }
    #     
    #     for (j in 1:nrow(gData)) {
    #       rowNMiss[j] <- sum(is.na(gData[j,]))
    #       rowPctMiss[j] <- sum(is.na(gData[j,]))/ncol(gData)
    #       rowTooManyMiss[j] <- rowPctMiss[j] >= cutOff
    #     }
    
    #check marker and genotype ids prior to merging into new cross object
    #     gDataNoMiss <- gData[which(rowTooManyMiss==FALSE),which(colTooManyMiss==FALSE)]
    #     cdsub =crossData[which(rowTooManyMiss==FALSE),which(colTooManyMiss==FALSE)]
    #     crossData[,which(rowTooManyMiss==FALSE)]
    #     crossData[which(rowTooManyMiss==FALSE),] #,which(colTooManyMiss==FALSE)]
    #     P_geno = make.unique(c(colnames(pData), colnames(gDataNoMiss), "rowGData"))[ncol(pData)+ncol(gDataNoMiss)+1]
    #     gData2 <- cbind(rownames(gDataNoMiss),gDataNoMiss)
    #     colnames(gData2)[1] <- P_geno
    #     pData2 <- cbind(rownames(pData),pData)
    #     colnames(pData2)[1] <- P_geno
    #     mData2 <- cbind(rownames(mData),mData)
    #     colnames(mData2)[1] <- "m.Id"
    #     
    #     #matching the three input files
    #     QTLdataL <- QTL.dataprep(P_data = pData2, G_data = gData2, M_data = mData2, P_geno)
    
  } else {
    
    crossDataNM <- fill.geno(crossData, method = "imp")
    nRowMiss <- NA
    missGeno <- NULL
    nColMiss <- NA
    missMarkers <- NULL
    
    # rm(fake.f2)
    #     g1 <- pull.geno(fake.f2)
    #     g2 <- pull.geno(fake.f22)
    
  }
  
  return(list(crossObj = crossData, nMissingRow = nRowMiss, missingRow = missGeno, 
              nMissingCol = nColMiss, missingCol = missMarkers))
              #, cross = crossType, bcGen = bcNum, fGen = fNum))
}