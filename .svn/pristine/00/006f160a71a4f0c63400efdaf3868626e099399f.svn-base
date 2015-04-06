####################################################################
#  QTLDataPrep
#' Function for preparing data file(s) for QTL analysis 
#
#  Parameters:
#' @param P_data name of phenotypic data set (R data format)
#' @param G_data genotypic data set
#' @param M_data map data set
#' @param P_geno name of genotype variable in P_data
#'
#' @return list with the elements G_diffGid, P_diffGid, M_diffMid, G_diffMid, isNewPhenoCreated, isNewMapCreated, isNewGenoCreated
#  where:
#  G_diffGid - list of genotypes in G_data w/c are not in P_data
#  P_diffGid - list of genotypes in P_data w/c are not in G_data
#  M_diffMid - list of markers in M_data w/c are not in G_data
#  G_diffMid - list of markers in G_data w/c are not in M_data
#  isNewPhenoCreated - logical; whether a new phenotype file is created
#  isNewMapCreated - logical; whether a new map file is created
#  isNewGenoCreated - logical; whether a new genotype file is created
#
#' @author Author: Rose Imee Zhella Morantte
#-------------------------------------------------

QTLDataPrep <- function(P_data, G_data, M_data, P_geno) UseMethod("QTLDataPrep")

QTLDataPrep.default <- function(P_data, G_data, M_data, P_geno) {


  #trim the strings of genotype and Markers ID
  P_data[,match(P_geno, names(P_data))] <- trimStrings(as.matrix(P_data[match(P_geno, names(P_data))]))
  G_data[,1] <- trimStrings(G_data[,1])
  M_data[,1] <- trimStrings(M_data[,1])
  
  ###################################
  #P_data vs G_data
  
  #get genotype and marker "variable" in the data sets
  colnames(G_data)[1] <- P_geno
  P_gid <- unique(P_data[,match(P_geno,names(P_data))])
  M_mid <- M_data[1]
  colnames(M_mid) <- c("1")
  G_mid <- colnames(G_data)[-1] #G_dataMat[1,]
  G_midt <- as.data.frame(G_mid) #t(G_dataMat)[,1])
  rownames(G_midt) <- NULL
  G_gidt <- data.frame(G_data[,1]) #G_gidt <-data.frame(I(G_data[,1]))# G_gidt <- as.data.frame(t(G_dataMat)[1,])
  colnames(G_gidt) <- P_geno #"G.id" #replaced "V1"

  ##check if there are genotypes in G_data w/c are not in P_data; for displaying (if any)
  G_diffGid <- as.character(as.matrix(setdiff(paste(G_gidt[,1]), P_data[,match(P_geno, names(P_data))])))
  G_diffGidNoSpace<-G_diffGid[which(G_diffGid!="")]
  
  ##check if there are genotypes in P_data w/c are not in G_data; for displaying (if any)
  P_diffGid <- as.character(as.matrix(setdiff(P_data[,match(P_geno, names(P_data))], paste(G_gidt[,1]))))
  P_diffGidNoSpace <- P_diffGid[which(P_diffGid!="")]
  
  ##reduce (if needed) P_data, sort genotypes as in G_data
  P_dataRed <- merge(P_data, G_gidt, by = P_geno, sort = FALSE)
  
  ##reduce (if needed) G_data
  G_dataRed <- merge(G_data, P_gid, by = P_geno, sort = FALSE)
  
  isNewPhenoCreated<-FALSE
  if (length(P_diffGidNoSpace)!=0) {
    ##save new P_data as csv file
    write.table(P_dataRed,file=paste(getwd(),"/newPhenoData.csv", sep=""), quote = FALSE, sep = ",", row.names = FALSE, col.names=TRUE)
    isNewPhenoCreated<-TRUE
  }
  
  ###################################
  #G_data vs M_data
  
  colnames(M_data) <- c("V1","V2_1","V3_1")
  G_datat <- as.data.frame(t(G_dataRed)) # t(G_dataRed) ###no -1 row?
  G_datat <- cbind(G_datat, rownames(G_datat))
  ncolGdatat <- dim(G_datat)[2]
  colnames(G_datat)[ncolGdatat] <- "mID"
  
  ##check if there are markers in M_data w/c are not in G_data; for displaying (if any)
  M_diffMid <- as.character(as.matrix(setdiff(M_data[,1], G_datat[,"mID"])))
  M_diffMidNoSpace<-M_diffMid[which(M_diffMid!="")]
  
  #reduce, if needed, M_data
  M_dataRed <- merge(M_data, G_midt, by.x = "V1", by.y = names(G_midt)[1], sort = FALSE)
  
  ##check if there are markers in G_data w/c are not in M_data; for displaying (if any)
  G_diffMid<-as.character(as.matrix(setdiff(G_datat[-1,"mID"], M_data[,1])))
  G_diffMidNoSpace<-G_diffMid[which(G_diffMid!="")]
  
  #reduce G_data
  G_dataRed <- merge(M_mid, G_datat[-1,], by.x = "1", by.y = "mID", sort = FALSE)
  G_dataRed2 <- as.data.frame(t(G_dataRed))
  rownames(G_dataRed2)[1] <- P_geno
  rownames(G_dataRed2)[2:ncolGdatat] <- t(G_datat[1,c(1:ncolGdatat-1)])
  
  isNewMapCreated<-FALSE
  if (length(M_diffMidNoSpace)!=0) {
    ##save new M_data
    write.table(M_dataRed,file=paste(getwd(), "/newMapData.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE)
    isNewMapCreated<-TRUE
  }
  
  isNewGenoCreated<-FALSE
  if (length(G_diffMidNoSpace)!=0 || length(G_diffGidNoSpace)!=0) {
    #save new G_data 
    write.table(G_dataRed2,file=paste(getwd(), "/newGenoData.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE)
    isNewGenoCreated<-TRUE
  }
  
  return(list(G_diffGid = G_diffGid,
              P_diffGid = P_diffGid,
              M_diffMid = M_diffMid,
              G_diffMid = G_diffMid,
              isNewPhenoCreated =isNewPhenoCreated,
              isNewMapCreated =isNewMapCreated,
              isNewGenoCreated =isNewGenoCreated))
}