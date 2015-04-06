#-------------------------------------------------
# QTLdataPrep.R
# DATA PREPARATION for QTL analysis 
# Author: Rose Imee Zhella Morantte
#
# REQUIRED input: 
#	P_data - phenotypic data set (R data format)
#	Gdatafile - text file containing genotypic data
#	Mdatafile - text file containing map information
# 	P_geno - name of genotype variable in P_data
#-------------------------------------------------

QTL.dataprep <- function(P_data, Gdatafile, Mdatafile, P_geno) UseMethod("QTL.dataprep")

QTL.dataprep.default <- function(P_data, Gdatafile, Mdatafile, P_geno) {

	#read in datasets
	G_data <-read.table(Gdatafile,header=FALSE, na.strings = c("NA","."), blank.lines.skip=TRUE, sep = "\t")
	M_data <-read.table(Mdatafile,header=FALSE, na.strings = c("NA","."), blank.lines.skip=TRUE, sep = "\t")
  
  #trim the strings of genotype and Markers ID
	P_data[match(P_geno, names(P_data))]<-trimStrings(as.matrix(P_data[match(P_geno, names(P_data))]))
	G_data[,1]<-trimStrings(G_data[,1])
	M_data[,1]<-trimStrings(M_data[,1])
  
	G_dataMat <- as.matrix(G_data) 
	
	###################################
	#P_data vs G_data
	
	#get genotype and marker "variable" in the data sets
	P_gid <- unique(P_data[match(P_geno,names(P_data))])
	M_mid <- M_data[1]
	colnames(M_mid) <- c("1")
	G_mid <- G_dataMat[1,]
	G_midt <- as.data.frame(t(G_dataMat)[,1])
	rownames(G_midt) <- NULL
	G_gidt <- as.data.frame(t(G_dataMat)[1,])
	colnames(G_gidt) <- "G.id" #replaced "V1"

	##check if there are genotypes in G_data w/c are not in P_data; for displaying (if any)
	G_diffGid <- as.character(as.matrix(setdiff(G_data[,1],as.matrix(P_data[match(P_geno, names(P_data))]))))
	G_diffGidNoSpace<-G_diffGid[which(G_diffGid!="")]
  
	##check if there are genotypes in P_data w/c are not in G_data; for displaying (if any)
	P_diffGid <- as.character(as.matrix(setdiff(as.matrix(P_data[match(P_geno, names(P_data))]), G_data[,1])))
	P_diffGidNoSpace<-P_diffGid[which(P_diffGid!="")]
  	
  ##reduce (if needed) P_data, sort genotypes as in G_data
  P_dataRed <- merge(P_data, G_gidt, by.x = P_geno, by.y = "G.id", sort = TRUE)

  ##reduce (if needed) G_data
  G_dataRed <- merge(G_dataMat, P_gid, by.x = "V1", by.y = P_geno, sort = FALSE)
  G_dataRed <- rbind(G_mid, G_dataRed)
	
  isNewPhenoCreated<-FALSE
	if (length(P_diffGidNoSpace)!=0) {
  	##save new P_data as csv file
  	write.table(P_dataRed,file=paste(getwd(),"/newPhenoData.csv", sep=""), quote = FALSE, sep = ",", row.names = FALSE, col.names=TRUE)
  	isNewPhenoCreated<-TRUE
	}
  
	###################################
	#G_data vs M_data
	
	colnames(M_data) <- c("V1","V2_1","V3_1")
	G_datat <- as.data.frame(t(G_dataRed))
	rownames(G_datat) <- NULL
	
	##check if there are markers in M_data w/c are not in G_data; for displaying (if any)
	M_diffMid <- as.character(as.matrix(setdiff(M_data[,1],G_datat[,1])))
	M_diffMidNoSpace<-M_diffMid[which(M_diffMid!="")]
	
	#reduce, if needed, M_data
	M_dataRed <- merge(M_data, G_midt, by.x = "V1", by.y = names(G_midt)[1], sort = FALSE)
	
	##check if there are markers in G_data w/c are not in M_data; for displaying (if any)
	G_diffMid<-as.character(as.matrix(setdiff(G_datat[,1], M_data[,1])))
	G_diffMidNoSpace<-G_diffMid[which(G_diffMid!="")]
	#rownames(G_diffMid) <- NULL
  
	#reduce G_data
	G_datatRed <- merge(M_mid, G_datat, by.x = "1", by.y = names(G_datat)[1], sort = FALSE)
	G_datatRed <- rbind(G_datat[1,], G_datatRed)
	G_dataRed2 <- t(G_datatRed)

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