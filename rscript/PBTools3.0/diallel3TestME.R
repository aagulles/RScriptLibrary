#----------------------------------------------------------------
# This function performs multi-environment analysis for Diallel 3 experiments using Griffing Method                          
# (F1's + Reciprocals)                       
#
# ARGUMENTS:
# design - experiment design
# data - name of data frame
# respvar - vector of response variables
# p1 - string, name of p1 factor
# p2 - string, name of p2 factor
# rep - string, name of rep factor
# block - string, name of block factor
# row - string, name of row factor
# column - string, name of column factor
# cross - logical
# individual - string, name of individual factor
# environment - string, name of environment factor
# alpha - level of significance
#          
# SAS Script Created by: MS Kang                            
# R Script Created by: Violeta Bartolome   09.2011
# R Script Modified by: Nellwyn Sales, Guoyou Ye 
#----------------------------------------------------------------

diallel3TestME <-function(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, p1, p2, rep=NULL, block=NULL, row=NULL, column=NULL, cross=TRUE, individual=NULL, environment, alpha=0.05) UseMethod("diallel3TestME")
  
diallel3TestME.default <-function(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, p1, p2, rep=NULL, block=NULL, row=NULL, column=NULL, cross=TRUE, individual=NULL, environment, alpha=0.05) {
  
  options(show.signif.stars=FALSE)
  data <- eval(parse(text = data))
  library(doBy)
	
  # --- trim the strings --- #
  respvar <- trimStrings(respvar)
  p1 <- trimStrings(p1)
  p2 <- trimStrings(p2)
  if (!is.null(block)) {block <- trimStrings(block) }
  if (!is.null(rep)) {rep <- trimStrings(rep) }
  if (!is.null(row)) {row <- trimStrings(row) }
  if (!is.null(column)) {column <- trimStrings(column) }
  if (!is.null(individual)) {individual <- trimStrings(individual) }
  environment <-trimStrings(environment)
  alpha <- trimStrings(alpha)
	
  # --- create titles --- #
  if (design == "CRD") { designName<-"CRD"}
  if (design == "RCB") { designName<-"RCB"}
  if (design == "Alpha") { designName<-"ALPHA-LATTICE"}
  if (design == "RowColumn") { designName<-"ROW-COLUMN"}
  
  if (cross) {parentsType<-"CROSS"
  } else {parentsType<-"SELF"}
  cat("\nMULTIPLE ENVIRONMENT ANALYSIS\n")
  cat("\nDIALLEL ANALYSIS: GRIFFING METHOD III IN ",designName, " (", parentsType, ")\n", sep="")
  
  # --- define factors --- #
  data[,match(p1, names(data))] <- factor(trimStrings(data[,match(p1, names(data))]))
  data[,match(p2, names(data))] <- factor(trimStrings(data[,match(p2, names(data))]))
  data[,match(environment, names(data))] <- factor(trimStrings(data[,match(environment, names(data))]))
  if (design == "RCB") {data[,match(block, names(data))] <- factor(trimStrings(data[,match(block, names(data))]))  }
  if (design == "Alpha") {
    data[,match(rep, names(data))] <- factor(trimStrings(data[,match(rep, names(data))]))  
    data[,match(block, names(data))] <- factor(trimStrings(data[,match(block, names(data))]))
  }
  if (design == "RowColumn") {
    data[,match(rep, names(data))] <- factor(trimStrings(data[,match(rep, names(data))]))  
    data[,match(row, names(data))] <- factor(trimStrings(data[,match(row, names(data))]))
    data[,match(column, names(data))] <- factor(trimStrings(data[,match(column, names(data))]))
  }
	
  #--- CREATION OF Matrices for GCA, SCA, REC, MAT, NORM ---#
  p <- length(unique(c(levels(data[,match(p1, names(data))]), levels(data[,match(p2, names(data))]))))
  i <- rep(1:p, each=p)
  j <- rep(1:p,p)
  
  ID <- data.frame(i,j)
  
  #--- for method 3 ---#
  M3 <- subset(ID, i != j)
  t <- nrow(M3)
  
  #--- GCA estimable effects ---#
  GCA <- as.matrix(rep(0,t*(p-1)),nrow=t, ncol=(p-1))
  dim(GCA) <- c(t,p-1)
  for (o in 1:t)
    for (n in 1:(p-1))  { 
      GCA[o,n] <- (as.logical(M3[o,"i"]==n)-
        as.logical(M3[o,"i"]==p))+
        (as.logical(M3[o,"j"]==n)-
        as.logical(M3[o,"j"]==p))  }
  colnames(GCA) <- paste("GCA", 1:(p-1), sep="")
  
  #--- SCA estimable effects ---#
  sca.col <- as.data.frame(subset(ID, (i<=(p-3) & (j >= i+1 & j<=(p-1)))))
  ncol <- nrow(sca.col)
  
  SCA <- as.matrix(rep(0,t*ncol),nrow=t, ncol=ncol)
  dim(SCA) <- c(t,ncol)
  
  for (o in 1:t)
    for (n in 1:ncol)  {
      SCA[o,n] <- as.logical(M3[o,"i"]==sca.col[n, "i"])*
        as.logical(M3[o,"j"]==sca.col[n, "j"])-
        (as.logical(M3[o,"i"]==sca.col[n, "i"])+
        as.logical(M3[o,"i"]==sca.col[n, "j"]))*
        as.logical(M3[o,"j"]==p)+
        as.logical(M3[o,"j"]==sca.col[n, "i"])*
        as.logical(M3[o,"i"]==sca.col[n, "j"])-
        as.logical(M3[o,"i"]==p)*
        (as.logical(M3[o,"j"]==sca.col[n, "i"])+
        as.logical(M3[o,"j"]==sca.col[n, "j"]))
      if (((M3[o,"i"]>=(p-2))&(M3[o,"j"]>=(p-1)))|((M3[o,"i"]>=(p-1))&(M3[o,"j"]>=(p-2))))   {
        SCA[o,n] <- -1*as.logical(M3[o,"i"]==(p-2))*
          as.logical(M3[o,"j"]==(p-1))+
          as.logical(M3[o,"i"]>=(p-2))*
          as.logical(M3[o,"j"]==p)*
          as.logical(M3[o,"i"]!=sca.col[n,"j"])-
          as.logical(M3[o,"j"]==(p-2))*
          as.logical(M3[o,"i"]==(p-1))+
          as.logical(M3[o,"j"]>=(p-2))*
          as.logical(M3[o,"i"]==p)*
          as.logical(M3[o,"j"]!=sca.col[n,"j"])  
      }
    }
  colnames(SCA) <- paste("SCA",sca.col$"i", sca.col$j, sep="")
  
  #--- Reciprocal estimable effects ---#
  rec.col <- as.data.frame(subset(ID, (i<=(p-1) & (j >= i+1 & j<=p))))
  ncol <- nrow(rec.col)
  
  REC <- as.matrix(rep(0,t*ncol),nrow=t, ncol=ncol)
  dim(REC) <- c(t,ncol)
  
  for (o in 1:t)
    for (n in 1:ncol)  {
      REC[o,n] <- as.logical(M3[o,"i"]==rec.col[n,"i"])*
        as.logical(M3[o,"j"]==rec.col[n,"j"])-
        as.logical(M3[o,"j"]==rec.col[n,"i"])*
        as.logical(M3[o,"i"]==rec.col[n,"j"])   }
  colnames(REC) <- paste("REC", rec.col$i,rec.col$j,sep="")
  
  #--- Maternal estimable effects ---#
  MAT <- as.matrix(rep(0,t*(p-1)),nrow=t, ncol=(p-1))
  dim(MAT) <- c(t,p-1)
  
  for (o in 1:t)
    for (n in 1:(p-1))  { 
      MAT[o,n] <- as.logical(M3[o,"i"]==n)+
        as.logical(M3[o,"j"]==p)-
        as.logical(M3[o,"j"]==n)-
        as.logical(M3[o,"i"]==p)   
    }
  colnames(MAT) <- paste("MAT", 1:(p-1), sep="")
  
  #--- NonMaternal estimable effects ---#
  nonm.col <- as.data.frame(subset(ID, (i<=(p-2) & (j >= i+1 & j<=(p-1)))))
  ncol <- nrow(nonm.col)
  
  NONM <- as.matrix(rep(0,t*ncol),nrow=t, ncol=ncol)
  dim(NONM) <- c(t,ncol)
  
  for (o in 1:t)
    for (n in 1:ncol)  {
      NONM[o,n] <- (as.logical(M3[o,"i"]==nonm.col[n,"i"])*
        as.logical(M3[o,"j"]==nonm.col[n,"j"]))-
        as.logical(M3[o,"i"]==nonm.col[n,"j"])*
        as.logical(M3[o,"j"]==nonm.col[n,"i"])+
        ((as.logical(M3[o,"i"]==nonm.col[n,"j"])-
        as.logical(M3[o,"i"]==nonm.col[n,"i"]))*
        as.logical(M3[o,"j"]==p))+
        (as.logical(M3[o,"i"]==p)*
        (as.logical(M3[o,"j"]==nonm.col[n,"i"])-
        as.logical(M3[o,"j"]==nonm.col[n,"j"])))     
    }
  colnames(NONM) <- paste("NONM", nonm.col$i,nonm.col$j,sep="")
  # end of CREATION OF GCA, SCA, REC, MAT, NONM 
  
	result <- list()
	for (z in (1:length(respvar))) {
		result[[z]] <- list()
		cat("\n-----------------------------")
		cat("\nRESPONSE VARIABLE: ", respvar[z], "\n", sep="")
		cat("-----------------------------\n")
		
		if (design == "CRD") {
			vec<-c(match(p1, names(data)), match(p2, names(data)), match(environment, names(data)), match(respvar[z], names(data)))
			temp2.data <-subset(data, select=vec)  
		}
		if (design == "RCB") {
			vec<-c(match(p1, names(data)), match(p2, names(data)), match(environment, names(data)), match(block, names(data)), match(respvar[z], names(data)))
			temp2.data <-subset(data, select=vec)
		}
		if (design == "Alpha") {
		  vec<-c(match(p1, names(data)), match(p2, names(data)), match(environment, names(data)), match(rep, names(data)), match(block, names(data)), match(respvar[z], names(data)))
		  temp2.data <-subset(data, select=vec)  
		}
		if (design == "RowColumn") {
		  vec<-c(match(p1, names(data)), match(p2, names(data)), match(environment, names(data)), match(rep, names(data)), match(row, names(data)), match(column, names(data)), match(respvar[z], names(data)))
		  temp2.data <-subset(data, select=vec)  
		}
		temp.data <- temp2.data[(is.na(temp2.data[,match(respvar[z], names(temp2.data))]) == FALSE),]

		if ((nrow(temp.data)/nrow(data)) < 0.80) {
		  cat("\nToo many missing observations. Cannot proceed with the analysis.\n\n")
		  result[[z]]$canEstimateMissing<-"FALSE"
		  next
		} else {
		  # --- create new column containing treatment combinations --- #
		  temp.data$crossEnv<-factor(paste(temp.data[,p1], ":", temp.data[,p2], ":", temp.data[,environment], sep=""))
		  temp.data<-temp.data[order(temp.data$crossEnv),]
		  
		  # --- compute harmonic mean that will be used later in the estimation of genetic variances --- #
		  lengthPerCross<-tapply(temp.data[,respvar[z]], temp.data$crossEnv, length)
		  repHarmonicMean<-1/mean((1/lengthPerCross), na.rm=TRUE)
		  
		  e <- nlevels(temp.data[,match(environment, names(temp.data))])
		  if (design == "CRD") {
		    # --- add column Rep --- #
		    temp.data<-data.frame(temp.data, Rep=sequence(lengthPerCross))
		    
		    nlevelsRep<-max(lengthPerCross, na.rm=TRUE) 
		  }
		  if (design == "RCB") {nlevelsRep <- nlevels(temp.data[,match(block, names(temp.data))]) }
		  if (design == "Alpha" || design == "RowColumn") {nlevelsRep <- nlevels(temp.data[,match(rep, names(temp.data))]) }
		  
		  nBalance<-(p*p*e*nlevelsRep)-(p*nlevelsRep*e)
		  temp.data<-temp.data[-c(match("crossEnv", names(temp.data)))]
      
		  # --- check if max of lengthPerCross is equal to nlevelsRep --- #
		  if (max(lengthPerCross, na.rm=TRUE) > nlevelsRep) {
		    if (design == "RCB") {
		      blockLabelError <- paste("The number of levels of the blocking factor is", nlevelsRep, "but at least one treatment combination is replicated", max(lengthPerCross, na.rm=TRUE), "times in one environment.")
		      blockLabelError2 <- paste("Please check if the column for block is properly labeled.")
		    } else {
		      blockLabelError <- paste("The number of levels of the replicate factor is", nlevelsRep, "but at least one treatment combination is replicated", max(lengthPerCross, na.rm=TRUE), "times in one environment.")
		      blockLabelError2 <- paste("Please check if the column for replicate is properly labeled.")
		    }
		    result[[z]]$canEstimateMissing <-"FALSE"
		    cat("\n ERROR:", blockLabelError)
		    cat("\n       ", blockLabelError2, "\n\n\n")
		    next
		  }
		  
		  # --- data summary --- #
		  #funcTrialSum <- class.information2(names(temp.data),temp.data)
		  funcTrialSum <- class.information(names(temp.data),temp.data)
		  cat("\nDATA SUMMARY: ","\n\n", sep="")
		  print(funcTrialSum)
		  cat("\n Number of observations read: ",nrow(data), sep="")
		  cat("\n Number of observations used: ",nrow(temp.data), sep="")
		  missingObs<-nBalance-nrow(temp.data)
		  
		  if (nBalance<nrow(temp.data)) {
		    cat("\n\n***\nERROR: The number of observations read in the data exceeds the size of a balanced data.\n       Please check if the column for block/replicate is properly labeled.\n***\n\n")
		  } else {
		    # --- check if raw data is balanced. If not, generate estimates for missing values --- #
		    canEstimateMissing <- TRUE
		    cat("\n\n\nANOVA TABLE: (Crosses = ", p1,":", p2, ")\n\n", sep="")
		    if ((nrow(temp.data)/nBalance) >= 0.90) {
		      
		      # estimate missing values per environment
		      dataForAnova<-NULL
		      for (x in (1:e)) {
		        if (design == "CRD") {
		          tempSplit<-temp.data[temp.data[,environment]==levels(temp.data[,match(environment, names(temp.data))])[x],]
		          tempDataForAnova<-tempSplit[,c(p1, p2, "Rep", respvar[z])]
		          tempDataForAnova[,match("Rep", names(tempDataForAnova))] <- factor(trimStrings(tempDataForAnova[,match("Rep", names(tempDataForAnova))]))
		          
		          # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
		          outRecode<-recodeDiallelData(design="diallel3", data=tempDataForAnova, p1=p1, p2=p2, rep="Rep", respvar=respvar[z])
		          balancedData<-outRecode$tempData
		          codingGuide<-outRecode$newCoding
		          
		          # --- estimate missing values --- #
		          tempDataForAnova<-estimateMissingData(design="CRD", data=balancedData, respvar[z], "newCodeP1", "newCodeP2", "Rep")
		          tempDataForAnova<-tempDataForAnova[,c("newCodeP1", "newCodeP2", "Rep", respvar[z])]
		        }
		        if (design == "RCB") {
		          tempSplit<-temp.data[temp.data[,environment]==levels(temp.data[,match(environment, names(temp.data))])[x],]
		          tempDataForAnova<-tempSplit[,c(p1, p2, block, respvar[z])]
		          
		          # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
		          outRecode<-recodeDiallelData(design="diallel3", data=tempDataForAnova, p1=p1, p2=p2, rep=block, respvar=respvar[z])
		          balancedData<-outRecode$tempData
		          codingGuide<-outRecode$newCoding
		          
		          # --- estimate missing values --- #
		          tempDataForAnova<-estimateMissingData(design="RCB", data=balancedData, respvar[z], "newCodeP1", "newCodeP2", block)
		          tempDataForAnova<-tempDataForAnova[,c("newCodeP1", "newCodeP2", block, respvar[z])]
		        }
		        if (design == "Alpha") {
		          temp2Split<-temp2.data[temp2.data[,environment]==levels(temp2.data[,match(environment, names(temp2.data))])[x],]
		          tempSplit<-temp.data[temp.data[,environment]==levels(temp.data[,match(environment, names(temp.data))])[x],]
		          if (nrow(tempSplit)==nBalance/e) {
		            tempDataForAnova<-tempSplit[,c(p1, p2, rep, block, respvar[z])]
		            
		            # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
		            outRecode<-recodeDiallelData(design="diallel3", data=tempDataForAnova, p1=p1, p2=p2, rep=rep, block=block, respvar=respvar[z])
		            tempDataForAnova<-outRecode$tempData
		            codingGuide<-outRecode$newCoding
		            
		          } else {
		            if (nrow(temp2Split)==nBalance/e) {
		              tempDataForAnova<-temp2Split[,c(p1, p2, rep, block, respvar[z])]
		              tempDataForAnova<-estimateNA(design="Alpha", fullData=tempDataForAnova, respvar[z], p1, p2, rep, block, row=NULL, column=NULL)
		              
		              # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
		              outRecode<-recodeDiallelData(design="diallel3", data=tempDataForAnova, p1=p1, p2=p2, rep=rep, block=block, respvar=respvar[z])
		              tempDataForAnova<-outRecode$tempData
		              codingGuide<-outRecode$newCoding
		              
		            } else {
		              tempDataForAnova<-temp2Split[,c(p1, p2, rep, block, respvar[z])]
		              canEstimateMissing<-FALSE
		              stop("The dataset has missing row(s). PBTools cannot supply the value(s) of block factor for the missing row(s).")
		            }
		          }
		        }
		        if (design == "RowColumn") {
		          temp2Split<-temp2.data[temp2.data[,environment]==levels(temp2.data[,match(environment, names(temp2.data))])[x],]
		          tempSplit<-temp.data[temp.data[,environment]==levels(temp.data[,match(environment, names(temp.data))])[x],]
		          if (nrow(tempSplit)==nBalance/e) {
		            tempDataForAnova<-tempSplit[,c(p1, p2, rep, row, column, respvar[z])]
		            
		            # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
		            outRecode<-recodeDiallelData(design="diallel3", data=tempDataForAnova, p1=p1, p2=p2, rep=rep, block=NULL, row=row, column=column, respvar=respvar[z])
		            tempDataForAnova<-outRecode$tempData
		            codingGuide<-outRecode$newCoding
		            
		          } else {
		            if (nrow(temp2Split)==nBalance/e) {
		              tempDataForAnova<-temp2Split[,c(p1, p2, rep, row, column, respvar[z])]
		              tempDataForAnova<-estimateNA(design="RowColumn", fullData=tempDataForAnova, respvar[z], p1, p2, rep, block=NULL, row, column)
		              
		              # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
		              outRecode<-recodeDiallelData(design="diallel3", data=tempDataForAnova, p1=p1, p2=p2, rep=rep, block=NULL, row=row, column=column, respvar=respvar[z])
		              tempDataForAnova<-outRecode$tempData
		              codingGuide<-outRecode$newCoding
		              
		            } else {
		              tempDataForAnova<-temp2Split[,c(p1, p2, rep, row, column, respvar[z])]
		              canEstimateMissing<-FALSE
		              stop("The dataset has missing row(s). PBTools cannot supply the values of row and column factors for the missing row(s).")
		            }
		          }
		        }
		        tempDataForAnova$EnvIndex<-levels(temp.data[,match(environment, names(temp.data))])[x]
		        dataForAnova<-rbind(dataForAnova,tempDataForAnova)
		      }
		      colnames(dataForAnova)[match("EnvIndex",names(dataForAnova))]<-environment
		      dataForAnova[,match(environment, names(dataForAnova))] <- factor(dataForAnova[,match(environment, names(dataForAnova))])
		      
		      if (nrow(temp.data) == nBalance) {
		        anovaRemark <- "REMARK: Raw dataset is balanced."
		      } else {
		        if (canEstimateMissing) {anovaRemark <- "REMARK: Raw data and estimates of the missing values are used."
		        } else {anovaRemark <- "REMARK: Cannot estimate missing values because some treatment combinations are not specified in the data."}
		      }
          
		      result[[z]]$canEstimateMissing <- toString(canEstimateMissing)
		      
		      r<-nlevelsRep
		      
		      if (canEstimateMissing==TRUE) {
		        # --- one-step ANOVA --- #
		        #if (design == "CRD") {myformula1 <- paste(respvar[z], " ~ ", environment, " + ", p1, ":", p2, " + ", environment, ":", p1, ":", p2, sep = "") }
		        #if (design == "RCB") {myformula1 <- paste(respvar[z], " ~ ", environment, " + ", block, ":", environment," + ", p1, ":", p2, " + ", environment, ":", p1, ":", p2, sep = "") }
		        
		        # --- two-step ANOVA --- #
		        if (design == "CRD") {
		          # --- step 1 ANOVA ---#
		          myformula1 <- paste(respvar[z], " ~ ", environment, sep = "")
		          out1 <- anova(aov(formula(myformula1), data = dataForAnova))
		          totalSS1<-sum(out1[,2])
		          totalDF1<-sum(out1[,1])
		          
		          # --- step 2 ANOVA ---#
		          # --- compute the means --- #
		          myformula2 <- paste(respvar[z], " ~ newCodeP1 + newCodeP2", " + ", environment, sep = "")
		          meanData<-summaryBy (formula(myformula2),data=dataForAnova,FUN=c(mean))
		          
		          myformula3 <- paste(respvar[z], ".mean ~ newCodeP1:newCodeP2", " + ", environment, sep = "")
		          out2<-anova(aov(formula(myformula3),data=meanData))
		          
		          DFvector2<-out2[,1]
		          SSvector2<-r*out2[,2]
		          totalDF2<-sum(out2[,1])
		          totalSS2<-sum(r*out2[,2])
		          
		          # --- construct ANOVA table --- #
		          ESS<-totalSS1-totalSS2
		          EDF<-totalDF1- totalDF2 - missingObs
		          Df<-c(DFvector2,EDF)
		          SS<-c(SSvector2,ESS)
		          MS<-SS/Df
		          F_value<-c(MS[1]/MS[4],MS[2]/MS[3],MS[3]/MS[4],MS[4]/MS[4])
		          p_value<-pf(F_value,Df,c(Df[4],Df[3],Df[4],Df[4]),lower.tail = F)
		          AOV<-cbind(Df,SS, MS,F_value,p_value)
		          
		          row.names(AOV)<-c(environment,"Crosses",paste("Crosses x ",environment,sep=""),"Residuals")
		          AOV[4,4:5]<-NA
		          AOV <- as.data.frame(AOV)
		        }      
		        
		        if (design == "RCB") {
		          # --- step 1 ANOVA ---#
		          myformula1 <- paste(respvar[z], " ~ ", block, ":", environment," + ", environment, sep = "")
		          out1 <- anova(aov(formula(myformula1), data = dataForAnova))
		          blockEnv <- paste(block, ":", environment,sep="")
		          repEnvDF <- out1[rownames(out1)==blockEnv,"Df"]
		          repEnvSS <- out1[rownames(out1)==blockEnv,"Sum Sq"]
		          totalSS1<-sum(out1[,2])
		          totalDF1<-sum(out1[,1])
		          
		          # --- step 2 ANOVA ---#
		          # --- compute the means --- #
		          myformula2 <- paste(respvar[z], " ~ newCodeP1 + newCodeP2"," + ", environment, sep = "")
		          meanData<-summaryBy (formula(myformula2),data=dataForAnova,FUN=c(mean))
		          
		          myformula3 <- paste(respvar[z], ".mean ~ newCodeP1:newCodeP2"," + ", environment, sep = "")
		          out2<-anova(aov(formula(myformula3),data=meanData))
		          
		          envDF<-out2[rownames(out2)==environment,"Df"]
		          envSS<-r*out2[rownames(out2)==environment,"Sum Sq"]
		          DFvector2<-out2[2:3,1]
		          SSvector2<-r*out2[2:3,2]
		          totalDF2<-sum(out2[,1])
		          totalSS2<-sum(r*out2[,2])
		          
		          # --- construct ANOVA table --- #
		          ESS<-totalSS1-totalSS2-repEnvSS
		          EDF<-totalDF1- totalDF2-repEnvDF - missingObs
		          Df<-c(envDF,repEnvDF,DFvector2,EDF)
		          SS<-c(envSS,repEnvSS,SSvector2,ESS)
		          MS<-SS/Df
		          F_value<-c(MS[1]/MS[2],MS[2]/MS[5],MS[3]/MS[4],MS[4]/MS[5],MS[5]/MS[5])
		          p_value<-pf(F_value,Df,c(Df[2],Df[5],Df[4],Df[5],Df[5]),lower.tail = F)
		          AOV<-cbind(Df,SS, MS,F_value,p_value)
		          
		          row.names(AOV)<-c(environment,paste(block,"(",environment,")",sep=""),"Crosses",paste("Crosses x ",environment,sep=""),"Residuals")
		          AOV[5,4:5]<-NA
		          AOV <- as.data.frame(AOV)
		        }
		        
		        if (design == "Alpha") {
		          # --- step 1 ANOVA ---#
		          myformula1 <- paste(respvar[i], " ~ ", environment, " + ", environment, ":", rep, " + ", environment, ":", rep, ":", block, sep = "")
		          out1 <- anova(aov(formula(myformula1), data = dataForAnova))
		          envRepDF <- out1[rownames(out1)==paste(environment, ":", rep,sep=""),"Df"]
		          envRepSS <- out1[rownames(out1)==paste(environment, ":", rep,sep=""),"Sum Sq"]
		          envRepBlockDF<-out1[rownames(out1)==paste(environment, ":", rep,":", block, sep=""),"Df"]
		          envRepBlockSS<-out1[rownames(out1)==paste(environment, ":", rep,":", block, sep=""),"Sum Sq"]
		          totalSS1<-sum(out1[,2])
		          totalDF1<-sum(out1[,1])
		          
		          # --- step 2 ANOVA ---#
		          # --- compute the means --- #
		          myformula2 <- paste(respvar[z], " ~ newCodeP1 + newCodeP2"," + ", environment, sep = "")
		          meanData<-summaryBy (formula(myformula2),data=dataForAnova,FUN=c(mean))
		          
		          myformula3 <- paste(respvar[z], ".mean ~ newCodeP1:newCodeP2", " + ", environment, sep = "")
		          out2<-anova(aov(formula(myformula3),data=meanData))
		          
		          envDF<-out2[rownames(out2)==environment,"Df"]
		          envSS<-r*out2[rownames(out2)==environment,"Sum Sq"]
		          DFvector2<-out2[2:3,1]
		          SSvector2<-r*out2[2:3,2]
		          totalDF2<-sum(out2[,1])
		          totalSS2<-sum(r*out2[,2])
		          
		          # --- construct ANOVA table --- #
		          ESS<-totalSS1-totalSS2-envRepSS-envRepBlockSS
		          EDF<-totalDF1-totalDF2-envRepDF-envRepBlockDF - missingObs
		          Df<-c(envDF,envRepDF,envRepBlockDF,DFvector2,EDF)
		          SS<-c(envSS,envRepSS,envRepBlockSS,SSvector2,ESS)
		          MS<-SS/Df
		          F_value<-c(MS[1]/MS[2],MS[2]/MS[6],MS[3]/MS[6],MS[4]/MS[5],MS[5]/MS[6],MS[6]/MS[6])
		          p_value<-pf(F_value,Df,c(Df[2],Df[6],Df[6],Df[5],Df[6],Df[6]),lower.tail = F)
		          AOV<-cbind(Df,SS, MS,F_value,p_value)
		          
		          row.names(AOV)<-c(environment,paste(environment,":",rep,sep=""),paste(environment,":",rep,":",block,sep=""),"Crosses",paste("Crosses x ",environment,sep=""),"Residuals")
		          AOV[6,4:5]<-NA
		          AOV <- as.data.frame(AOV)
		        }
		        
		        if (design == "RowColumn") {
		          # --- step 1 ANOVA ---#
		          myformula1 <- paste(respvar[i], " ~ ", environment, " + ", environment, ":", rep, " + ", environment, ":", rep, ":", row, " + ", environment, ":", rep, ":", column, sep = "") 
		          out1 <- anova(aov(formula(myformula1), data = dataForAnova))
		          envRepDF <- out1[rownames(out1)==paste(environment, ":", rep,sep=""),"Df"]
		          envRepSS <- out1[rownames(out1)==paste(environment, ":", rep,sep=""),"Sum Sq"]
		          envRepRowDF<-out1[rownames(out1)==paste(environment, ":", rep,":", row, sep=""),"Df"]
		          envRepRowSS<-out1[rownames(out1)==paste(environment, ":", rep,":", row, sep=""),"Sum Sq"]
		          envRepColumnDF<-out1[rownames(out1)==paste(environment, ":", rep,":", column, sep=""),"Df"]
		          envRepColumnSS<-out1[rownames(out1)==paste(environment, ":", rep,":", column, sep=""),"Sum Sq"]
		          totalSS1<-sum(out1[,2])
		          totalDF1<-sum(out1[,1])
		          
		          # --- step 2 ANOVA ---#
		          # --- compute the means --- #
		          myformula2 <- paste(respvar[z], " ~ newCodeP1 + newCodeP2", " + ", environment, sep = "")
		          meanData<-summaryBy (formula(myformula2),data=dataForAnova,FUN=c(mean))
		          
		          myformula3 <- paste(respvar[z], ".mean ~ newCodeP1:newCodeP2"," + ", environment, sep = "")
		          out2<-anova(aov(formula(myformula3),data=meanData))
		          
		          envDF<-out2[rownames(out2)==environment,"Df"]
		          envSS<-r*out2[rownames(out2)==environment,"Sum Sq"]
		          DFvector2<-out2[2:3,1]
		          SSvector2<-r*out2[2:3,2]
		          totalDF2<-sum(out2[,1])
		          totalSS2<-sum(r*out2[,2])
		          
		          # --- construct ANOVA table --- #
		          ESS<-totalSS1-totalSS2-envRepSS-envRepRowSS-envRepColumnSS
		          EDF<-totalDF1-totalDF2-envRepDF-envRepRowDF-envRepColumnDF - missingObs
		          Df<-c(envDF,envRepDF,envRepRowDF,envRepColumnDF,DFvector2,EDF)
		          SS<-c(envSS,envRepSS,envRepRowSS,envRepColumnSS,SSvector2,ESS)
		          MS<-SS/Df
		          F_value<-c(MS[1]/MS[2],MS[2]/MS[7],MS[3]/MS[7],MS[4]/MS[7],MS[5]/MS[6],MS[6]/MS[7],MS[7]/MS[7])
		          p_value<-pf(F_value,Df,c(Df[2],Df[7],Df[7],Df[7],Df[6],Df[7],Df[7]),lower.tail = F)
		          AOV<-cbind(Df,SS, MS,F_value,p_value)
		          
		          row.names(AOV)<-c(environment,paste(environment,":",rep,sep=""),paste(environment,":",rep,":",row,sep=""),paste(environment,":",rep,":",column,sep=""),"Crosses",paste("Crosses x ",environment,sep=""),"Residuals")
		          AOV[7,4:5]<-NA
		          AOV <- as.data.frame(AOV)
		        }
		        
		        AOV_print <- formatAnovaTable(AOV)
		        print(AOV_print)
		        cat("-------\n")
		        cat(anovaRemark)
		        cat("\n")
		        result[[z]]$anova.table <- AOV_print
		        
		        # compute Sum of Squares of GCA, SCA, REc, MAT, NORM
		        data.all.gca <- merge(meanData, as.data.frame(cbind(M3,GCA,SCA,REC)), by.x=c("newCodeP1","newCodeP2"), by.y=c("i","j"))
		        #data.all.gca <- data.all.gca[,-c(match("i", names(data.all.gca)),match("j", names(data.all.gca)))]
		        
		        #myformula4 <- paste(respvar[z],".mean ~ ", environment, " + . - ", p1, " - ", p2, " + (. - ", p1, " - ", p2, "):",environment, sep = "")
		        myformula4 <- paste(respvar[z],".mean ~ ", environment, " + . - newCodeP1 - newCodeP2 + (. - newCodeP1 - newCodeP2):", environment, sep = "")
            
		        #GCA
		        out.gca <- anova(aov(formula(myformula4), data = data.all.gca))
		        GCA.effect <- subset(out.gca, substr(rownames(out.gca),1,3)=="GCA")
		        GCA.df <- sum(GCA.effect$Df)
		        GCA.SS <- sum(GCA.effect$"Sum Sq")*r
		        GCA.MS <- GCA.SS/GCA.df
		        
		        envgca <- paste(environment, ":GCA",sep="")
		        GCAxE.effect <- subset(out.gca, substr(rownames(out.gca),1,(nchar(environment)+4))==envgca)
		        GCAxE.df <- sum(GCAxE.effect$Df)
		        GCAxE.SS <- sum(GCAxE.effect$"Sum Sq")*r
		        GCAxE.MS <- GCAxE.SS/GCAxE.df
		        
		        #SCA
		        SCA.effect <- subset(out.gca, substr(rownames(out.gca),1,3)=="SCA")
		        SCA.df <- sum(SCA.effect$Df)
		        SCA.SS <- sum(SCA.effect$"Sum Sq")*r
		        SCA.MS <- SCA.SS/SCA.df
		        
		        envsca <- paste(environment, ":SCA",sep="")
		        SCAxE.effect <- subset(out.gca, substr(rownames(out.gca),1,(nchar(environment)+4))==envsca)
		        SCAxE.df <- sum(SCAxE.effect$Df)
		        SCAxE.SS <- sum(SCAxE.effect$"Sum Sq")*r
		        SCAxE.MS <- SCAxE.SS/SCAxE.df
		        
		        #REC
		        REC.effect <- subset(out.gca, substr(rownames(out.gca),1,3)=="REC")
		        REC.df <- sum(REC.effect$Df)
		        #REC.SS <- sum(REC.effect$"Sum Sq")*r
		        REC.SS <- AOV["Crosses","SS"]-GCA.SS-SCA.SS
		        REC.MS <- REC.SS/REC.df
		        
		        envrec <- paste(environment, ":REC",sep="")
		        RECxE.effect <- subset(out.gca, substr(rownames(out.gca),1,(nchar(environment)+4))==envrec)
		        RECxE.df <- sum(RECxE.effect$Df)
		        #RECxE.SS <- sum(RECxE.effect$"Sum Sq")*r
		        RECxE.SS <-AOV[paste("Crosses x ",environment,sep=""),"SS"]-GCAxE.SS-SCAxE.SS
		        RECxE.MS <- RECxE.SS/RECxE.df
		        
		        # MAT
		        data.all.gca <- merge(meanData, as.data.frame(cbind(M3,GCA,SCA,MAT,NONM)), by.x=c("newCodeP1","newCodeP2"), by.y=c("i","j"))
		        #data.all.gca <- data.all.gca[,-c(match("i", names(data.all.gca)),match("j", names(data.all.gca)))]
		        
		        out.mat <- anova(aov(formula(myformula4), data = data.all.gca))
		        MAT.effect <- subset(out.mat, substr(rownames(out.mat),1,3)=="MAT")
		        MAT.df <- sum(MAT.effect$Df)
		        MAT.SS <- sum(MAT.effect$"Sum Sq")*r
		        MAT.MS <- MAT.SS/MAT.df
		        
		        envmat <- paste(environment, ":MAT",sep="")
		        MATxE.effect <- subset(out.mat, substr(rownames(out.mat),1,(nchar(environment)+4))==envmat)
		        MATxE.df <- sum(MATxE.effect$Df)
		        MATxE.SS <- sum(MATxE.effect$"Sum Sq")*r
		        MATxE.MS <- MATxE.SS/MATxE.df
		        
		        #NONM
		        NONM.effect <- subset(out.mat, substr(rownames(out.mat),1,4)=="NONM")
		        NONM.df <- sum(NONM.effect$Df)
		        #NONM.SS <- sum(NONM.effect$"Sum Sq")*r
		        NONM.SS <- REC.SS-MAT.SS
		        NONM.MS <- NONM.SS/NONM.df
		        
		        envnonm <- paste(environment, ":NONM",sep="")
		        NONMxE.effect <- subset(out.mat, substr(rownames(out.mat),1,(nchar(environment)+5))==envnonm)
		        NONMxE.df <- sum(NONMxE.effect$Df)
		        #NONMxE.SS <- sum(NONMxE.effect$"Sum Sq")*r
		        NONMxE.SS <- RECxE.SS-MATxE.SS
		        NONMxE.MS <- NONMxE.SS/NONMxE.df
		        
		        # --- construct ANOVA table --- #
		        DF=c(GCA.df,SCA.df,REC.df,MAT.df,NONM.df,GCAxE.df,SCAxE.df,RECxE.df,MATxE.df,NONMxE.df,EDF)
		        SS=c(GCA.SS,SCA.SS,REC.SS,MAT.SS,NONM.SS,GCAxE.SS,SCAxE.SS,RECxE.SS,MATxE.SS,NONMxE.SS,ESS)
		        MS<-SS/DF
		        F_value<-c(MS[1]/MS[6],MS[2]/MS[7],MS[3]/MS[8],MS[4]/MS[9],MS[5]/MS[10],MS[6]/MS[11],MS[7]/MS[11],MS[8]/MS[11],MS[9]/MS[11],MS[10]/MS[11],MS[11]/MS[11])
		        p_value<-pf(F_value,DF,c(DF[6],DF[7],DF[8],DF[9],DF[10],DF[11],DF[11],DF[11],DF[11],DF[11],DF[11]),lower.tail = F)
		        AOV2<-cbind(DF,SS, MS,F_value,p_value)
		        
		        row.names(AOV2)<-c("GCA", "SCA", "REC", "  MAT", "  NONM", "GCAxE", "SCAxE", "RECxE", "  MATxE", "  NONMxE", "Residuals")
		        AOV2[11,4:5]<-NA
		        AOV2<-as.data.frame(AOV2)
		        AOV2_print<-formatAnovaTable(AOV2)
		        
		        cat("\n\nANOVA TABLE:\n\n")
		        print(AOV2_print)
		        cat("-------\n")
		        cat(anovaRemark)
		        cat("\n")
		        result[[z]]$anova2.table <- AOV2_print
		        
		        # --- Estimates of GCA, SCA, and REC effects --- #
		        
		        myformula4<- paste(respvar[z], "  ~ newCodeP1 + newCodeP2", sep = "")
		        meandata <- summaryBy(formula(myformula4), data=dataForAnova)
		        
		        # serial to parallel of meandata
		        mdata <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
		        dim(mdata) <- c(p,p)
		        
		        for (I in 1:p)  {
		          for (J in 1:p)   {
		            if (I != J) mdata[I,J] <- meandata[(meandata[,match("newCodeP1", names(meandata))]==I & meandata[,match("newCodeP2", names(meandata))]==J),3]
		          } 
		        }
		        colnames(mdata) <- rownames(mdata) <- seq(1:p)
		        
		        # --- printing the matrix of means --- #
		        mdata2 <- format(round(mdata, 4), digits=4, nsmall=4)
		        
		        # recode back to the user's notation
		        rownames(mdata2) <- colnames(mdata2) <- codingGuide$levelsParents[match(colnames(mdata2), codingGuide$newCodeParents)]
		        colnames(mdata2) <- paste(p2,"=",colnames(mdata2), sep="")
		        rownames(mdata2) <- paste(p1,"=",rownames(mdata2), sep="")
            mdata3 <-noquote(format(gsub(" 0.0000", "", mdata2),justify="right"))
		        cat("\n\nMATRIX OF MEANS:\n\n")
		        print(mdata3)
		        cat("-------\n")
		        cat(anovaRemark)
		        cat("\n")
		        result[[z]]$means.matrix <- mdata3
		        
		        XI <- rowSums(mdata)
		        XJ <- colSums(mdata)		
		        SUMX <- sum(mdata)
		        
		        #----- Estimates of GCA Effects -----#
		        G_SCA <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
		        dim(G_SCA) <- c(p,p)
		        for (I in 1:p)  {
		          for (J in 1:p)   {
		            if (I == J) G_SCA[I,J] <- (1/(2*p*(p-2)))*(p*(XI[I]+XJ[J])-(2*SUMX))
		          } 
		        }
		        
		        #----- Estimates of SCA Effects -----#
		        for (I in 1:p)  {
		          for (J in 1:p)   {
		            if (I < J)   {
		              B1 <-(1/2)*(mdata[I,J]+mdata[J,I])
		              B2 <- (1/(2*(p-2)))* (XI[I] + XJ[I] + XI[J] + XJ[J])
		              B3 <- (1/((p-1)*(p-2)))*SUMX
		              G_SCA[I,J] <- B1 - B2 + B3
		            }
		          } 
		        }
		        
		        #----- Estimates of Reciprocal Effects -----#
		        for (I in 1:p)  {
		          for (J in 1:p)   {
		            if (I > J)   G_SCA[I,J] <- (1/2)*(mdata[J,I] - mdata[I,J])
		          } 
		        }
		        
		        colnames(G_SCA) <- rownames(G_SCA) <- seq(1:p)
		        G_SCA2 <-noquote(format(round(G_SCA,4), digits=4, nsmall=4))
		        
		        # recode back to the user's notation
		        rownames(G_SCA2) <- colnames(G_SCA2) <- codingGuide$levelsParents[match(colnames(G_SCA2), codingGuide$newCodeParents)]
		        colnames(G_SCA2) <- paste(p2,"=",colnames(G_SCA2), sep="")
		        rownames(G_SCA2) <- paste(p1,"=",rownames(G_SCA2), sep="")
		        cat("\n\nGENERAL COMBINING ABILITY EFFECTS (diagonal), SPECIFIC COMBINING ABILITY EFFECTS\n")
		        cat("(above diagonal) AND RECIPROCAL EFFECTS (below diagonal): (AVERAGED OVER ENVIRONMENTS)\n\n")
		        print(G_SCA2)
		        result[[z]]$gcasca.matrix <- G_SCA2
		        
		        # --- estimates of standard errors --- #
		        EMS<-ESS/EDF
		        MEPRIME <-EMS/r  
		        Ve <- MEPRIME
		        SE_GI <- sqrt((p-1)*Ve/(2*e*p*(p-2)))
		        LSD_GI <- NA
		        SE_SIJ <- sqrt((p-3)*Ve/(2*e*(p-1)))
		        LSD_SIJ <- NA
		        SE_RIJ <- sqrt(Ve/(2*e))
		        LSD_RIJ <- NA
		        
		        SE_Gdiff <- sqrt(Ve/(e*(p-2)))
		        LSD_Gdiff <- qt(.975,EDF)*SE_Gdiff
		        SE_SIJ_SIK <- sqrt((p-3)*Ve/(e*(p-2)))
		        LSD_SIJ_SIK <- qt(.975,EDF)*SE_SIJ_SIK
		        SE_SIJ_SKL <- sqrt((p-4)*Ve/(e*(p-2)))
		        LSD_SIJ_SKL <- qt(.975,EDF)*SE_SIJ_SKL
		        
		        STDERR <- round(rbind(SE_GI,SE_SIJ, SE_RIJ, SE_Gdiff, SE_SIJ_SIK, SE_SIJ_SKL), digits=4)
		        LSD <- round(rbind(LSD_GI,LSD_SIJ, LSD_RIJ, LSD_Gdiff, LSD_SIJ_SIK, LSD_SIJ_SKL),digits=4)
		        VAREST <- as.table(cbind(STDERR, LSD))
		        
		        rownames(VAREST) <- c('Gi', 'Sij', 'Rij', 'Gi-Gj', 'Sij-Sik', 'Sij-Skl')
		        colnames(VAREST) <- c("Std. Error", "LSD")
		        cat("\n\nTABLE OF STANDARD ERRORS AND LSDs:\n\n")
		        print(VAREST)
		        result[[z]]$stderror.table <- VAREST
		        
		        #--- Variance Component ---#
		        
		        Vg <- (GCA.MS-GCAxE.MS-SCA.MS+SCAxE.MS)/(r*e*(p-2)) #removed 2 in the formula (2*r*e*(p-2))
		        Vs <- (SCA.MS-SCAxE.MS)/(2*r*e)
		        Vr <- (REC.MS-RECxE.MS)/(2*r*e)
		        Vge <- (GCAxE.MS-SCAxE.MS)/(2*r*(p-2))
		        Vse <- (SCAxE.MS-EMS)/(2*r)
		        Vre <- (RECxE.MS-EMS)/(2*r)
		        
		        if (Vg<0) Vg <- 0
		        if (Vs<0) Vs <- 0
		        if (Vr<0) Vr <- 0
		        if (Vge<0) Vge <- 0
		        if (Vse<0) Vse <- 0
		        if (Vre<0) Vre <- 0
		        if (EMS < 0) EMS <- 0
            
		        VC <- format(round(rbind(Vg, Vs, Vr, Vge, Vse, Vre, EMS), digits=4), digits=4, nsmall=4, scientific=FALSE)
		        colnames(VC) <- c("Estimate")
		        rownames(VC) <- c("GCA", "SCA", "REC","GCAxE", "SCAxE", "RECxE", "Residuals")
		        TABLE <- as.table(VC)
		        cat("\n\nESTIMATES OF VARIANCE COMPONENTS:\n\n")
		        print(TABLE)
		        result[[z]]$var.components <-TABLE
		        
		        
		        #---- Genetic Variance components ----#
		        if (cross) {F<-0
		        } else {F<-1}
		        
		        VA <- (4/(1+F))*Vg
		        VD <- (4/(1+F)^2)*Vs
		        VAE <- (4/(1+F))*Vge
		        VDE <- (4/(1+F)^2)*Vse
		        VE <- EMS
            
		        if (VA < 0 || VA < 1e-10) VA <- 0
		        if (VD < 0 || VD < 1e-10) VD <- 0
		        
		        VP <- VA + VD + VAE + VDE+ VE
		        h2B <- (VA + VD) / VP               # individual based
		        h2N <- VA / VP                 	    # individual based
		        Dominance.ratio <- sqrt(2*VD/VA)   
		        
		        VA_p<-formatNumericValue(VA)
		        VAE_p<-formatNumericValue(VAE)
		        VD_p<-formatNumericValue(VD)
		        VDE_p<-formatNumericValue(VDE)
		        h2N_p<-formatNumericValue(h2N)
		        h2B_p<-formatNumericValue(h2B)
		        Dominance.ratio_p<-formatNumericValue(Dominance.ratio)
		        
		        Estimate <- rbind(VA_p, VAE_p,  VD_p, VDE_p, h2N_p, h2B_p, Dominance.ratio_p)
		        with.colheader<-format(rbind("Estimate", Estimate), justify="right")
		        colnames(with.colheader) <- c("")
		        rownames(with.colheader) <- c("", "VA", "VAxE", "VD", "VDxE", "h2-narrow sense", "H2-broad sense", "Dominance Ratio")
		        TABLE2 <- as.table(with.colheader)
		        cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:\n")
		        print(TABLE2)
		        result[[z]]$genvar.components <-TABLE2 
		      } else {
		        cat("\n",anovaRemark)
		      }
		    } else {anovaRemark <- "ERROR: Too many missing values. Cannot perform the analysis." 
		            cat(" ",anovaRemark)
		            result[[z]]$canEstimateMissing<-"FALSE"
		    }
		    cat("\n")
		  } ## end of else (if (nBalance<nrow(temp.data)))
		} ## end of else (if ((nrow(temp.data)/nrow(data)) < 0.80))
	} ## end of loop (z)
  cat("\n==============================================================\n")
  #detach("package:doBy")
  return(list(output = result))
  
}

