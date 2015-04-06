#----------------------------------------------------------------
# This function performs analysis for Diallel I experiments using Griffing Method          
# (Parents + F1's + reciprocals)
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
# Script Created by: Violeta Bartolome   07.2011
# Script Modified by: Nellwyn Sales and Guoyou Ye 
#----------------------------------------------------------------

diallel1Test <-function(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, p1, p2, rep=NULL, block=NULL, row=NULL, column=NULL, cross=TRUE, individual=NULL, environment=NULL, alpha=0.05) UseMethod("diallel1Test")

diallel1Test.default <-function(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, p1, p2, rep=NULL, block=NULL, row=NULL, column=NULL, cross=TRUE, individual=NULL, environment=NULL, alpha=0.05) {
	
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
	if (!is.null(environment)) {environment <-trimStrings(environment) }
	alpha <- trimStrings(alpha)
	
	# --- create titles --- #
	if (design == "CRD") { designName<-"CRD"}
	if (design == "RCB") { designName<-"RCB"}
	if (design == "Alpha") { designName<-"ALPHA-LATTICE"}
	if (design == "RowColumn") { designName<-"ROW-COLUMN"}
  
	if (cross) {parentsType<-"CROSS"
	} else {parentsType<-"SELF"}
	cat("\nDIALLEL ANALYSIS: GRIFFING METHOD I IN ",designName, " (", parentsType, ")\n", sep="")
	
	# --- get number of environment levels --- #
	if (!is.null(environment)) {
	  data[,match(environment, names(data))] <- factor(trimStrings(data[,match(environment, names(data))]))
	  envNumLevels<-nlevels(data[,match(environment, names(data))])
	} else {envNumLevels<-1}
	
	result <- list()
	for (i in (1:length(respvar))) {
	  result[[i]] <- list()
	  cat("\n-----------------------------")
	  cat("\nRESPONSE VARIABLE: ", respvar[i], "\n", sep="")
	  cat("-----------------------------\n")
	  for (j in (1:envNumLevels)) {
	    result[[i]]$site[[j]] <- list()
	    if (!is.null(environment)) {
	      crdVars<-c(respvar[i], p1, p2, environment)
	      rcbVars<-c(respvar[i], p1, p2, block, environment)
	      alphaVars<-c(respvar[i], p1, p2, rep, block, environment)
	      rowcolVars<-c(respvar[i], p1, p2, rep, row, column, environment)
	    } else {
	      crdVars<-c(respvar[i], p1, p2)
	      rcbVars<-c(respvar[i], p1, p2, block)
	      alphaVars<-c(respvar[i], p1, p2, rep, block)
	      rowcolVars<-c(respvar[i], p1, p2, rep, row, column)
	    }
	    if (design == "CRD") {temp.data <- data[sort(match(crdVars, names(data)))]}
	    if (design == "RCB") {temp.data <- data[sort(match(rcbVars, names(data)))]}
	    if (design == "Alpha") {temp.data <- data[sort(match(alphaVars, names(data)))]}
	    if (design == "RowColumn") {temp.data <- data[sort(match(rowcolVars, names(data)))]}
      
	    if (!is.null(environment)) {
	      cat("-----------------------------")
	      cat("\nANALYSIS FOR: ",environment, " = " ,levels(temp.data[,match(environment, names(temp.data))])[j],"\n", sep="")
	      cat("-----------------------------\n")
        
	      result[[i]]$site[[j]]$env <- levels(temp.data[,match(environment, names(temp.data))])[j]
        
	      temp.data <- temp.data[temp.data[,match(environment, names(temp.data))] == levels(temp.data[,match(environment, names(temp.data))])[j],]
        temp.data[,match(environment, names(temp.data))] <- factor(trimStrings(temp.data[,match(environment, names(temp.data))]))
	    }
	    
	    # --- define factors and number of levels --- #
	    obsread<-nrow(temp.data)
			temp.data[,match(p1, names(temp.data))] <- factor(trimStrings(temp.data[,match(p1, names(temp.data))]))
			temp.data[,match(p2, names(temp.data))] <- factor(trimStrings(temp.data[,match(p2, names(temp.data))]))
	    p <- length(unique(c(levels(temp.data[,match(p1, names(temp.data))]), levels(temp.data[,match(p2, names(temp.data))]))))
      
	    # --- create new column containing treatment combinations --- #
	    temp.data$cross<-factor(paste(temp.data[,p1], ":", temp.data[,p2], sep=""))
	    temp.data<-temp.data[order(temp.data$cross),]
	    
	    # --- compute harmonic mean that will be used later in the estimation of genetic variances --- #
	    lengthPerCross<-tapply(temp.data[,respvar[i]], temp.data$cross, length)
	    repHarmonicMean<-1/mean((1/lengthPerCross), na.rm=TRUE)
	    
	    if (design == "CRD") {
	      # --- add column Rep --- #
	      temp.data<-data.frame(temp.data, Rep=sequence(lengthPerCross))
	      
	      nlevelsRep<-max(lengthPerCross, na.rm=TRUE)
	    }
	    if (design == "RCB") {
	      temp.data[,match(block, names(temp.data))] <- factor(trimStrings(temp.data[,match(block, names(temp.data))]))
	      nlevelsRep <- nlevels(temp.data[,match(block, names(temp.data))])
	    }
	    if (design == "Alpha") {
	      temp.data[,match(rep, names(temp.data))] <- factor(trimStrings(temp.data[,match(rep, names(temp.data))]))
	      temp.data[,match(block, names(temp.data))] <- factor(trimStrings(temp.data[,match(block, names(temp.data))]))
	      nlevelsRep <- nlevels(temp.data[,match(rep, names(temp.data))])
	      if (!is.null(environment)) {
	        blockSizeFrame<-as.data.frame.table(tapply(temp.data[,respvar[i]], temp.data[,c(environment, rep, block)], length))
	      } else {
	        blockSizeFrame<-as.data.frame.table(tapply(temp.data[,respvar[i]], temp.data[,c(rep, block)], length))
	      }
	      blockSize<-max(blockSizeFrame$Freq, na.rm=TRUE)
	      nBlocksWithinRep<-(p*p)/blockSize
	    }
	    if (design == "RowColumn") {
	      temp.data[,match(rep, names(temp.data))] <- factor(trimStrings(temp.data[,match(rep, names(temp.data))]))
	      temp.data[,match(row, names(temp.data))] <- factor(trimStrings(temp.data[,match(row, names(temp.data))]))
	      temp.data[,match(column, names(temp.data))] <- factor(trimStrings(temp.data[,match(column, names(temp.data))]))
	      nlevelsRep <- nlevels(temp.data[,match(rep, names(temp.data))])
	      
	      if (!is.null(environment)) {
	        rowWithinRepFrame<-as.data.frame.table(tapply(temp.data[,respvar[i]], temp.data[,c(environment, rep, row)], length))
	      } else {
	        rowWithinRepFrame<-as.data.frame.table(tapply(temp.data[,respvar[i]], temp.data[,c(rep, row)], length))
	      }
	      rowWithinRep<-max(rowWithinRepFrame$Freq, na.rm=TRUE)
	      
	      if (!is.null(environment)) {
	        columnWithinRepFrame<-as.data.frame.table(tapply(temp.data[,respvar[i]], temp.data[,c(environment, rep, column)], length))
	      } else {
	        columnWithinRepFrame<-as.data.frame.table(tapply(temp.data[,respvar[i]], temp.data[,c(rep, column)], length))
	      }
	      columnWithinRep<-max(columnWithinRepFrame$Freq, na.rm=TRUE)
	    }
	    temp.data<-temp.data[-c(match("cross", names(temp.data)))]
	    nBalance=p*p*nlevelsRep
      
      # --- check if max of lengthPerCross is equal to nlevelsRep --- #
	    result[[i]]$site[[j]]$maxLengthPerCross <- max(lengthPerCross, na.rm=TRUE)
	    result[[i]]$site[[j]]$nlevelsRep <- nlevelsRep
      
      if (max(lengthPerCross, na.rm=TRUE) > nlevelsRep) {
        if (design == "RCB") {
          blockLabelError <- paste("The number of levels of the blocking factor is", nlevelsRep, "but at least one treatment combination is replicated", max(lengthPerCross, na.rm=TRUE), "times.")
          blockLabelError2 <- paste("Please check if the column for block is properly labeled.")
        } else {
          blockLabelError <- paste("The number of levels of the replicate factor is", nlevelsRep, "but at least one treatment combination is replicated", max(lengthPerCross, na.rm=TRUE), "times.")
          blockLabelError2 <- paste("Please check if the column for replicate is properly labeled.")
        }
        result[[i]]$site[[j]]$blockLabelError <- blockLabelError
        result[[i]]$site[[j]]$blockLabelError2 <- blockLabelError2
        cat("\n ERROR:", blockLabelError)
        cat("\n       ", blockLabelError2, "\n\n\n")
        break
      }
      
      if (design == "Alpha" || design == "RowColumn") {
        temp.data.withNA <- temp.data
      }
	    
	    # --- remove rows with missing observations --- #
	    temp.data <- temp.data[(is.na(temp.data[,match(respvar[i], names(temp.data))]) == FALSE),]
	    obsused<-nrow(temp.data)
      
	    responseRate<-(obsused/nBalance)
	    
	    if (responseRate < 0.80) {
	      cat("\nToo many missing observations. Cannot proceed with the analysis.\n\n\n")
	      result[[i]]$site[[j]]$tooManyNAWarning <- "YES"
	      next
	    } else {
	      result[[i]]$site[[j]]$tooManyNAWarning <- "NO"
        
	      # --- data summary --- #
	      funcTrialSum <- class.information2(names(temp.data),temp.data)
	      cat("\nDATA SUMMARY: ","\n\n", sep="")
	      print(funcTrialSum)
	      cat("\nNumber of observations read: ",obsread, sep="")
	      cat("\nNumber of observations used: ",obsused, sep="")
	      missingObs<-nBalance-obsused
        
	      result[[i]]$site[[j]]$funcTrialSum <- funcTrialSum
	      result[[i]]$site[[j]]$obsread <- obsread
	      result[[i]]$site[[j]]$obsused <- obsused
	      
	      # --- call recodeDiallelData to recode p1 and p2 to standard notation and generate balancedData --- #
	      if (design == "CRD") {
	        outRecode<-recodeDiallelData(design="diallel1", data=temp.data, p1=p1, p2=p2, rep="Rep", respvar=respvar[i])
	        balancedData<-outRecode$tempData
	        codingGuide<-outRecode$newCoding
	      }
	      if (design == "RCB") {
	        outRecode<-recodeDiallelData(design="diallel1", data=temp.data, p1=p1, p2=p2, rep=block, respvar=respvar[i])
	        balancedData<-outRecode$tempData
	        codingGuide<-outRecode$newCoding
	      }
	      if (design == "Alpha") {
          if (nrow(temp.data.withNA)==nBalance) {
            outRecode<-recodeDiallelData(design="diallel1", data=temp.data.withNA, p1=p1, p2=p2, rep=rep, block=block, respvar=respvar[i])
            balancedData<-outRecode$tempData
            codingGuide<-outRecode$newCoding
          } else stop("The dataset has missing row(s). PBTools cannot supply the value(s) of block factor for the missing row(s).")
	      }
	      if (design == "RowColumn") {
	        if (nrow(temp.data.withNA)==nBalance) {
	          outRecode<-recodeDiallelData(design="diallel1", data=temp.data.withNA, p1=p1, p2=p2, rep=rep, row=row, column=column, respvar=respvar[i])
	          balancedData<-outRecode$tempData
	          codingGuide<-outRecode$newCoding
	        } else stop("The dataset has missing row(s). PBTools cannot supply the values of row and column factors for the missing row(s).")
	      }
        
	      if (nBalance<nrow(temp.data)) {
	        cat("\n\n***\nERROR: The number of observations read in the data exceeds the size of a balanced data.\n       Please check if the column for block/replicate is properly labeled.\n***\n\n")
	        result[[i]]$site[[j]]$exceededWarning <- "YES"
	      } else {
	        result[[i]]$site[[j]]$exceededWarning <- "NO"
                    
	        # --- ANOVA for Diallel Method 1 experiment --- #
	        estimatedMissing <- FALSE
          
	        result[[i]]$site[[j]]$responseRate <- (nrow(temp.data)/nBalance)
                    
	        if ((nrow(temp.data)/nBalance) >= 0.90) {
	          if (nrow(temp.data) == nBalance) {
	            anovaRemark1 <- "Raw dataset is balanced."
	            dataForAnova<-balancedData   
	          } else {
	            if (design == "CRD") {
	              #dataForAnova<-estimateMissingData(design="CRD", data=balancedData, respvar[i], "newCodeP1", "newCodeP2", "Rep")
	              anovaRemark1 <- "Raw dataset is unbalanced."
	            }
	            if (design == "RCB") {
	              dataForAnova<-estimateMissingData(design="RCB", data=balancedData, respvar[i], "newCodeP1", "newCodeP2", block)
	              anovaRemark1 <- "Raw data and estimates of the missing values are used."
	            }
	            if (design == "Alpha") {
	              dataForAnova<-estimateNA(design="Alpha", fullData=balancedData, respvar[i], "newCodeP1", "newCodeP2", rep, block, row=NULL, column=NULL)
	              anovaRemark1 <- "Raw data and estimates of the missing values are used."
	            }
	            if (design == "RowColumn") {
	              dataForAnova<-estimateNA(design="RowColumn", fullData=balancedData, respvar[i], "newCodeP1", "newCodeP2", rep, block=NULL, row=row, column=column)
	              anovaRemark1 <- "Raw data and estimates of the missing values are used."
	            }
	            estimatedMissing <- TRUE
	          }
	        } else {
	          if (design == "CRD") {
	            anovaRemark1 <- "Raw dataset is unbalanced."
	          } else {
	            anovaRemark1 <- "ERROR: Too many missing values. Cannot perform the analysis."
	          }
	        }
	        result[[i]]$site[[j]]$anovaRemark1 <- anovaRemark1
         
	        if (design == "CRD" || design == "RCB") {
	          cat("\n\nANOVA TABLE FOR THE EXPERIMENT: \n", sep="")
	          if (design == "CRD") {
	            myformula0 <- paste(respvar[i], " ~ + newCodeP1*newCodeP2", sep = "")
	            anovaFixed<-summary(aov(formula(myformula0), data=balancedData))
              
	            rownames(anovaFixed[[1]])[match("newCodeP1:newCodeP2", trimStrings(rownames(anovaFixed[[1]])))]<-paste(p1,":",p2, sep="")
	            rownames(anovaFixed[[1]])[match("newCodeP1", trimStrings(rownames(anovaFixed[[1]])))]<-p1
	            rownames(anovaFixed[[1]])[match("newCodeP2", trimStrings(rownames(anovaFixed[[1]])))]<-p2
	            rownames(anovaFixed[[1]]) <- trimStrings(rownames(anovaFixed[[1]]))
              
	            anovaFormat<-formatAnovaTable(anovaFixed)
	            print(anovaFormat)
	            cat("-------\n")
	            cat(paste("REMARK: ",anovaRemark1, sep=""))
	            result[[i]]$site[[j]]$diallel1.anova <- anovaFormat
              
	          }
            
	          if (design == "RCB") {
	            if ((nrow(temp.data)/nBalance) >= 0.90) {
	              myformula0 <- paste(respvar[i], " ~ ", block," + newCodeP1*newCodeP2", sep = "")
	              model0 <- aov(formula(myformula0), data = dataForAnova)
	              anovaFixed<-summary(model0)
                
	              rownames(anovaFixed[[1]])[match("newCodeP1:newCodeP2", trimStrings(rownames(anovaFixed[[1]])))]<-paste(p1,":",p2, sep="")
	              rownames(anovaFixed[[1]])[match("newCodeP1", trimStrings(rownames(anovaFixed[[1]])))]<-p1
	              rownames(anovaFixed[[1]])[match("newCodeP2", trimStrings(rownames(anovaFixed[[1]])))]<-p2
	              rownames(anovaFixed[[1]]) <- trimStrings(rownames(anovaFixed[[1]]))
	              
	              #rerun aov using temp.data to get the original df's
	              anovaFixed.temp<-summary(aov(formula(myformula0), data=balancedData))
	              anovaFixed.temp[[1]]$"Df"[length(anovaFixed.temp[[1]]$"Df")]<-anovaFixed[[1]]$"Df"[length(anovaFixed[[1]]$"Df")]-missingObs
	              anovaFormat<-adjustAnovaDf(anovaFixed, anovaFixed.temp[[1]]$"Df")
	              anovaFormat<-formatAnovaTable(anovaFormat)
                
	              print(anovaFormat)
	              cat("-------\n")
	              cat(paste("REMARK: ",anovaRemark1, sep=""))
	              result[[i]]$site[[j]]$diallel1.anova <- anovaFormat
	            } else {
	              cat(anovaRemark1)
	            }
	          }
	        }
	        
	        # --- testing for genotypic effect --- #
	        if (design == "CRD") {
	          Crosses <- balancedData[,match("newCodeP1", names(balancedData))]:balancedData[,match("newCodeP2", names(balancedData))]
	        } else {
	          if ((nrow(temp.data)/nBalance) >= 0.90) {
	            Crosses <- dataForAnova[,match("newCodeP1", names(dataForAnova))]:dataForAnova[,match("newCodeP2", names(dataForAnova))]
	          }
	        }
	        
	        pValue <- 0
	        cat("\n\n\nTESTING FOR THE SIGNIFICANCE OF CROSS EFFECT: (Crosses = ", p1,":", p2, ")\n", sep="")
	        
	        # --- assign value for r --- #
	        r<-nlevelsRep
	        
	        if (design == "CRD") {
	          myformula1 <- paste(respvar[i], " ~ Crosses",sep = "") 
	          model1 <- aov(formula(myformula1), data=balancedData)
	          anova.table<-summary(model1)
	          pValue<-anova.table[[1]]$"Pr(>F)"[[1]]
	          
	          # print anova table
	          anova_print <- formatAnovaTable(anova.table)
	          cat("\nFormula: ", myformula1, "\n\n", sep="")
	          
	          result[[i]]$site[[j]]$genoEffect.anova <-anova_print
	          result[[i]]$site[[j]]$pValue <- pValue
            
	          #get MSE
	          EMS <- as.numeric(anova.table[[1]]$"Mean Sq"[2])
	          EDF <- as.numeric(anova.table[[1]]$"Df"[2])
            
	          print(anova_print)
	          result[[i]]$site[[j]]$formula1 <- myformula1
	        }
	        if (design == "RCB") {
	          if ((nrow(temp.data)/nBalance) >= 0.90) {
	            myformula1 <- paste(respvar[i], " ~ Crosses + (1|", block, ")", sep = "") 
	            myformula2 <- paste(respvar[i], " ~ (1|", block, ")", sep = "")
	            library(lme4)
	            model1 <- lmer(formula(myformula1), data=dataForAnova) 
	            model2 <- lmer(formula(myformula2), data=dataForAnova)  
	            
	            anova_print<-modelComparisonTable(model1, model2)
	            pValue<-as.numeric(toString(anova_print[2,"Pr(>Chisq)"]))
	            
	            cat("\nFormula for Model 1: ", myformula1, sep="")
	            cat("\nFormula for Model 2: ", myformula2,"\n\n", sep="")
	            
	            # get MSE
	            varcomp <- summary(model1)@REmat
	            EMS <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
	            EDF <- (((p*p)-1)*(r-1)) - missingObs
	            
	            print(anova_print)
	            result[[i]]$site[[j]]$formula1 <- myformula1
	          } else {
	            cat(anovaRemark1)
	          }
	        }
	        if (design == "Alpha") {
	          if ((nrow(temp.data)/nBalance) >= 0.90) {
	            myformula1 <- paste(respvar[i], " ~ Crosses + (1|", rep, ") + (1|", block, ":", rep, ")", sep = "") 
	            myformula2 <- paste(respvar[i], " ~ (1|", rep, ") + (1|", block, ":", rep, ")", sep = "")
	            library(lme4)
	            model1 <- lmer(formula(myformula1), data=dataForAnova) 
	            model2 <- lmer(formula(myformula2), data=dataForAnova)  
	            
	            anova_print<-modelComparisonTable(model1, model2)
	            pValue<-as.numeric(toString(anova_print[2,"Pr(>Chisq)"]))
	            
	            cat("\nFormula for Model 1: ", myformula1, sep="")
	            cat("\nFormula for Model 2: ", myformula2,"\n\n", sep="")
	            
	            #get MSE
	            varcomp <- summary(model1)@REmat
	            EMS <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
	            EDF <- (((r*p*p)-1)-(r-1)-(r*(nBlocksWithinRep-1))-((p*p)-1)) - missingObs
	            
	            print(anova_print)
	            result[[i]]$site[[j]]$formula1 <- myformula1
	          } else {
	            cat(anovaRemark1)
	          }
	        }
	        if (design == "RowColumn") {
	          if ((nrow(temp.data)/nBalance) >= 0.90) {
	            myformula1 <- paste(respvar[i], " ~ Crosses + (1|", rep, ") + (1|", row, ":", rep, ") + (1|", column, ":", rep, ")", sep = "") 
	            myformula2 <- paste(respvar[i], " ~ (1|", rep, ") + (1|", row, ":", rep, ") + (1|", column, ":", rep, ")", sep = "")
	            library(lme4)
	            model1 <- lmer(formula(myformula1), data=dataForAnova) 
	            model2 <- lmer(formula(myformula2), data=dataForAnova)  
	            
	            anova_print<-modelComparisonTable(model1, model2)
	            pValue<-as.numeric(toString(anova_print[2,"Pr(>Chisq)"]))
	            
	            cat("\nFormula for Model 1: ", myformula1, sep="")
	            cat("\nFormula for Model 2: ", myformula2,"\n\n", sep="")
	            
	            #get MSE
	            varcomp <- summary(model1)@REmat
	            EMS <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
	            numberTrt<-p*p
	            EDF <- (((numberTrt*r)-1)-(numberTrt-1)-(r-1)-((rowWithinRep-1)*r)-((columnWithinRep-1)*r)) - missingObs
	            
	            print(anova_print)
	            result[[i]]$site[[j]]$formula1 <- myformula1
	          } else {
	            cat(anovaRemark1)
	          }
	        }
	        
          
	        if (design == "CRD") {
	          result[[i]]$site[[j]]$data <- cbind(balancedData, Crosses)
	        } else {
	          if ((nrow(temp.data)/nBalance) >= 0.90) {
	            result[[i]]$site[[j]]$data <- cbind(dataForAnova, Crosses)
	          } 
	        }
	        
	        # --- mean data of full diallel --- #
	        myformula2<- paste(respvar[i], " ~ newCodeP1 + newCodeP2", sep = "")
          
	        if (design == "CRD") {
	          balancedDataNoNa <- balancedData[(is.na(balancedData[,respvar[i]]) == FALSE),]
	          meandata <- summaryBy(formula(myformula2), data=balancedDataNoNa)
	          
	          # --- check if there is missing cross --- #
	          if (nrow(meandata)<(p*p)) {
	            meandata <- summaryBy(formula(myformula2), data=balancedData, na.rm=TRUE)
	            meandata[,3][is.nan(meandata[,3])]<-NA
	            meansComplete <- FALSE
	          } else {
	            meansComplete <- TRUE
	          }
	        } else {
	          if ((nrow(temp.data)/nBalance) >= 0.90) {
	            meandata <- summaryBy(formula(myformula2), data=dataForAnova)
	            meansComplete <- TRUE
	          } else {
	            meansComplete <- FALSE
	          }
	        }
	        
	        result[[i]]$site[[j]]$anovaRemark2 <- anovaRemark1
	        result[[i]]$site[[j]]$meansComplete <-toString(meansComplete)
          
	        if (design == "CRD" || meansComplete) {
	          # --- serial to parallel of meandata --- #
	          mdata <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
	          dim(mdata) <- c(p,p)
	          
	          for (I in 1:p)  {
	            for (J in 1:p)   {
	              mdata[I,J] <- meandata[(meandata[,"newCodeP1"]==I & meandata[,"newCodeP2"]==J),3]
	            } 
	          }
	          colnames(mdata) <- rownames(mdata) <- seq(1:p)
	          
	          # --- printing the matrix of means --- #
	          mdata2 <- format(round(mdata,4), digits=4, nsmall=4)
	          
	          # recode back to the user's notation
	          rownames(mdata2) <- colnames(mdata2) <- codingGuide$levelsParents[match(colnames(mdata2), codingGuide$newCodeParents)]
	          colnames(mdata2) <- paste(p2,"=",colnames(mdata2), sep="")
	          rownames(mdata2) <- paste(p1,"=",rownames(mdata2), sep="")
	          mdata_print <-noquote(format(gsub(" 0.0000", "", mdata2),justify="right"))
	          cat("\n\nMATRIX OF MEANS:\n\n")
	          print(mdata_print)
	          result[[i]]$site[[j]]$means.matrix <-mdata_print
	          #result[[i]]$site[[j]]$estimatedMissing <-toString(estimatedMissing)
	          
	          cat("-------\n")
	          cat(paste("REMARK: ",anovaRemark1, sep=""))
	        }
          
          if (meansComplete) {
            # --- check if genotypic effect is significant. If significant, proceed to diallel analysis --- #
            alpha <- as.numeric(alpha)
            #if (pValue < alpha) {
            cat("\n\nANALYSIS OF VARIANCE:\n")
            if (design == "CRD" || (nrow(temp.data)/nBalance) >= 0.90) {
              
              XI <- rowSums(mdata)
              XJ <- colSums(mdata)
              SUMX <- sum(mdata)
              
              MEPRIME <-EMS/r            #  MSE PRIME (ME')  
              
              # --- GENERAL COMBINING ABILITY (GCA) SUM OF SQUARES ---- #
              SG <- round(((1/(2*p))*sum((XI+XJ)^2)) - ((2/(p^2))*(SUMX^2)), 6)
              
              # --- SPECIFIC COMBINING ABILITY (SCA) SUM OF SQUARES --- #
              B1 <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
              dim(B1) <- c(p,p)
              
              for (I in 1:p)  {
                for (J in 1:p)   {
                  B1[I,J] <- mdata[I,J]*(mdata[I,J] + mdata[J,I])
                } 
              }
              SS <- round(((1/2)*sum(B1))-((1/(2*p))*sum((XI+XJ)^2))+((1/(p^2))*(SUMX^2)),6)
              
              # --- RECIPROCAL SUM OF SQUARES ---- #
              B2 <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
              dim(B2) <- c(p,p)
              
              for (I in 1:p)  {
                for (J in 1:p)   {
                  if (I>J) B2[I,J] <- mdata[I,J] - mdata[J,I]
                } 
              }
              
              SR <- round((1/2)*sum(B2^2),6)
              
              # --- ERROR PRIME SUM OF SQUARES --- #
              SE <- round(MEPRIME * EDF,6)
              
              # --- COMPUTATION OF MEAN SQUARE AND F-VALUES ---- #
              DG <- p-1
              DS <- p*(p-1)/2
              DR <- DS
              DE <- EDF
              
              MG <- SG/DG
              MS <- SS/DS
              MR <- SR/DR
              
              # --- Error MS for GCA (MS*) --- #
              
              A <- p*(p-1)/(p^2-p+1)
              MS.star <- (1-A) * MEPRIME + A*MS
              c.prime <- p^2 - p + 1
              k <- ((p^2)/(2*c.prime))*((MS-MEPRIME)/MEPRIME) 
              #df.star <- ((EDF*p^2)*(p-1)*(p+(2*(p-1)*k))^2)/(((p^2)*(p-1)*((1-A)^2))+(2*EDF*A^2*(p^2+2*c.prime*k)))
              df.star <- (EDF*(p^2)*(p-1)*((p+(2*(p-1)*k))^2))/(((p^2)*(p-1)*((1-A)^2))+(2*EDF*(A^2)*((p^2)+(2*c.prime*k))))
              FG <- MG/MS.star
              FS <- MS/MEPRIME
              FR <- MR/MEPRIME
              FE <- NA
              
              PG <- 1-pf(FG, DG, df.star)
              PS <- 1-pf(FS, DS, EDF)
              PR <- 1-pf(FR, DR, EDF)
              PE <- NA
              
              # --- printing the anova table --- #
              DF <- c(DG,DS,DR,DE)
              SSq <- c(SG, SS, SR, SE)
              MSq <- c(MG, MS, MR, MEPRIME)
              F_value<-c(FG, FS, FR, FE)
              p_value <- c(PG,PS,PR,PE)
              AOV<-cbind(DF,SSq, MSq,F_value,p_value)
              row.names(AOV)<-c("GCA", "SCA", "Reciprocal", "Error")
              AOV<-as.data.frame(AOV)
              AOV_print<-formatAnovaTable(AOV)
              
              print(AOV_print)
              cat(" -----\n")
              cat(" NOTE: MS* = ",round(MS.star,digits=6), "   Error used for GCA MS with df = ",round(df.star, digits=1), "\n", sep="")
              
              result[[i]]$site[[j]]$gcasca.anova <-AOV_print  
              result[[i]]$site[[j]]$gcasca.anovaRemark <- paste("MS* = ",round(MS.star,digits=6), "   Error used for GCA MS with df = ",round(df.star, digits=1), "\n", sep="")
              
              #--- Estimation of variance components ---#
              
              Ve <- MEPRIME
              SEe <- sqrt((2/EDF)*(MEPRIME^2))
              
              Vr <- (1/2)*(MR -MEPRIME)
              if (Vr < 0) Vr <- 0
              SEr <- sqrt((1/(p*(p-1)))*MR^2 + (1/(2*EDF))*MEPRIME^2)
              
              Vs <- (p^2/(2*(p^2-p+1)))*(MS-MEPRIME)
              if (Vs < 0) Vs <- 0
              SEs <- sqrt((p^3/((p-1)*c.prime^2))*(MS^2) + ((p^4)/(2*(c.prime^2)*EDF))*(MEPRIME^2))
              
              Vg <- (1/(2*p))*(MG - (MEPRIME + p*(p-1)*MS)/(p^2 - p +1)) 
              if (Vg < 0) Vg <- 0 
              SEg <- sqrt((1/(2*p^2*(p-1)))*MG^2 + ((p-1)/(p*c.prime^2))*MS^2 + (1/(2*p^2*c.prime^2*EDF))*MEPRIME^2)
              
              VC <- round(rbind(Vg, Vs, Vr, Ve), digits=4)
              StdErr <- round(rbind(SEg, SEs, SEr, SEe), digits=4)
              TABLE <- (cbind(VC, StdErr))
              colnames(TABLE) <- c("Estimate", "Std. Error")
              rownames(TABLE) <- c("GCA", "SCA", "Reciprocal", "Error")
              TABLE_print<-as.table(TABLE)
              cat("\n\nESTIMATES OF VARIANCE COMPONENTS:\n\n")
              print(TABLE_print)
              result[[i]]$site[[j]]$var.components <-TABLE_print
              
              #---- Genetic Variance components ----#
              if (cross) {F<-0
              } else {F<-1}
              
              VA <- (4/(1+F))*Vg
              VVA <- 4*SEg^2
              VD <- (4/(1+F)^2)*Vs
              VVD <- SEs^2
              if (VA < 0 || VA < 1e-10) VA <- 0
              if (VD < 0 || VD < 1e-10) VD <- 0
              VE <- Ve
              VP <- VA + VD + VE
              h2B <- (VA + VD) / VP                 # individual based
              # Vh2B <- (2*(1-h2B)^2*(1+(r-1)*h2B)^2)/(r*(r-1)*(p^2-1))
              
              h2N <- VA / VP                 # individual based
              # Vh2N <- (2*(1-h2N)^2*(1+(r-1)*h2N)^2)/(r*(r-1)*(p^2-1))
              Dominance.ratio <- sqrt(2*VD/VA) 
              
              # --- format values to print in the table --- #
              VA_p<-formatNumericValue(VA)
              VD_p<-formatNumericValue(VD)
              h2N_p<-formatNumericValue(h2N)
              h2B_p<-formatNumericValue(h2B)
              Dominance.ratio_p<-formatNumericValue(Dominance.ratio)
              Estimate <- rbind(VA_p, VD_p, h2N_p, h2B_p, Dominance.ratio_p)
              
              with.colheader<-format(rbind("Estimate", Estimate), justify="right")
              colnames(with.colheader) <- c("")
              rownames(with.colheader) <- c("", "VA", "VD", "h2-narrow sense", "H2-broad sense", "Dominance Ratio")
              TABLE2 <- as.table(with.colheader)
              cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:\n")
              print(TABLE2)
              result[[i]]$site[[j]]$genvar.components <-TABLE2 
              
              #----  ESTIMATES OF GCA EFFECTS  ----#
              G_SCA <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
              dim(G_SCA) <- c(p,p)
              
              for (I in 1:p)  {
                for (J in 1:p)   {
                  if (I == J) G_SCA[I,J] <- ((1/(2*p)) * (XI[I] + XJ[J])) - (1/(p^2)) * SUMX
                } 
              }
              
              #----  ESTIMATES OF SCA EFFECTS  ----#
              for (I in 1:p)  {
                for (J in 1:p)  {
                  if (I<J)   G_SCA[I,J] <- ((1/2)* (mdata[I,J]+mdata[J,I]))-((1/(2*p))*(XI[I]+XJ[I]+XI[J]+XJ[J]))+(1/(p^2)) * SUMX
                }
              } 
              
              #----  ESTIMATES OF RECIPROCAL EFFECTS  ----#
              for (I in 1:p)  {
                for (J in 1:p)  {
                  if (I>J)     G_SCA[I,J] <- (1/2)*(mdata[J,I] - mdata[I,J])
                }
              }
              
              rownames(G_SCA) <- colnames(G_SCA) <- seq(1:p)
              # recode back to the user's notation
              G_SCA_2 <-format(round(G_SCA,4), digits=4, nsmall=4, justify="right")
              rownames(G_SCA_2) <- colnames(G_SCA_2) <- codingGuide$levelsParents[match(colnames(G_SCA_2), codingGuide$newCodeParents)]
              colnames(G_SCA_2) <- paste(p2,"=",colnames(G_SCA_2), sep="")
              rownames(G_SCA_2) <- paste(p1,"=",rownames(G_SCA_2), sep="")
              
              G_SCA_print <- noquote(gsub(" 0.0000", "", G_SCA_2))
              cat("\n\nGENERAL COMBINING ABILITY EFFECTS (diagonal), SPECIFIC COMBINING\n")
              cat("ABILITY EFFECTS (above diagonal) AND RECIPROCAL EFFECTS (below diagonal):\n\n")
              print(G_SCA_print)
              result[[i]]$site[[j]]$gcasca.matrix <- G_SCA_print
              
              #--- estimates of standard errors ---#
              
              SE_GI <- sqrt(((p-1)/(2*(p^2)))*MEPRIME)
              LSD_GI <- NA
              SE_SII <- sqrt((((p-1)^2)/(p^2))*MEPRIME)
              LSD_SII <- NA
              SE_SIJ <- sqrt((1/(2*(p^2)))*((p^2)-(2*p)+2)*MEPRIME)
              LSD_SIJ <- NA
              SE_RIJ <- sqrt((1/2)*MEPRIME)
              LSD_RIJ <- NA
              
              #SE_Gdiff <- sqrt((1/(p^2))*MEPRIME)
              SE_Gdiff <- sqrt((1/p)*MEPRIME)   #from book of Singh and Chaudhary
              LSD_Gdiff <- qt(.975,EDF)*SE_Gdiff
              
              SE_SII_SJJ <- sqrt(((2*(p-2))/p)*MEPRIME)
              LSD_SII_SJJ <- qt(.975,EDF)*SE_SII_SJJ
              
              SE_SII_SIJ <- sqrt(((3*p-2)/(2*p))*MEPRIME)
              LSD_SII_SIJ <- qt(.975,EDF)*SE_SII_SIJ
              
              SE_SII_SJK <- sqrt((3*(p-2)/(2*p))*MEPRIME)
              LSD_SII_SJK <- qt(.975,EDF)*SE_SII_SJK
              
              SE_SIJ_SIK <- sqrt(((p-1)/p)*MEPRIME)
              LSD_SIJ_SIK <- qt(.975,EDF)*SE_SIJ_SIK
              
              SE_SIJ_SKL <- sqrt(((p-2)/p)*MEPRIME)
              LSD_SIJ_SKL <- qt(.975,EDF)*SE_SIJ_SKL
              
              SE_RIJ_RKL <- sqrt(MEPRIME)
              LSD_RIJ_RKL <- qt(.975,EDF)*SE_RIJ_RKL
              
              STDERR <- round(rbind(SE_GI, SE_SII, SE_SIJ, SE_RIJ,SE_Gdiff, SE_SII_SJJ, SE_SII_SIJ, SE_SII_SJK, 
                                    SE_SIJ_SIK, SE_SIJ_SKL, SE_RIJ_RKL), digits=4)
              
              LSD <- round(rbind(LSD_GI, LSD_SII, LSD_SIJ, LSD_RIJ,LSD_Gdiff, LSD_SII_SJJ, LSD_SII_SIJ, LSD_SII_SJK, 
                                 LSD_SIJ_SIK, LSD_SIJ_SKL, LSD_RIJ_RKL),digits=4)
              
              VAREST <- as.table(cbind(STDERR, LSD))
              
              rownames(VAREST) <- c('Gi', 'Sii', 'Sij', 'Rij', 'Gi-Gj', 'Sii-Sjj', 'Sii-Sij', 'Sii-Sjk', 'Sij-Sik', 'Sij-Skl', 'Rij-Rkl')
              colnames(VAREST) <- c("Std. Error", "LSD")
              cat("\n\nTABLE OF STANDARD ERRORS AND LSDs:\n\n")
              print(VAREST)
              cat("\n\n")
              result[[i]]$site[[j]]$stderror.table <-VAREST
            } else {
              cat("\n ERROR: Too many missing values. Cannot perform test for significance of GCA and SCA effects") 
              cat("\n        and estimation of genetic variance components.\n\n")
            }
            #} ## end of if (pValue < alpha)
          } ## end of if (meansComplete)
	      } ## end of else (if (nBalance<nrow(temp.data)))
	    } ## end of else (if (responseRate < 0.80))
		} ## end of for loop (j)
	}## end of loop (i)
	cat("\n==============================================================\n")
	detach("package:doBy")
	return(list(output = result))
}

