#----------------------------------------------------------------
# Half-sib mating, design NC III - Multi- environment       
# F2 lines (males) are randomly selected and back-crossed   
# to parents (females)                                      
# Estimates genetic variance components                     
# 
# ARGUMENTS:
# design - experiment design
# data - name of data frame
# respvar - vector of response variables
# tester - string, name of tester factor
# f2lines - string, name of F2 lines factor
# rep - string, name of rep factor
# block - string, name of block factor
# row - string, name of row factor
# column - string, name of column factor
# individual - string, name of individual factor
# environment - string, name of environment factor
#                                                     
# Script Created by: Violeta Bartolome   08.2011
# Script Modified by: Nellwyn Sales and Guoyou Ye
#----------------------------------------------------------------

nc3TestME <-function(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, tester, f2lines, rep=NULL, block=NULL, row=NULL, column=NULL, individual=NULL, environment=NULL) UseMethod("nc3TestME")
  
nc3TestME.default <-function(design = c("CRD", "RCB", "Alpha", "RowColumn"), data, respvar, tester, f2lines, rep=NULL, block=NULL, row=NULL, column=NULL, individual=NULL, environment=NULL) {
	
  options(show.signif.stars=FALSE)
  data <- eval(parse(text = data))
  library(lme4)
	
  # --- trim the strings --- #
	respvar <- trimStrings(respvar)
	tester <- trimStrings(tester)
	f2lines <- trimStrings(f2lines)
  if (!is.null(block)) {block <- trimStrings(block) }
  if (!is.null(rep)) {rep <- trimStrings(rep) }
  if (!is.null(row)) {row <- trimStrings(row) }
  if (!is.null(column)) {column <- trimStrings(column) }
  if (!is.null(individual)) {individual <- trimStrings(individual) }
	environment <-trimStrings(environment)
  
  # --- create titles --- #
  if (design == "CRD") { designName<-"CRD"}
  if (design == "RCB") { designName<-"RCB"}
  if (design == "Alpha") { designName<-"ALPHA-LATTICE"}
  if (design == "RowColumn") { designName<-"ROW-COLUMN"}
  
  cat("\nMULTIPLE ENVIRONMENT ANALYSIS\n")
  cat("\nDESIGN: NORTH CAROLINA EXPERIMENT III IN ",designName, "\n", sep="")
	 
	result <- list()
	for (i in (1:length(respvar))) {
		result[[i]] <- list()
		temp.data <- data[(is.na(data[,match(respvar[i], names(data))]) == FALSE),]
    
		cat("\n-----------------------------")
		cat("\nRESPONSE VARIABLE: ", respvar[i], "\n", sep="")
		cat("-----------------------------\n")
		
		responseRate<-(nrow(temp.data)/nrow(data))
		
		if (responseRate < 0.80) {
		  cat("\nToo many missing observations. Cannot proceed with the analysis.\n\n")
		  next
		} else {
		  # --- define factors and get number of levels --- #
		  temp.data[,match(tester, names(temp.data))] <- factor(trimStrings(temp.data[,match(tester, names(temp.data))]))
		  temp.data[,match(f2lines, names(temp.data))] <- factor(trimStrings(temp.data[,match(f2lines, names(temp.data))]))
		  temp.data[,match(environment, names(temp.data))] <- factor(trimStrings(temp.data[,match(environment, names(temp.data))]))
		  
		  f<-nlevels(temp.data[,match(tester, names(temp.data))])
		  m<-nlevels(temp.data[,match(f2lines, names(temp.data))])
		  e<-nlevels(temp.data[,match(environment, names(temp.data))])
		  
		  # --- compute for harmonic mean --- #
		  lengthPerCross<-tapply(temp.data[,respvar[i]], temp.data[,c(tester, f2lines, environment)], length)
		  repHarmonicMean<-1/mean(1/lengthPerCross, na.rm=TRUE)
		  
		  if (design == "CRD") {
		    nlevelsRep<-max(lengthPerCross, na.rm=TRUE)
		  }
		  if (design == "RCB") {
		    temp.data[,match(block, names(temp.data))] <- factor(trimStrings(temp.data[,match(block, names(temp.data))]))  
		    nlevelsRep<-nlevels(temp.data[,match(block, names(temp.data))])
		  }
		  if (design == "Alpha") {
		    temp.data[,match(rep, names(temp.data))] <- factor(trimStrings(temp.data[,match(rep, names(temp.data))]))
		    temp.data[,match(block, names(temp.data))] <- factor(trimStrings(temp.data[,match(block, names(temp.data))]))
		    nlevelsRep<-nlevels(temp.data[,match(rep, names(temp.data))])
		  }
		  if (design == "RowColumn") {
		    temp.data[,match(rep, names(temp.data))] <- factor(trimStrings(temp.data[,match(rep, names(temp.data))]))
		    temp.data[,match(row, names(temp.data))] <- factor(trimStrings(temp.data[,match(row, names(temp.data))]))
		    temp.data[,match(column, names(temp.data))] <- factor(trimStrings(temp.data[,match(column, names(temp.data))]))
		    nlevelsRep<-nlevels(temp.data[,match(rep, names(temp.data))])
		  }
		  
		  nBalance<-2*m*e*nlevelsRep
		  
		  # --- data summary --- #
		  funcTrialSum <- class.information2(names(temp.data),temp.data)
		  cat("\nDATA SUMMARY: ","\n\n", sep="")
		  print(funcTrialSum)
		  cat("\n Number of observations read: ",nrow(data), sep="")
		  cat("\n Number of observations used: ",nrow(temp.data), sep="")
		  missingObs<-nBalance-nrow(temp.data)
		  
		  # --- if design is CRD or RCB, check if raw data is balanced. --- #
		  # --- give warning if the number of obs in the dataset exceeds the balanced data size --- #
		  if (nBalance<nrow(temp.data)) {
		    cat("\n\n***\nERROR: The number of observations read in the data exceeds the size of a balanced data.\n***\n\n")
		  } else {
		    if (design == "CRD" || design == "RCB") {
		      cat("\n\n\nANOVA TABLE:\n\n")
		      if ((nrow(temp.data)/nBalance) >= 0.90) {
		        if (nrow(temp.data) == nBalance) {
		          anovaRemark <- "REMARK: Raw dataset is balanced."
		          dataForAnova<-temp.data  
		        } else {
		          #for CRD, no estimation of missing values
		          if (design == "CRD") {
		            dataForAnova <- temp.data
		            anovaRemark <- "REMARK: Raw dataset is unbalanced."
		          }
		          #for RCB, estimate missing values, if there are any
		          if (design == "RCB") {
		            dataForAnova<-NULL
		            for (x in (1:e)) {
		              tempSplit<-temp.data[temp.data[,environment]==levels(temp.data[,match(environment, names(temp.data))])[x],]
		              tempDataForAnova<-tempSplit[,c(f2lines, tester, block, respvar[i])]
		              balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], f2lines, tester, block)
		              tempDataForAnova<-estimateMissingData(design="RCB", data=balancedData, respvar[i], f2lines, tester, block)
		              tempDataForAnova<-tempDataForAnova[,c(f2lines, tester, block, respvar[i])]
		              tempDataForAnova$EnvIndex<-levels(temp.data[,match(environment, names(temp.data))])[x]
		              dataForAnova<-rbind(dataForAnova,tempDataForAnova)
		            }
		            colnames(dataForAnova)[match("EnvIndex",names(dataForAnova))]<-environment
		            dataForAnova[,match(environment, names(dataForAnova))] <- factor(dataForAnova[,match(environment, names(dataForAnova))])
		            
		            anovaRemark <- "REMARK: Raw data and estimates of the missing values are used."
		          }
		        }
		        
		        # --- ANOVA for NC3 experiment --- #
		        if (design == "CRD") {
		          myformula <- paste(respvar[i], " ~ ", environment, " * ", f2lines, "*", tester, sep = "") 
		          anovaFactorial<-summary(aov(formula(myformula), data=dataForAnova))
		          
		          #rerun aov using temp.data to get the original df's
		          anovaFactorial.temp<-summary(aov(formula(myformula), data=temp.data))
		          anovaFactorial<-adjustAnovaDf(anovaFactorial, anovaFactorial.temp[[1]]$"Df")
		        }
		        if (design == "RCB") {
		          myformula <- paste(respvar[i], " ~ ", environment, ":", block, " + ", f2lines, "*", tester, "*", environment, sep = "")
		          anovaFactorial<-summary(aov(formula(myformula), data=dataForAnova))
		          
		          #rerun aov using temp.data to get the original df's
		          anovaFactorial.temp<-summary(aov(formula(myformula), data=temp.data))
		          anovaFactorial.temp[[1]]$"Df"[length(anovaFactorial.temp[[1]]$"Df")]<-anovaFactorial[[1]]$"Df"[length(anovaFactorial[[1]]$"Df")]-missingObs
		          anovaFactorial<-adjustAnovaDf(anovaFactorial, anovaFactorial.temp[[1]]$"Df")
		          
		          #rearrange the rows of anovaFactorial
		          index<-match(paste(environment, ":", block, sep=""), trimStrings(rownames(anovaFactorial)))
		          indexEnv<-match(paste(environment), trimStrings(rownames(anovaFactorial)))
		          anovaFactorial <- rbind(anovaFactorial[c(indexEnv,index),], anovaFactorial[-I(match(c(indexEnv, index), 1:nrow(anovaFactorial))),])
		          
		          #recompute f value and pvalue of environment
		          anovaFactorial[1,"F value"] <- anovaFactorial[1, "Mean Sq"]/anovaFactorial[2, "Mean Sq"]
		          anovaFactorial[1,"Pr(>F)"] <-  pf(anovaFactorial[1,"F value"], anovaFactorial[1,"Df"], anovaFactorial[2,"Df"], lower.tail = FALSE)
		        }
		        anovaFormat<-formatAnovaTable(anovaFactorial)
		        print(anovaFormat)
		        cat("-------\n")
		        cat(anovaRemark)
		        result[[i]]$nc3.anova <- anovaFormat
		        
		      } else {anovaRemark <- "ERROR: Too many missing values. Cannot perform ANOVA for balanced data." 
		              cat(anovaRemark)
		      }
		    }
		    
		    # --- LMER for the design --- #
		    if (design == "CRD") {myformula1 <- paste(respvar[i], " ~ 1 + (", tester, ") + (1|", environment, ") + (1|", f2lines, ") + (1|", tester, ":", f2lines, ") + (1|", environment, ":", tester, ") + (1|", environment, ":", f2lines, ") + (1|", environment, ":", tester, ":", f2lines,")", sep = "") }
		    if (design == "RCB") {myformula1 <- paste(respvar[i], " ~ 1 + (", tester, ") + (1|", environment, ") + (1|", environment, ":", block, ") + (1|", f2lines, ") + (1|", tester, ":", f2lines, ") + (1|", environment, ":", tester, ") + (1|", environment, ":", f2lines, ") + (1|", environment, ":", tester, ":", f2lines,")", sep = "") }
		    if (design == "Alpha") {myformula1 <- paste(respvar[i], " ~ 1 + (", tester, ") + (1|", environment, ") + (1|", environment, ":", rep, ") + (1|", environment, ":", rep, ":", block, ") + (1|", f2lines, ") + (1|", tester, ":", f2lines, ") + (1|", environment, ":", tester, ") + (1|", environment, ":", f2lines, ") + (1|", environment, ":", tester, ":", f2lines,")", sep = "") }
		    if (design == "RowColumn") {myformula1 <- paste(respvar[i], " ~ 1 + (", tester, ") + (1|", environment, ") + (1|", environment, ":", rep, ") + (1|", environment, ":", rep, ":", row, ") + (1|", environment, ":", rep, ":", column, ") + (1|", f2lines, ") + (1|", tester, ":", f2lines, ") + (1|", environment, ":", tester, ") + (1|", environment, ":", f2lines, ") + (1|", environment, ":", tester, ":", f2lines,")", sep = "") }
		    
		    model <- lmer(formula(myformula1), data = temp.data)
		    result[[i]]$lmer.result <- summary(model)
		    
		    # --- edit format of lmer output before printing --- #
		    fixeff.table<-round(summary(model)@coefs, digits=4)
        rownames(fixeff.table) <- paste(" ", rownames(fixeff.table), sep="")
		    remat<-summary(model)@REmat
		    Groups<-remat[,1]
		    Variance<-formatC(as.numeric(remat[,3]), format="f")
		    Std.Deviation<-formatC(as.numeric(remat[,4]), format="f")
		    Variance2<-format(rbind("Variance", cbind(Variance)), justify="right")
		    Std.Deviation2<-format(rbind("Std. Deviation", cbind(Std.Deviation)), justify="right")
		    Groups2<-format(rbind("Groups",cbind(Groups)), justify="left")
		    rematNew<-noquote(cbind(Groups2, Variance2, Std.Deviation2))
		    colnames(rematNew)<-c("", "", "")
		    rownames(rematNew)<-rep("",nrow(rematNew))
		    cat("\n\n\nLINEAR MIXED MODEL FIT BY RESTRICTED MAXIMUM LIKELIHOOD:\n\n")
		    cat(" Formula: ", myformula1,"\n\n")
		    print(summary(model)@AICtab) 
		    #cat("\n Fixed Effects:\n")
		    #print(fixeff.table)
		    cat("\n Random Effects:")
		    print(rematNew)
		    
		    #--- Estimates of genetic variance components ---#
		    varcomp <- summary(model)@REmat
		    Ve <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
		    Vf_m <- as.numeric(varcomp[varcomp[,1] == paste(tester,":",f2lines,sep=""), "Variance"])
		    Vm <- as.numeric(varcomp[varcomp[,1] == f2lines, "Variance"])
		    Vef_m <- as.numeric(varcomp[varcomp[,1] == paste(environment,":",tester,":",f2lines,sep=""), "Variance"])
		    Vem <- as.numeric(varcomp[varcomp[,1] == paste(environment,":",f2lines,sep=""), "Variance"])
		    Vef <- as.numeric(varcomp[varcomp[,1] == paste(environment,":",tester,sep=""), "Variance"])
		    
		    F<-0
		    r<-repHarmonicMean
		    
		    VA <- 4*Vm
		    VD <- 2*Vf_m
		    VAE <- 4*Vem
		    VDE <- 2*Vef_m
		    if (VD < 0 || VD < 1e-10) VD <- 0
		    if (VA < 0 || VA < 1e-10) VA <- 0
		    if (VAE < 0 || VAE < 1e-10) VAE <- 0
		    if (VDE < 0 || VDE < 1e-10) VDE <- 0
		    VE <- Ve
		    VP <- VA + VD + VAE + VDE + VE
		    h2B <- (VA + VD) / VP            	    # individual based
		    h2N <- VA / VP             		        # individual based
		    h2f <- Vm /(Vm + VE/(r*e))   		      # family heritability
		    Dominance.ratio <- sqrt(2*VD/VA)
		    
		    VA_p<-formatNumericValue(VA)
		    VAE_p<-formatNumericValue(VAE)
		    VD_p<-formatNumericValue(VD)
		    VDE_p<-formatNumericValue(VDE)
		    h2N_p<-formatNumericValue(h2N)
		    h2B_p<-formatNumericValue(h2B)
		    h2f_p<-formatNumericValue(h2f)
		    Dominance.ratio_p<-formatNumericValue(Dominance.ratio)
		    
		    Estimate <- rbind(VA_p, VAE_p,  VD_p, VDE_p, h2N_p, h2B_p, h2f_p, Dominance.ratio_p)
		    with.colheader<-format(rbind("Estimate", Estimate), justify="right")
		    colnames(with.colheader) <- c("")
		    rownames(with.colheader) <- c("", " VA", " VAxE", " VD", " VDxE", " h2-narrow sense", " H2-broad sense", " h2-family", " Dominance Ratio")
		    TABLE <- as.table(with.colheader)
		    cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:\n")
		    print(TABLE)
		    result[[i]]$genvar.components <- TABLE
		    cat("\n")
		  }
		}
	}## end of loop (i)
  cat("\n==============================================================\n")
	detach("package:lme4")
  return(list(output = result))
}

