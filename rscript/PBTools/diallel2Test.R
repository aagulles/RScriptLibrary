diallel2Test <-
function(design = c("CRD", "RCB"), data, respvar, p1, p2, rep=NULL, block=NULL, cross=TRUE, individual=NULL, environment=NULL, alpha=0.05) {
	
	options(show.signif.stars=FALSE)
	data <- eval(parse(text = data))
	
	#trim the strings
	respvar <- trim.strings(respvar)
	p1 <- trim.strings(p1)
	p2 <- trim.strings(p2)
	block <- trim.strings(block)
	rep <- trim.strings(rep)
	individual <- trim.strings(individual)
	if (!is.null(environment)) {environment <-trim.strings(environment)}
	alpha <- trim.strings(alpha)
	
	# --- create titles --- #
	if (cross) {parentsType<-"CROSS"
	} else {parentsType<-"SELF"}
	cat("\nDIALLEL ANALYSIS: GRIFFING METHOD II IN ",design, " (", parentsType, ")\n", sep="")
	
	# --- get number of environment levels --- #
	if (!is.null(environment)) {
	  data[,match(environment, names(data))] <- factor(trim.strings(data[,match(environment, names(data))]))
	  envNumLevels<-nlevels(data[,match(environment, names(data))])
	} else {envNumLevels<-1}
  
	result <- list()
	for (i in (1:length(respvar))) {
	  result[[i]] <- list()
	  cat("\n\nRESPONSE VARIABLE: ", respvar[i], "\n", sep="")
	  for (j in (1:envNumLevels)) {
	    result[[i]]$site[[j]] <- list()
	    if (!is.null(environment)) {
	      crdVars<-c(respvar[i], p1, p2, rep, environment)
	      rcbVars<-c(respvar[i], p1, p2, block, environment)
	    } else {
	      crdVars<-c(respvar[i], p1, p2, rep)
	      rcbVars<-c(respvar[i], p1, p2, block)
	    }
	    if (design == "CRD") {temp.data <- data[sort(match(crdVars, names(data)))]}
	    if (design == "RCB") {temp.data <- data[sort(match(rcbVars, names(data)))]}
	    if (!is.null(environment)) {
	      cat("\n-----------------------------")
	      cat("\nANALYSIS FOR: ",environment, " = " ,levels(temp.data[,match(environment, names(temp.data))])[j],"\n", sep="")
	      cat("-----------------------------\n")
	      temp.data <- subset(temp.data, temp.data[,match(environment, names(temp.data))] == levels(temp.data[,match(environment, names(temp.data))])[j])
	      temp.data[,match(environment, names(temp.data))] <- factor(trim.strings(temp.data[,match(environment, names(temp.data))]))
	    }
	    
	    # --- define factors --- #
			temp.data[,match(p1, names(temp.data))] <- factor(trim.strings(temp.data[,match(p1, names(temp.data))]))
			temp.data[,match(p2, names(temp.data))] <- factor(trim.strings(temp.data[,match(p2, names(temp.data))]))
			if (design == "CRD") {temp.data[,match(rep, names(temp.data))] <- factor(trim.strings(temp.data[,match(rep, names(temp.data))])) }
			if (design == "RCB") {temp.data[,match(block, names(temp.data))] <- factor(trim.strings(temp.data[,match(block, names(temp.data))]))	}
			
			# --- check if raw data is balanced. If not, generate estimates for missing values --- #
			temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[i], names(temp.data))]) == FALSE))
			if (design == "CRD") {
				tempDataForAnova<-temp.data[,c(p1, p2, rep, respvar[i])]
				balancedData<-generateBalancedData(design="DIALLEL2", data=tempDataForAnova, respvar[i], p1, p2, rep)
			}
			if (design == "RCB") {
				tempDataForAnova<-temp.data[,c(p1, p2, block, respvar[i])]
				balancedData<-generateBalancedData(design="DIALLEL2", data=tempDataForAnova, respvar[i], p1, p2, block)
			}
      
	    # --- data summary --- #
	    funcTrialSum <- class.information2(names(temp.data),temp.data)
	    cat("\nDATA SUMMARY: ","\n", sep="")
	    print(funcTrialSum)
	    cat("\nNumber of observations read: ",nrow(temp.data), sep="")
	    cat("\nNumber of missing observations: ",nrow(balancedData)-nrow(temp.data), sep="")
	    
	    # --- ANOVA for Diallel Method 2 experiment --- #
	    estimatedMissing <- FALSE
	    #cat("\n\n\nANOVA TABLE FOR THE EXPERIMENT\n")
	    if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
				if (nrow(temp.data) == nrow(balancedData)) {
				  anovaRemark <- "REMARK: Raw dataset is balanced."
					dataForAnova<-tempDataForAnova  
				} else {
					if (design == "CRD") {dataForAnova<-estimateMissingData(design="CRD", data=balancedData, respvar[i], p1, p2, rep)  }
					if (design == "RCB") {dataForAnova<-estimateMissingData(design="RCB", data=balancedData, respvar[i], p1, p2, block)  }
					anovaRemark <- "REMARK: Raw data and estimates of the missing values are used."
					estimatedMissing <- TRUE
				}  
			} else {anovaRemark <- "ERROR: Too many missing values. Cannot perform ANOVA for balanced data." 
			        #cat("\n",anovaRemark)
			}		
      
			# --- testing for genotypic effect --- #
			pValue=0;
			if (design == "CRD") {
        myformula1 <- paste(respvar[i], " ~ ", p1,":",p2, sep = "")
        model <- aov(formula(myformula1), data = temp.data)
        anovatable<-summary(model)
        
        # --- format anova table --- #
        p<-formatC(as.numeric(format(anovatable[[1]][1,5], scientific=FALSE)), format="f")
        anova_new<-cbind(round(anovatable[[1]][1:4],digits=2), rbind(p, " "))
        anova_print<-replace(anova_new, is.na(anova_new), " ")
        pValue <- anovatable[[1]]$"Pr(>F)"[[1]] 
      }
			if (design == "RCB") {
        myformula1 <- paste(respvar[i], " ~ ", block," + ", p1,":",p2, sep = "") 
        model <- aov(formula(myformula1), data = temp.data)
        anovatable<-summary(model)
        
        # --- format anova table --- #
        pval1<-formatC(as.numeric(format(anovatable[[1]][1,5], scientific=FALSE)), format="f")
        pval2<-formatC(as.numeric(format(anovatable[[1]][2,5], scientific=FALSE)), format="f")
        anova_new<-cbind(round(anovatable[[1]][1:4],digits=2), rbind(pval1, pval2, " "))
        anova_print<-replace(anova_new, is.na(anova_new), " ")
        pValue <- anovatable[[1]]$"Pr(>F)"[[2]]
      }

			colnames(anova_print) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
			cat("\n\n\nANOVA TABLE FOR TESTING SIGNIFICANCE OF GENOTYPIC EFFECT \n", sep="")
			cat(" Formula: ", myformula1,"\n\n")
			print(anova_print)
			result[[i]]$site[[j]]$genoEffect.anova <-anova_print
		
			out <- anova(model)
			EMS <- out[rownames(out)=="Residuals","Mean Sq"]
			EDF <- out[rownames(out)=="Residuals","Df"]
		
			# --- mean data of diallel --- #   
			library(doBy)
			myformula2<- paste(respvar[i], " ~ ", p1," + ",p2, sep = "")
              
			if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
			  meandata <- summaryBy(formula(myformula2), data=dataForAnova)
			  
        # --- serial to parallel of meandata --- #
			  p <- nlevels(temp.data[,match(p1, names(temp.data))])
			  if (design == "CRD") {r <- nlevels(temp.data[,match(rep, names(temp.data))])}
			  if (design == "RCB") {r <- nlevels(temp.data[,match(block, names(temp.data))]) }
			  n <- 1
			  
			  mdata <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
			  dim(mdata) <- c(p,p)
			  
			  for (I in 1:p)  {
			    for (J in 1:p)   {
			      if (I <= J) mdata[I,J] <- meandata[(meandata[,match(p1, names(meandata))]==I & meandata[,match(p2, names(meandata))]==J),3]
			    } 
			  }
			  
			  colnames(mdata) <- rownames(mdata) <- levels(temp.data[,match(p1, names(temp.data))])
			  
			  # --- printing the matrix of means --- #
			  mdata2 <- format(round(mdata, 4),nsmall=4)
			  mdata3 <-noquote(format(gsub(" 0.0000", "", mdata2),justify="right"))
			  cat("\n\nMATRIX OF MEANS:\n")
			  print(mdata3)
			  result[[i]]$site[[j]]$means.matrix  <- mdata3
			  if (estimatedMissing) {
			    cat("-------\n")
			    cat(anovaRemark,"\n")
			  }
      } 

			# --- check if genotypic effect is significant. If significant, proceed to diallel analysis --- #
			alpha <- as.numeric(alpha)
			if (pValue < alpha) {
			  cat("\n\nANALYSIS OF VARIANCE:")
			  if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
			    mirror <- mdata
			    
			    for (I in 1:p)  {
			      for (J in 1:p)   {
			        if (I > J) mirror[I,J] <- mirror[J,I]
			      } 
			    }
			    
			    XI <- rowSums(mirror)
			    minus <- c(0, rep(p))
			    plus <- c(0, rep(p))
			    for (I in 1:p) {
			      minus[I] <- XI[I] - mdata[I,I]
			      plus[I] <- XI[I] + mdata[I,I]
			    }
			    
			    SUMX <- sum(mdata)
			    MEPRIME <-EMS/r            #  ERROR PRIME (ME')  
			    
			    #----  GENERAL COMBINING ABILITY (GCA) SUM OF SQUARES  ----#
			    SG <- (1/(p+2))*(sum(plus^2)-((4/p)*SUMX^2))
			    
			    #----- SPECIFIC COMBINING ABILITY (SCA) SUM OF SQUARES  ----#
			    SS <- sum(mdata^2) - ((1/(p+2))*sum(plus^2))+((2*SUMX^2)/((p+1)*(p+2)))
			    
			    #----  ERROR PRIME SUM OF SQUARES  ----#
			    SE <- MEPRIME * EDF
			    
			    #----  COMPUTATION OF MEAN SQUARE AND F-VALUES  ----#
			    DG <- p-1
			    DS <- p*(p-1)/2
			    DE <- EDF
			    
			    MG <- SG/DG
			    MS <- SS/DS
			    
			    FG <- MG/MS
			    FS <- MS/MEPRIME
			    FE <- NA
			    
			    PG <- 1 - pf(FG, DG, DS)
			    PS <- 1 - pf(FS, DS, EDF)
			    PE <- NA
			    
			    SV <- c("SV","GCA", "SCA", "Error")
			    DF <- format(rbind("Df",DG,DS,DE), justify="right")
			    SSq <- format(rbind("Sum Sq",format(round(rbind(SG, SS, SE),digits=2))), justify="right")
			    MSq <- format(rbind("Mean Sq",format(round(rbind(MG, MS, MEPRIME),digits=2))), justify="right")
			    Fvalue <- format(round(rbind(FG, FS),digits=2))
			    Fvalue2 <- format(rbind("F value", Fvalue, " "), justify="right")
			    PG2<-formatC(as.numeric(format(PG,scientific=FALSE)),format="f")
			    PS2<-formatC(as.numeric(format(PS,scientific=FALSE)),format="f")
			    P <- format(rbind("Pr(>F)",PG2,PS2, " "), justify="right")
			    AOV <- noquote(cbind(SV,DF,SSq,MSq,Fvalue2,P))
			    colnames(AOV) <- c("", "", "", "", "", "")
			    rownames(AOV) <- rep("",nrow(AOV))
			    print(AOV)
			    result[[i]]$site[[j]]$gcasca.anova <-AOV
			    
			    #--- Estimation of varinace components ---#
			    
			    Ve <- MEPRIME
			    SEe <- sqrt((2/EDF)*(MEPRIME^2))
			    
			    Vs <- (MS-MEPRIME)
			    if (Vs < 0) Vs <- 0
			    SEs <- sqrt((4/(p*(p-1)))*MS^2 + (2/EDF)*MEPRIME)
			    
			    Vg <- (1/(p+2))*(MG - MS)	 
			    if (Vg < 0) Vg <- 0
			    SEg <- sqrt((2/((p+1)*(p+2)))*MG^2 + (4/(p*(p+1)*(p+2)))*MS^2)
			    
			    VC <- round(rbind(Vg, Vs, Ve), digits=4)
			    StdErr <- round(rbind(SEg, SEs, SEe), digits=4)
			    TABLE <- cbind(VC, StdErr)
			    colnames(TABLE) <- c("Estimate", "Std. Error")
			    rownames(TABLE) <- c(" GCA", " SCA", " Error")
			    TABLE_print<-as.table(TABLE)
			    cat("\n\nESTIMATES OF VARIANCE COMPONENTS:\n")
			    print(TABLE_print)
			    result[[i]]$site[[j]]$var.components <-TABLE_print
			    
			    #---- Genetic Variance components ----#
			    if (cross) {F<-0}
			    else {F<-1}
          
			    VA <- (4/(1+F))*Vg
			    VVA <- 4*SEg^2
			    
			    VD <- (4/(1+F)^2)*Vs
			    if (VD < 0) VD <- 0
			    VVD <- SEs^2
			    
			    VE <- Ve
			    VP <- VA + VD + VE
			    h2B <- (VA + VD) / VP                 # individual based
			    # Vh2B <- (2*(1-h2B)^2*(1+(r-1)*h2B)^2)/(r*(r-1)*p*(p+1)-1)
			    
			    h2N <- VA / VP                 # individual based
			    # Vh2N <- (2*(1-h2N)^2*(1+(r-1)*h2N)^2)/(r*(r-1)*p*(p+1)-1)
			    Dominance.ratio <- sqrt(2*VD/VA)   
			    
			    Estimate <- round(rbind(VA, VD, h2N, h2B, Dominance.ratio), digits=4)
			    colnames(Estimate) <- c("Estimate")
			    rownames(Estimate) <- c(" VA", " VD", " h2-narrow sense", " H2-broad sense", " Dominance Ratio")
			    TABLE2 <- as.table(Estimate)
			    cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:\n")
			    print(TABLE2)
			    result[[i]]$site[[j]]$genvar.components <-TABLE2 
			    
			    #----  ESTIMATES OF GCA EFFECTS  ----#
			    G_SCA <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
			    dim(G_SCA) <- c(p,p)
			    
			    for (I in 1:p)  {
			      for (J in 1:p)   {
			        if (I == J) G_SCA[I,J] <- (1/(p+2))* (XI[I] + mdata[I,J] - (2/p)*SUMX)
			      } 
			    }
			    
			    #----  ESTIMATES OF SCA EFFECTS  ----#
			    for (I in 1:p)  {
			      for (J in 1:p)   {
			        if (I < J)   {
			          B1 <- mdata[I,J]
			          B2 <- (1/(p+2))* (plus[I] + plus[J])
			          B3 <- (2/((p+1)*(p+2)))*SUMX
			          G_SCA[I,J] <- B1 - B2 + B3
			        }
			      } 
			    }
			    
			    colnames(G_SCA) <- rownames(G_SCA) <- levels(temp.data[,match(p1, names(temp.data))])
			    G_SCA_2 <-format(round(G_SCA, 4),nsmall=4, justify="right")
			    G_SCA_3 <- noquote(gsub(" 0.0000", "", G_SCA_2))
			    cat("\n\nGENERAL COMBINING ABILITY EFFECTS (diagonal), SPECIFIC COMBINING\n")
			    cat("ABILITY EFFECTS (above diagonal):\n")
			    print(G_SCA_3)
			    result[[i]]$site[[j]]$gcasca.matrix <-G_SCA_3
			    
			    #--- estimates of standard errors ---#
			    
			    SE_GI <- sqrt((p-1)*Ve/(p*(p+2)))
			    LSD_GI <- NA
			    
			    #------------------------------------------#
			    # Note: in the book of singh and chaudhary #
			    # the formula for SE_SII and SE_SIJ are    #
			    # interchanged.                            #
			    #------------------------------------------#
			    
			    SE_SII <- sqrt(p*(p-1)*Ve/((p+1)*(p+2)))
			    LSD_SII <- NA
			    
			    SE_SIJ <- sqrt(((p^2+p+2)*Ve)/((p+1)*(p+2)))
			    LSD_SIJ <- NA
			    SE_Gdiff <- sqrt(2*Ve/(p+2))
			    LSD_Gdiff <- qt(.975,EDF)*SE_Gdiff
			    
			    SE_SII_SJJ <- sqrt(2*(p-2)*MEPRIME/(p+2))
			    LSD_SII_SJJ <- qt(.975,EDF)*SE_SII_SJJ
			    
			    SE_SIJ_SIK <- sqrt(2*(p+1)*Ve/(p+2))
			    LSD_SIJ_SIK <- qt(.975,EDF)*SE_SIJ_SIK
			    
			    SE_SIJ_SKL <- sqrt(2*p*MEPRIME/(p+2))
			    LSD_SIJ_SKL <- qt(.975,EDF)*SE_SIJ_SKL
			    
			    STDERR <- round(rbind(SE_GI,SE_SII, SE_SIJ, SE_Gdiff, SE_SII_SJJ, SE_SIJ_SIK, SE_SIJ_SKL), digits=4)
			    LSD <- round(rbind(LSD_GI,LSD_SII, LSD_SIJ, LSD_Gdiff, LSD_SII_SJJ, LSD_SIJ_SIK, LSD_SIJ_SKL), digits=4)
			    VAREST <- as.table(cbind(STDERR, LSD))
			    
			    rownames(VAREST) <- c('Gi', 'Sii', 'Sij', 'Gi-Gj', 'Sii-Sjj', 'Sij-Sik', 'Sij-Skl')
			    colnames(VAREST) <- c("Std. Error", "LSD")
			    cat("\n\nTABLE OF STANDARD ERRORS AND LSDs:\n")
			    print(VAREST)
			    result[[i]]$site[[j]]$stderror.table <-VAREST
			  }else {
			    cat("\n\n ERROR: Too many missing values. Cannot perform test for significance of GCA and SCA effects") 
			    cat("\n        and estimation of genetic variance components.\n\n")
			  }
			} ## end of if statement
		} ## end of for loop (j)
	cat("\n==============================================================\n")
	}## end of loop (i)
	detach("package:doBy")
	return(list(output = result))
}

