diallel1Test <-
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
	cat("\nDIALLEL ANALYSIS: GRIFFING METHOD I IN ",design, " (", parentsType, ")\n", sep="")
	
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
				balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], p1, p2, rep)
			}
			if (design == "RCB") {
				tempDataForAnova<-temp.data[,c(p1, p2, block, respvar[i])]
				balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], p1, p2, block)
			}
			
	    # --- data summary --- #
	    funcTrialSum <- class.information2(names(temp.data),temp.data)
	    cat("\nDATA SUMMARY: ","\n", sep="")
	    print(funcTrialSum)
	    cat("\nNumber of observations read: ",nrow(temp.data), sep="")
	    cat("\nNumber of missing observations: ",nrow(balancedData)-nrow(temp.data), sep="")
	    
	    # --- ANOVA for Diallel Method 1 experiment --- #
      estimatedMissing <- FALSE
	    cat("\n\n\nANOVA TABLE FOR THE EXPERIMENT\n")
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

				if (design == "CRD") {myformula <- paste(respvar[i], " ~ ", p1, "*", p2, sep = "")  }
				if (design == "RCB") {myformula <- paste(respvar[i], " ~ ", block, " + ", p1, "*", p2, sep = "")  }
				anova.factorial<-summary(aov(formula(myformula), data=dataForAnova))		
				print(anova.factorial)
				cat("-------\n")
				cat(anovaRemark)
				result[[i]]$site[[j]]$diallel1.anova <- anova.factorial

			} else {anovaRemark <- "ERROR: Too many missing values. Cannot perform ANOVA for balanced data." 
			        cat("\n",anovaRemark)
			}
					
			# --- testing for genotypic effect ---#
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
			
			# --- mean data of full diallel --- #
			library(doBy)
			myformula2<- paste(respvar[i], " ~ ", p1," + ",p2, sep = "")
      
	    if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
	      meandata <- summaryBy(formula(myformula2), data=dataForAnova)
	    } else {
	      meandata <- summaryBy(formula(myformula2), data=temp.data)
	    }
	    respvardotmean<-paste(respvar[i],".mean", sep = "")
	    mdata <- reshape(meandata, v.names=respvardotmean, idvar=p1, timevar=p2, direction="wide", sep=".")
	    
	    # --- printing the matrix of means --- #
	    respvardotmeandot<-paste(respvardotmean,".", sep="")
	    p2equals<-paste(p2,"=", sep="")	
	    colnames(mdata) <- gsub(respvardotmeandot, p2equals, colnames(mdata))
	    mdata_print<-round(data.matrix(mdata), digits=4)
	    rownames(mdata_print)<-rep("",nrow(mdata_print))
	    cat("\n\nMATRIX OF MEANS:\n")
	    print(mdata_print)
	    result[[i]]$site[[j]]$means.matrix <-mdata_print
      if (estimatedMissing) {
	      cat("-------\n")
	      cat(anovaRemark,"\n")
      }
	    mdata[,match(p1, names(mdata))] <- NULL

	    # --- check if genotypic effect is significant. If significant, proceed to diallel analysis --- #
	    alpha <- as.numeric(alpha)
	    if (pValue < alpha) {
	      cat("\n\nANALYSIS OF VARIANCE:")
	      if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
	        p <- nlevels(temp.data[,match(p1, names(temp.data))])
	        if (design == "CRD") {r <- nlevels(temp.data[,match(rep, names(temp.data))])}
	        if (design == "RCB") {r <- nlevels(temp.data[,match(block, names(temp.data))]) }
	        
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
	        df.star <- ((EDF*p^2)*(p-1)*(p+2*(p-1)*k)^2)/((p^2)*(p-1)*((1-A)^2)+(2*EDF*A^2*(p^2+2*c.prime*k)))
	        
	        FG <- MG/MS.star
	        FS <- MS/MEPRIME
	        FR <- MR/MEPRIME
	        FE <- NA
	        
	        PG <- 1-pf(FG, DG, df.star)
	        PS <- 1-pf(FS, DS, EDF)
	        PR <- 1-pf(FR, DR, EDF)
	        PE <- NA
	        
	        # --- printing the anova table --- #
	        SV <- c("SV","GCA", "SCA", "Reciprocal", "Error")
	        DF <- format(rbind("Df",DG,DS,DR,DE), justify="right")
	        SSq <- format(rbind("Sum Sq",format(round(rbind(SG, SS, SR, SE),digits=2))), justify="right")
	        MSq <- format(rbind("Mean Sq",format(round(rbind(MG, MS, MR, MEPRIME),digits=2))), justify="right")
	        Fvalue<-format(round(rbind(FG, FS, FR),digits=2))
	        Fvalue2 <- format(rbind("F value", Fvalue, " "), justify="right")
	        PG2<-formatC(as.numeric(format(PG,scientific=FALSE)),format="f")
	        PS2<-formatC(as.numeric(format(PS,scientific=FALSE)),format="f")
	        PR2<-formatC(as.numeric(format(PR,scientific=FALSE)),format="f")
	        P <- format(rbind("Pr(>F)",PG2,PS2,PR2, " "), justify="right")
	        AOV <- noquote(cbind(SV,DF,SSq,MSq,Fvalue2,P))
	        colnames(AOV) <- c("", "", "", "", "", "")
	        rownames(AOV) <- rep("",nrow(AOV))
	        print(AOV)
	        cat(" -----\n")
	        cat(" NOTE: MS* = ",round(MS.star,digits=6), "   Error used for GCA MS with df = ",round(df.star, digits=1), "\n", sep="")
	        result[[i]]$site[[j]]$gcasca.anova <-AOV  
          
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
	        colnames(TABLE) <- c("Variance Component", "Std. Error")
	        rownames(TABLE) <- c(" GCA", " SCA", " Reciprocal", " Error")
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
	        # Vh2B <- (2*(1-h2B)^2*(1+(r-1)*h2B)^2)/(r*(r-1)*(p^2-1))
	        
	        h2N <- VA / VP                 # individual based
	        # Vh2N <- (2*(1-h2N)^2*(1+(r-1)*h2N)^2)/(r*(r-1)*(p^2-1))
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
	        
	        colnames(G_SCA) <- rownames(G_SCA) <- levels(temp.data[,match(p1, names(temp.data))])
	        G_SCA_print <- round(G_SCA, 4)
	        cat("\n\nGENERAL COMBINING ABILITY EFFECTS (diagonal), SPECIFIC COMBINING\n")
	        cat("ABILITY EFFECTS (above diagonal) AND RECIPROCAL EFFECTS (below diagonal):\n")
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
	        
	        rownames(VAREST) <- c(' Gi', ' Sii', ' Sij', ' Rij', ' Gi-Gj', ' Sii-Sjj', ' Sii-Sij', ' Sii-Sjk', ' Sij-Sik',
	                              ' Sij-Skl', ' Rij-Rkl')
	        colnames(VAREST) <- c("Std. Error", "LSD")
	        cat("\n\nTABLE OF STANDARD ERRORS AND LSDs:\n")
	        print(VAREST)
	        result[[i]]$site[[j]]$stderror.table <-VAREST
	      } else {
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

