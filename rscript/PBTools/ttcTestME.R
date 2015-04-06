ttcTestME <-
function(design = c("CRD", "RCB"), data, respvar, tester, f2lines, rep=NULL, block=NULL, individual=NULL, environment=NULL, codeP1, codeP2, codeF1, alpha=0.05) {
	
	options(show.signif.stars=FALSE)
	options(scipen=100)
	data <- eval(parse(text = data))
	
	# --- trim the strings --- #
	respvar <- trim.strings(respvar)
	tester <- trim.strings(tester)
	f2lines <- trim.strings(f2lines)
	block <- trim.strings(block)
	rep <- trim.strings(rep)
	individual <- trim.strings(individual)
	environment <-trim.strings(environment)
	codeP1 <- trim.strings(codeP1)
	codeP2 <- trim.strings(codeP2)
	codeF1 <- trim.strings(codeF1)
	alpha <- trim.strings(alpha)
	
	# --- create titles --- #
	cat("\nMULTIPLE ENVIRONMENT ANALYSIS\n")
	cat("\nDESIGN: TRIPLE TEST CROSS (NO PARENTS) IN ",design, "\n", sep="")
  
	#define as factors
	data[,match(tester, names(data))] <- factor(trim.strings(data[,match(tester, names(data))]))
	data[,match(f2lines, names(data))] <- factor(trim.strings(data[,match(f2lines, names(data))]))
	data[,match(environment, names(data))] <- factor(trim.strings(data[,match(environment, names(data))]))
	if (design == "CRD") {data[,match(rep, names(data))] <- factor(trim.strings(data[,match(rep, names(data))])) }
	if (design == "RCB") {data[,match(block, names(data))] <- factor(trim.strings(data[,match(block, names(data))]))	}
	
	result <- list()
	for (i in (1:length(respvar))) {
		result[[i]] <- list()
		temp.data <- subset(data, subset = (is.na(data[,match(respvar[i], names(data))]) == FALSE))
		cat("\n\nRESPONSE VARIABLE: ", respvar[i], "\n", sep="")
		
		# --- data summary --- #
		funcTrialSum <- class.information2(names(temp.data),temp.data)
		cat("\nDATA SUMMARY: ","\n", sep="")
		print(funcTrialSum)
		cat("\n Number of observations read: ",nrow(temp.data), sep="")
		#cat("\n Number of missing observations: ",nrow(balancedData)-nrow(temp.data), sep="")
    
		# --- test significance of tester:f2lines --- #
		Trt <- temp.data[,match(tester, names(temp.data))]:temp.data[,match(f2lines, names(temp.data))]	
		library(lme4)
		pValue <- 0
		cat("\n\n\nANOVA TABLE: (Trt = ", tester,":", f2lines, ")", sep="")
    
		# --- CONSTRUCT THE MODEL ---#
		if (design == "CRD") {
			myformula1 <- paste(respvar[i], " ~ Trt + (1|",environment,") + (1|Trt:",environment,")",sep = "")
			myformula2 <- paste(respvar[i], " ~ (1|",environment,") + (1|Trt:",environment,")",sep = "")
		}
		if (design == "RCB") {
			myformula1 <- paste(respvar[i], " ~ Trt + (1|",environment,") + (1|",block,":",environment,") + (1|Trt:",environment,")",sep = "")
			myformula2 <- paste(respvar[i], " ~ (1|",environment,") + (1|",block,":",environment,") + (1|Trt:",environment,")",sep = "")
		}	
		
		model1 <- lmer(formula(myformula1), data=temp.data)
		model2 <- lmer(formula(myformula2), data=temp.data)
		anova.table <- anova(model1,model2)
		pValue <- anova.table$"Pr(>Chisq)"[[2]]
		
		# print anova table
		p<-formatC(as.numeric(format(anova.table[2,7], scientific=FALSE)), format="f")
		anova_new<-cbind(round(anova.table[1:6],digits=2), rbind(" ",p))
		rownames(anova_new)<-c(" model2", " model1")
		colnames(anova_new)<-c("Df", "AIC", "BIC", "logLik", "Chisq", "Chi Df","Pr(>Chisq)" )
		anova_print<-replace(anova_new, is.na(anova_new), " ")
		cat("\n Formula for Model 1: ", myformula1, sep="")
		cat("\n Formula for Model 2: ", myformula2,"\n\n", sep="")
		print(anova_print)
		result[[i]]$testerF2.test <- anova_print
		
		# --- check if tester:f2lines is significant. if significant, proceed to test for epistasis --- #
		alpha <- as.numeric(alpha)
		if (pValue < alpha) {
		  cat("\n\nANOVA FOR TESTING EPISTASIS:")
		
		  # --- Create variables: s = P1 + P2 + F1; d = P1 - P2; sd = P1 + P2 - 2F1 --- #
		  if (design == "CRD") {
		    data2 <- reshape(temp.data, v.names=respvar[i],idvar=c(environment,f2lines,rep),timevar=tester,direction="wide",sep=".")
		  }
		  if (design == "RCB") {
		    data2 <- reshape(temp.data, v.names=respvar[i],idvar=c(environment,f2lines,block),timevar=tester,direction="wide",sep=".")
		  }
		  respvardot<-paste(respvar[i],".", sep = "")
		  colnames(data2) <- gsub(respvardot, "", colnames(data2))
		  
		  data2$s <- data2[,match(codeP1, names(data2))] + data2[,match(codeP2, names(data2))] + data2[,match(codeF1, names(data2))]
		  data2$d <- data2[,match(codeP1, names(data2))] - data2[,match(codeP2, names(data2))]
		  data2$sd <- data2[,match(codeP1, names(data2))] + data2[,match(codeP2, names(data2))] - 2*data2[,match(codeF1, names(data2))]
		  
		  m <- nlevels(factor(trim.strings(data2[,match(f2lines, names(data2))])))
		  if (design == "CRD") {
		    r <- nlevels(factor(trim.strings(data2[,match(rep, names(data2))])))
		  }
		  if (design == "RCB") {
		    r <- nlevels(factor(trim.strings(data2[,match(block, names(data2))])))
		  }
		  e <- nlevels(factor(trim.strings(data2[,match(environment, names(data2))])))
		  
		  #--- Test for epistasis ---#
		  library(doBy)
		  if (design == "CRD") {sdrepplusmale<-paste("sd~",rep,"+",f2lines,sep="")}
		  if (design == "RCB") {sdrepplusmale<-paste("sd~",block,"+",f2lines,sep="")}
		  sd.ij <- summaryBy(formula(sdrepplusmale), FUN=sum, data=data2)
		  sdmale<-paste("sd~",f2lines,sep="")
		  sd.i. <- summaryBy(formula(sdmale), FUN=sum, data=data2)
		  sdenvmale<-paste("sd~",environment,"+",f2lines,sep="")
		  sdhi. <- summaryBy(formula(sdenvmale), FUN=sum, data=data2)
		  
		  #--- ANOVA for sd ---#
		  if (design == "CRD") {ANOVA.sd <- ANOVA(data2$sd, data2, cont=6, DESIGN="CRD", TRT=data2[,match(f2lines, names(data2))], Block=NULL, Env=data2[,match(environment, names(data2))], Envname=environment)}
		  if (design == "RCB") {ANOVA.sd <- ANOVA(data2$sd, data2, cont=6, DESIGN="RCB", TRT=data2[,match(f2lines, names(data2))], Block=data2[,match(block, names(data2))], Env=data2[,match(environment, names(data2))], Envname=environment)}
		  
		  SSaxa <- sum(sd.ij$sd.sum)^2/(6*m*r*e)
		  DFaxa <- 1
		  MSaxa <- SSaxa/DFaxa
		  
		  SSad.dd <- ANOVA.sd$SSg
		  DFad.dd <- ANOVA.sd$DFg
		  MSad.dd <- ANOVA.sd$MSg
		  
		  SSaaxe <- ((1/(6*m*r*e))*sum(data2$sd)^2) - SSaxa - ANOVA.sd$SSenv
		  DFaaxe <- e-1
		  MSaaxe <- SSaaxe/DFaaxe
		  
		  SSad.ddxe <- ANOVA.sd$SSgxenv
		  DFad.ddxe <- ANOVA.sd$DFgxenv
		  MSad.ddxe <- ANOVA.sd$MSgxenv
		  
		  SSE <- (1/(6*r))*(sum(sdhi.$sd.sum^2)) - 2*ANOVA.sd$SSenv
		  DFE <- m*e
		  MSE <- SSE/DFE
		  
		  Faxa <- MSaxa/MSaaxe
		  Fad.dd <- MSad.dd/MSad.ddxe
		  Faaxe <- MSaaxe/ANOVA.sd$MSe
		  Fad.ddxe <- MSad.ddxe/ANOVA.sd$MSe
		  FE <- MSE/ANOVA.sd$MSe
		  Fe_sd <- NA
		  
		  Paxa <- 1-pf(Faxa,DFaxa,DFaaxe)
		  Pad.dd <- 1-pf(Fad.dd,DFad.dd,DFad.ddxe)
		  Paaxe <- 1-pf(Faaxe,DFaaxe,ANOVA.sd$DFe)
		  Pad.ddxe <- 1-pf(Fad.dd,DFad.dd,ANOVA.sd$DFe)
		  
		  PE <- 1-pf(FE,DFE,ANOVA.sd$DFe)
		  Pe_sd <- NA
		  
		  SV <- c("SV","AxA", "AxD and DxD", "AAxE", "ADxE or DDxE", "Total", "Residual")
		  DF <- format(rbind("Df",format(round(rbind(DFaxa, DFad.dd, DFaaxe, DFad.ddxe, DFE, ANOVA.sd$DFe),digits=0),nsmall=0)),justify="right")
		  SS <- format(rbind("Sum Sq",format(round(rbind(SSaxa, SSad.dd, SSaaxe, SSad.ddxe, SSE, ANOVA.sd$SSe),digits=4),nsmall=4)),justify="right")
		  MS <- format(rbind("Mean Sq",format(round(rbind(MSaxa, MSad.dd, MSaaxe, MSad.ddxe, MSE, ANOVA.sd$MSe),digits=4),nsmall=4)),justify="right")
		  Fvalue <- format(rbind("F value",format(round(rbind(Faxa, Fad.dd, Faaxe, Fad.ddxe, FE),digits=2),nsmall=2),""),justify="right")
		  Paxa2<-formatC(as.numeric(format(Paxa, scientific=FALSE)), format="f")
		  Pad.dd2<-formatC(as.numeric(format(Pad.dd, scientific=FALSE)), format="f")
		  Paaxe2<-formatC(as.numeric(format(Paaxe, scientific=FALSE)), format="f")
		  Pad.ddxe2<-formatC(as.numeric(format(Pad.ddxe, scientific=FALSE)), format="f")
		  PE2<-formatC(as.numeric(format(PE, scientific=FALSE)), format="f")
		  Prob <-format(rbind("Pr(>F)",Paxa2, Pad.dd2, Paaxe2, Pad.ddxe2, PE2, ""),justify="right")
		  
		  AOV <- as.table(cbind(SV, DF, SS, MS, Fvalue, Prob))
		  colnames(AOV) <- rep("",ncol(AOV))
		  rownames(AOV) <- rep("",nrow(AOV))
		  print(AOV)
		  result[[i]]$epistasis.test<- AOV
		  # --- NOTE: There is epistatis if either the AxA or AxD and DxD or both are significant. --- #
		  
		  #--- estimates of variance components --- #
		  
		  #---  ANOVA FOR s ---#
		  if (design == "CRD") {ANOVA.s <- ANOVA(data2$s, data2, cont=3, DESIGN="CRD", TRT=data2[,match(f2lines, names(data2))], Block=NULL, Env=data2[,match(environment, names(data2))], Envname=environment)}
		  if (design == "RCB") {ANOVA.s <- ANOVA(data2$s, data2, cont=3, DESIGN="RCB", TRT=data2[,match(f2lines, names(data2))], Block=data2[,match(block, names(data2))], Env=data2[,match(environment, names(data2))], Envname=environment)}
		  
		  Fg_s <- ANOVA.s$MSg/ANOVA.s$MSe
		  Fgxe_s <- ANOVA.s$MSgxenv/ANOVA.s$MSe
		  Fe_s <- NA
		  
		  Pg_s <- 1-(pf(Fg_s,ANOVA.s$DFg,ANOVA.s$DFe))
		  Pgxe_s <- 1-(pf(Fgxe_s,ANOVA.s$DFgxenv,ANOVA.s$DFe))
		  Pe_s <- NA
		  
		  #---  ANOVA FOR d ---#
		  if (design == "CRD") {ANOVA.d <- ANOVA(data2$d, data2, cont=2, DESIGN="CRD", TRT=data2[,match(f2lines, names(data2))],Block=NULL, Env=data2[,match(environment, names(data2))], Envname=environment)}
		  if (design == "RCB") {ANOVA.d <- ANOVA(data2$d, data2, cont=2, DESIGN="RCB", TRT=data2[,match(f2lines, names(data2))],Block=data2[,match(block, names(data2))], Env=data2[,match(environment, names(data2))],Envname=environment)}
		  
		  Fg_d <- ANOVA.d$MSg/ANOVA.d$MSe
		  Fgxe_d <- ANOVA.d$MSgxenv/ANOVA.d$MSe
		  Fe_d <- NA
		  
		  Pg_d <- 1-(pf(Fg_d,ANOVA.d$DFg,ANOVA.d$DFe))
		  Pgxe_d <- 1-(pf(Fgxe_d,ANOVA.d$DFgxenv,ANOVA.d$DFe))
		  Pe_d <- NA
		  
		  SV <- c("SV","s", "sxe", "e(s)", "d", "dxe", "e(d)")
		  DF <- format(rbind("Df",format(round(rbind(ANOVA.s$DFg, ANOVA.s$DFgxe, ANOVA.s$DFe, ANOVA.d$DFg, ANOVA.d$DFgxe, ANOVA.d$DFe),digits=0),nsmall=0)),justify="right")
		  SS <- format(rbind("Sum Sq",format(round(rbind(ANOVA.s$SSg, ANOVA.s$SSgxe, ANOVA.s$SSe, ANOVA.d$SSg, ANOVA.d$SSgxe, ANOVA.d$SSe),digits=4),nsmall=4)),justify="right")
		  MS <- format(rbind("Mean Sq",format(round(rbind(ANOVA.s$MSg, ANOVA.s$MSgxe, ANOVA.s$MSe, ANOVA.d$MSg, ANOVA.d$MSgxe, ANOVA.d$MSe),digits=4),nsmall=4)),justify="right")
		  Fvalue <- format(rbind("F value",format(round(Fg_s,digits=2),nsmall=2), format(round(Fgxe_s,digits=2),nsmall=2),"", format(round(Fg_d,digits=2),nsmall=2), format(round(Fgxe_d,digits=2),nsmall=2),""),justify="right")
		  Pg_s2 <- formatC(as.numeric(format(Pg_s,scientific=FALSE)),format="f")
		  Pgxe_s2 <- formatC(as.numeric(format(Pgxe_s,scientific=FALSE)),format="f")
		  Pg_d2 <- formatC(as.numeric(format(Pg_d,scientific=FALSE)),format="f")
		  Pgxe_d2 <- formatC(as.numeric(format(Pgxe_d,scientific=FALSE)),format="f")
		  Prob <-format(rbind("Pr(>F)",Pg_s2, Pgxe_s2, "", Pg_d2, Pgxe_d2, ""),justify="right")
		  
		  AOV2 <- as.table(cbind(SV, DF, SS, MS, Fvalue, Prob))
		  colnames(AOV2) <- rep("",ncol(AOV2))
		  rownames(AOV2) <- rep("",nrow(AOV2))
		  cat("\n\nANOVA TABLE")
      print(noquote(AOV2))
		  result[[i]]$sde.anova <-noquote(AOV2)
		  
		  sigma2A <- (ANOVA.s$MSg - ANOVA.s$MSe)/(3*r)
		  sigma2D <- (ANOVA.d$MSg - ANOVA.d$MSe)/(2*r)
		  
		  #--- Genetic variance estimates ---#
		  varcomp <- summary(model1)@REmat
		  Ve <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
		  VA <- 4*sigma2A
		  VVA <- (32/(9*(r^2)))*(((ANOVA.s$MSg^2)/(m-1+2))+((ANOVA.s$MSe^2)/((m-1)*(r-1)+2)))
		  
		  VD <- 2*sigma2D
		  if (VD < 0) VD <- 0
		  VVD <- (2/(r^2))*(((ANOVA.d$MSg^2)/(m-1+2)) + ((ANOVA.d$MSe^2)/((m-1)*(r-1)+2)))
		  
		  VE <- Ve
		  VP <- VA + VD + VE
		  
		  h2B <- (VA + VD) / VP                 # individual based
		  # Vh2B <- (2*(1-h2B)^2*(1+(r-1)*h2B)^2)/(r*(r-1)*(m-1))
		  
		  h2N <- VA / VP      	              # individual based
		  # Vh2N <- (2*(1-h2N)^2*(1+(r-1)*h2N)^2)/(r*(r-1)*(m-1))
		  Dominance.ratio <- sqrt(2*VD/VA)
		  
		  Estimate <- round(rbind(VA, VD, h2N, h2B, Dominance.ratio), digits=4)
		  Estimate <- as.table(Estimate)
		  colnames(Estimate ) <- c("Estimate")
		  rownames(Estimate ) <- c(" VA", " VD", " h2-narrow sense"," H2-broad sense", " Dominance Ratio")
		  cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:\n")
      print(Estimate)
		  result[[i]]$genvar.components <- Estimate
    }
		cat("\n")
		cat("\n==============================================================\n")
	}## end of loop (i)
	detach("package:lme4")
	return(list(output = result))
}

