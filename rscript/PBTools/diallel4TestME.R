diallel4TestME <-
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
  environment <-trim.strings(environment)
  alpha <- trim.strings(alpha)
	
  # --- create titles --- #
  if (cross) {parentsType<-"CROSS"
  } else {parentsType<-"SELF"}
  cat("\nMULTIPLE ENVIRONMENT ANALYSIS\n")
  cat("\nDIALLEL ANALYSIS: GRIFFING METHOD IV IN ",design, " (", parentsType, ")\n", sep="")
  
	#define factors
	data[,match(p1, names(data))] <- factor(trim.strings(data[,match(p1, names(data))]))
	data[,match(p2, names(data))] <- factor(trim.strings(data[,match(p2, names(data))]))
	data[,match(environment, names(data))] <- factor(trim.strings(data[,match(environment, names(data))]))
	if (design == "CRD") {data[,match(rep, names(data))] <- factor(trim.strings(data[,match(rep, names(data))])) }
	if (design == "RCB") {data[,match(block, names(data))] <- factor(trim.strings(data[,match(block, names(data))]))	}
	
	result <- list()
	for (z in (1:length(respvar))) {
		result[[z]] <- list()
		cat("\n\nRESPONSE VARIABLE: ", respvar[z], "\n", sep="")
		
		if (design == "CRD") {
			vec<-c(match(p1, names(data)), match(p2, names(data)), match(environment, names(data)), match(rep, names(data)), match(respvar[z], names(data)))
			temp2.data <-subset(data, select=vec)  
		}
		if (design == "RCB") {
			vec<-c(match(p1, names(data)), match(p2, names(data)), match(environment, names(data)), match(block, names(data)), match(respvar[z], names(data)))
			temp2.data <-subset(data, select=vec)
		}
		temp.data <- subset(temp2.data, subset = (is.na(data[,match(respvar[z], names(data))]) == FALSE))
		
		# --- data summary --- #
		funcTrialSum <- class.information2(names(temp.data),temp.data)
		cat("\nDATA SUMMARY: ","\n", sep="")
		print(funcTrialSum)
		cat("\n Number of observations read: ",nrow(temp.data), sep="")
		#cat("\n Number of missing observations: ",nrow(balancedData)-nrow(temp.data), sep="")
    
		name1 <- levels(temp.data[,match(p1, names(temp.data))])
		name2 <- levels(temp.data[,match(p2, names(temp.data))])
		name <- unique(c(name1, name2))
		
		p <- nlevels(factor(name))
		e <- nlevels(temp.data[,match(environment, names(temp.data))])
		if (design == "CRD") {r <- nlevels(temp.data[,match(rep, names(temp.data))]) }
		if (design == "RCB") {r <- nlevels(temp.data[,match(block, names(temp.data))]) }
		
		# --- constructing the ANOVA table ---#
		if (design == "CRD") {myformula1 <- paste(respvar[z], " ~ ", environment, " + ", p1, ":", p2, " + ", environment, ":", p1, ":", p2, sep = "") }
		if (design == "RCB") {myformula1 <- paste(respvar[z], " ~ ", environment, " + ", block, ":", environment," + ", p1, ":", p2, " + ", environment, ":", p1, ":", p2, sep = "") }
		model <- aov(formula(myformula1), data = temp.data)
		out1 <- anova(model)
		
		EnvDF <- out1[rownames(out1)==environment,"Df"]
		EnvSS <- out1[rownames(out1)==environment,"Sum Sq"]
		EnvMS <- out1[rownames(out1)==environment,"Mean Sq"]
		
		p1p2 <- paste(p1,":",p2,sep="")
		HybridDF <- out1[rownames(out1)==p1p2,"Df"]
		HybridSS <- out1[rownames(out1)==p1p2,"Sum Sq"]
		HybridMS <- out1[rownames(out1)==p1p2,"Mean Sq"]
		
		envp1p2 <- paste(environment, ":", p1, ":", p2,sep="")
		HybridxEDF <- out1[rownames(out1)==envp1p2,"Df"]
		HybridxESS <- out1[rownames(out1)==envp1p2,"Sum Sq"]
		HybridxEMS <- out1[rownames(out1)==envp1p2,"Mean Sq"]
		
		EDF <- out1[rownames(out1)=="Residuals","Df"]
		ESS <- out1[rownames(out1)=="Residuals","Sum Sq"]
		EMS <- out1[rownames(out1)=="Residuals","Mean Sq"]
		
		F.Hybrid <- HybridMS/HybridxEMS
		P.Hybrid <- 1-pf(F.Hybrid, HybridDF, HybridxEDF)
		F.HybridxE <- out1[rownames(out1)==envp1p2,"F value"]
		P.HybridxE <- out1[rownames(out1)==envp1p2,"Pr(>F)"]
		F.error <- NA
		P.error <- NA
		
		if (design == "RCB") {
			envblock <- paste(environment, ":", block,sep="")
			repEnvDF <- out1[rownames(out1)==envblock,"Df"]
			repEnvSS <- out1[rownames(out1)==envblock,"Sum Sq"]
			repEnvMS <- out1[rownames(out1)==envblock,"Mean Sq"]
			F.repEnv <- repEnvMS/EMS
			P.repEnv <- 1-pf(F.repEnv, repEnvDF, EDF)
			
			F.Env <- EnvMS/repEnvMS
			P.Env <- 1-pf(F.Env, EnvDF, repEnvDF)
			
			blockenv <- paste(block, "(", environment, ")",sep="")
			SV <- c("SV",environment, blockenv, "Hybrid", "HybridxE", "Error")
			DF <- format(rbind("Df",EnvDF,repEnvDF,HybridDF,HybridxEDF, EDF), justify="right")
			SSq <- format(rbind("Sum Sq",format(round(rbind(EnvSS,repEnvSS,HybridSS,HybridxESS, ESS),digits=2))), justify="right")
			MSq <- format(rbind("Mean Sq",format(round(rbind(EnvMS,repEnvMS,HybridMS,HybridxEMS, EMS),digits=2))), justify="right")
			Fvalue <- format(round(rbind(F.Env, F.repEnv, F.Hybrid, F.HybridxE),digits=2))
			Fvalue2 <- format(rbind("F value", Fvalue, " "), justify="right")
			P.Env2<-formatC(as.numeric(format(P.Env,scientific=FALSE)),format="f")
			P.repEnv2<-formatC(as.numeric(format(P.repEnv,scientific=FALSE)),format="f")
			P.Hybrid2<-formatC(as.numeric(format(P.Hybrid,scientific=FALSE)),format="f")
			P.HybridxE2<-formatC(as.numeric(format(P.HybridxE,scientific=FALSE)),format="f")
			P <- format(rbind("Pr(>F)",P.Env2,P.repEnv2,P.Hybrid2,P.HybridxE2, " "), justify="right")
		}
		
		if (design == "CRD") {
			F.Env <- out1[rownames(out1)==environment,"F value"]
			P.Env <- out1[rownames(out1)==environment,"Pr(>F)"]
			
			SV <- c("SV",environment, "Hybrid", "HybridxE", "Error")
			DF <- format(rbind("Df",EnvDF,HybridDF,HybridxEDF, EDF), justify="right")
			SSq <-  format(rbind("Sum Sq",format(round(rbind(EnvSS,HybridSS,HybridxESS, ESS),digits=2))), justify="right")
			MSq <- format(rbind("Mean Sq",format(round(rbind(EnvMS,HybridMS,HybridxEMS, EMS),digits=2))), justify="right")
			Fvalue <- format(round(rbind(F.Env, F.Hybrid, F.HybridxE),digits=2))
			Fvalue2 <- format(rbind("F value", Fvalue, " "), justify="right")
			P.Env2<-formatC(as.numeric(format(P.Env,scientific=FALSE)),format="f")
			P.Hybrid2<-formatC(as.numeric(format(P.Hybrid,scientific=FALSE)),format="f")
			P.HybridxE2<-formatC(as.numeric(format(P.HybridxE,scientific=FALSE)),format="f")
			P <- format(rbind("Pr(>F)",P.Env2,P.Hybrid2,P.HybridxE2, " "), justify="right")
		}
		
		AOV <- noquote(cbind(SV,DF,SSq,MSq,Fvalue2,P))
		colnames(AOV) <- c("", "", "", "", "", "")
		rownames(AOV) <- rep("",nrow(AOV))
		cat("\n\n\nANOVA TABLE")
		print(AOV)
		result[[z]]$anova.table <- AOV
		
		#--- CREATION OF ESTIMABLE EFFECTS ---#
		i <- rep(1:p, each=p)
		j <- rep(1:p,p)
		ID <- cbind(i,j)
		
		#--- for method 4 ---#
		M4 <- subset(ID, i < j)
		t <- nrow(M4)
		
		#--- GCA estimable effects ---#
		GCA <- as.matrix(rep(0,t*(p-1)),nrow=t, ncol=(p-1))
		dim(GCA) <- c(t,(p-1))
		for (o in 1:t)
			for (n in 1:(p-1))  { 
				GCA[o,n] <- (as.logical(M4[o,"i"]==n)-
							as.logical(M4[o,"i"]==p))+
						(as.logical(M4[o,"j"]==n)-
							as.logical(M4[o,"j"]==p))  
			}
		
		colnames(GCA) <- paste("GCA", 1:(p-1), sep="")
		
		#--- SCA estimable effects ---#
		sca.col <- as.data.frame(subset(ID, (i<=(p-3) & (j >= (i+1) & j<=(p-1)))))
		ncol <- nrow(sca.col)
		
		SCA <- as.matrix(rep(0,t*ncol),nrow=t, ncol=ncol)
		dim(SCA) <- c(t,ncol)
		
		for (o in 1:t)
			for (n in 1:ncol)  {
				SCA[o,n] <- as.logical(M4[o,"i"]==sca.col[n,"i"])*
						as.logical(M4[o,"j"]==sca.col[n,"j"])-
						(as.logical(M4[o,"i"]==sca.col[n,"i"])+
							as.logical(M4[o,"i"]==sca.col[n,"j"]))*
						as.logical(M4[o,"j"]==p);
				if (((M4[o,"i"]>=(p-2))&(M4[o,"j"]>=(p-1)))|((M4[o,"i"]>=(p-1))&(M4[o,"j"]>=(p-2))))  {
					SCA[o,n] <- -1*as.logical(M4[o,"i"]==(p-2))*
							as.logical(M4[o,"j"]==(p-1))+
							as.logical(M4[o,"i"]>=(p-2))*
							as.logical(M4[o,"j"]==p)*
							as.logical(M4[o,"i"] != sca.col[n,"j"])  
				}
			}
		colnames(SCA) <- paste("SCA",sca.col$"i", sca.col$j, sep="")
		
		esteffect <- as.data.frame(cbind(M4,GCA,SCA))
		data.all <- merge(temp.data,esteffect, by.x=c(p1,p2), by.y=c("i","j"))
		
		if (design == "CRD") {myformula2 <- paste(respvar[z], " ~ ", environment, " + ", rep, ":", environment," + . - ", p1, " - ", p2, " - ", rep, " + (. - ", p1, " - ", p2, " - ", rep, "):",environment, sep = "") }
		if (design == "RCB") {myformula2 <- paste(respvar[z], " ~ ", environment, " + ", block, ":", environment," + . - ", p1, " - ", p2, " - ", block, " + (. - ", p1, " - ", p2, " - ", block, "):",environment, sep = "") }
		model2 <- aov(formula(myformula2), data = data.all)
		out2 <- anova(model2)
		
		GCA.effect <- subset(out2, substr(rownames(out2),1,3)=="GCA")
		GCA.df <- sum(GCA.effect$Df)
		GCA.SS <- sum(GCA.effect$"Sum Sq")
		GCA.MS <- GCA.SS/GCA.df
		
		SCA.effect <- subset(out2, substr(rownames(out2),1,3)=="SCA")
		SCA.df <- sum(SCA.effect$Df)
		SCA.SS <- sum(SCA.effect$"Sum Sq")
		SCA.MS <- SCA.SS/SCA.df
		
		envgca <- paste(environment, ":GCA",sep="")
		GCAxE.effect <- subset(out2, substr(rownames(out2),1,7)==envgca)
		GCAxE.df <- sum(GCAxE.effect$Df)
		GCAxE.SS <- sum(GCAxE.effect$"Sum Sq")
		GCAxE.MS <- GCAxE.SS/GCAxE.df
		
		envsca <- paste(environment, ":SCA",sep="")
		SCAxE.effect <- subset(out2, substr(rownames(out2),1,7)==envsca)
		SCAxE.df <- sum(SCAxE.effect$Df)
		SCAxE.SS <- sum(SCAxE.effect$"Sum Sq")
		SCAxE.MS <- SCAxE.SS/SCAxE.df
		
		GCA.F <- GCA.MS/GCAxE.MS
		SCA.F <- SCA.MS/SCAxE.MS
		
		GCAxE.F <- GCAxE.MS/EMS
		SCAxE.F <- SCAxE.MS/EMS
		EF <- NA
		
		GCA.P <- 1 - pf(GCA.F, GCA.df, GCAxE.df)
		SCA.P <- 1 - pf(SCA.F, SCA.df, SCAxE.df)
		
		GCAxE.P <- 1 - pf(GCAxE.F, GCAxE.df, EDF)
		SCAxE.P <- 1 - pf(SCAxE.F, SCAxE.df, EDF)
		EP <- NA
		
		SV <- c("SV","GCA", "SCA", "GCAxE", "SCAxE", "Error")
		DF <- format(rbind("Df",GCA.df, SCA.df, GCAxE.df, SCAxE.df, EDF),justify="right")
		SSq <- format(rbind("Sum Sq",format(round(rbind(GCA.SS, SCA.SS, GCAxE.SS, SCAxE.SS, ESS),digits=2), nsmall=2)), justify="right")
		MSq <- format(rbind("Mean Sq",format(round(rbind(GCA.MS, SCA.MS, GCAxE.MS, SCAxE.MS, EMS),digits=2), nsmall=2)), justify="right")
		Fvalue <- format(rbind("F value",format(round(rbind(GCA.F, SCA.F, GCAxE.F, SCAxE.F),digits=2),nsmall=2)," "), justify="right")
		GCA.p2<-formatC(as.numeric(format(GCA.P,scientific=FALSE)),format="f")
		SCA.p2<-formatC(as.numeric(format(SCA.P,scientific=FALSE)),format="f")
		GCAxE.p2<-formatC(as.numeric(format(GCAxE.P,scientific=FALSE)),format="f")
		SCAxE.p2<-formatC(as.numeric(format(SCAxE.P,scientific=FALSE)),format="f")
		P <- format(rbind("Pr(>F)",GCA.p2, SCA.p2, GCAxE.p2, SCAxE.p2, " "), justify="right")
		
		AOV2 <- as.table(cbind (SV,DF,SSq,MSq,Fvalue,P))
		colnames(AOV2) <- c("", "", "", "", "", "")
		rownames(AOV2) <- rep("",nrow(AOV2))
		cat("\n\nANOVA TABLE")
		print(AOV2)
		result[[z]]$anova2.table <- AOV2
		
		# --- Estimates of GCA, SCA --- #
    
		library(doBy)
		myformula4<- paste(respvar[z], " ~ ", p1," + ",p2, sep = "")
		meandata <- summaryBy(formula(myformula4), data=temp.data)
		
		# serial to parallel of meandata
		mdata <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
		dim(mdata) <- c(p,p)
		
		for (I in 1:p)  {
			for (J in 1:p)   {
				if (I < J) mdata[I,J] <- meandata[(meandata[,match(p1, names(meandata))]==I & meandata[,match(p2, names(meandata))]==J),3]
			} 
		}
		
		# printing the matrix of means
		mdata2 <- format(round(mdata, 4),nsmall=4)
		rownames(mdata2) <- colnames(mdata2) <- levels(factor(name))
		mdata3 <-noquote(format(gsub(" 0.0000", "", mdata2),justify="right"))
		cat("\n\nMATRIX OF MEANS\n")
		print(mdata3)
		result[[z]]$means.matrix <- mdata3
		
		mirror <- mdata
		
		for (I in 1:p)  {
			for (J in 1:p)   {
				if (I > J) mirror[I,J] <- mirror[J,I]
			} 
		}
		
		XI <- colSums(mirror)
		SUMX <- sum(mdata)
		
		#----- Estimates of GCA Effects -----#
		G_SCA <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
		dim(G_SCA) <- c(p,p)
		for (I in 1:p)  {
			for (J in 1:p)   {
				if (I == J) G_SCA[I,J] <- (1/(p*(p-2))) * (p*XI[I]-(2*SUMX))
			} 
		}
		
		#----- Estimates of SCA Effects -----#
		for (I in 1:p)  {
			for (J in 1:p)   {
				if (I < J)   {
					B1 <- (1/(p-2))* (XI[I] + XI[J])
					B2 <- (2/((p-1)*(p-2)))*SUMX
					G_SCA[I,J] <- mdata[I,J] - B1 + B2
				}
			} 
		}
		colnames(G_SCA) <- rownames(G_SCA) <- levels(factor(name))

		G_SCA2 <-format(round(G_SCA, digits=4), nsmall=4)
		G_SCA3 <-noquote(format(gsub(" 0.0000", "", G_SCA2),justify="right"))
		cat("\n\nGENERAL COMBINING ABILITY EFFECTS (diagonal), SPECIFIC COMBINING\n")
		cat("ABILITY EFFECTS (above diagonal):\n")
		print(G_SCA3)
		result[[z]]$gcasca.matrix <- G_SCA3
		
		# --- estimates of standard errors --- #
		
		MEPRIME <-EMS/r      
		Ve <- MEPRIME
		SE_GI <- sqrt((p-1)*Ve/(p*e*(p-2)))
		LSD_GI <- NA
		SE_SIJ <- sqrt((p-3)*Ve/(e*(p-1)))
		LSD_SIJ <- NA
		
		SE_Gdiff <- sqrt(2*Ve/(e*(p-2)))
		LSD_Gdiff <- qt(.975,EDF)*SE_Gdiff
		SE_SIJ_SIK <- sqrt(2*(p-3)*Ve/(e*(p-2)))
		LSD_SIJ_SIK <- qt(.975,EDF)*SE_SIJ_SIK
		SE_SIJ_SKL <- sqrt(2*(p-4)*Ve/(e*(p-2)))
		LSD_SIJ_SKL <- qt(.975,EDF)*SE_SIJ_SKL
		
		STDERR <- round(rbind(SE_GI,SE_SIJ, SE_Gdiff, SE_SIJ_SIK, SE_SIJ_SKL), digits=4)
		LSD <- round(rbind(LSD_GI,LSD_SIJ, LSD_Gdiff, LSD_SIJ_SIK, LSD_SIJ_SKL),digits=4)
		VAREST <- as.table(cbind(STDERR, LSD))
		
		rownames(VAREST) <- c(' Gi', ' Sij', ' Gi-Gj', ' Sij-Sik', ' Sij-Skl')
		colnames(VAREST) <- c("Std. Error", "LSD")
		cat("\n\nTABLE OF STANDARD ERRORS AND LSDs:\n")
		print(VAREST)
		result[[z]]$stderror.table <- VAREST
		
		#--- Variance Component ---#
		
		Vg <- (GCA.MS-GCAxE.MS-SCA.MS+SCAxE.MS)/(r*e*(p-2))
		Vs <- (SCA.MS-SCAxE.MS)/(r*e)
		Vge <- (GCAxE.MS-SCAxE.MS)/(r*(p-2))
		Vse <- (SCAxE.MS-EMS)/r
		
		VC <- round(rbind(Vg, Vs, Vge, Vse, EMS), digits=4)
		colnames(VC) <- c("Estimate")
		rownames(VC) <- c(" GCA", " SCA", " GCAxE", " SCAxE", " Error")
		TABLE <- as.table(VC)
		cat("\n\nESTIMATES OF VARIANCE COMPONENTS:\n")
		print(TABLE)
		result[[z]]$var.components <-TABLE
		
		if (Vg<0) Vg <- 0
		if (Vs<0) Vs <- 0
		if (Vge<0) Vge <- 0
		if (Vse<0) Vse <- 0
		
		#---- Genetic Variance components ----#
		if (cross) {F<-0}
		else {F<-1}
    
    VA <- (4/(1+F))*Vg
		VD <- (4/(1+F)^2)*Vs
		if (VD < 0) VD <- 0
		
		VAE <- (4/(1+F))*Vge
		VDE <- (4/(1+F)^2)*Vse
		if (VDE < 0) VDE <- 0
		
		VE <- EMS
		VP <- VA + VD + VAE + VDE+ VE
		h2B <- (VA + VD) / VP               # individual based
		h2N <- VA / VP                 	# individual based
		Dominance.ratio <- sqrt(2*VD/VA)   
		
		Estimate <- round(rbind(VA, VD, VAE, VDE, h2N, h2B, Dominance.ratio), digits=4)
		colnames(Estimate) <- c("Estimate")
		rownames(Estimate) <- c(" VA", " VD", " VAE", " VDE", " h2-narrow sense", " H2-broad sense", " Dominance Ratio")
		TABLE2 <- as.table(Estimate)
		cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:\n")
		print(TABLE2)
		result[[z]]$genvar.components <-TABLE2
		cat("\n")
		cat("\n==============================================================\n")
	}## end of loop (z)
  detach("package:doBy")
  return(list(output = result))
}

