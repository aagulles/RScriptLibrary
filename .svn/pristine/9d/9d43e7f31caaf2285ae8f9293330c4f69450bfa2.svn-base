diallel3TestME <-
function(design = c("CRD", "RCB"), data, respvar, p1, p2, rep=NULL, block=NULL, cross=TRUE, individual=NULL, environment=NULL, alpha=0.05) {
  
  options(show.signif.stars=FALSE)
  data <- eval(parse(text = data))
	
  # --- trim the strings --- 
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
  cat("\nDIALLEL ANALYSIS: GRIFFING METHOD III IN ",design, " (", parentsType, ")\n", sep="")
  
	# --- define factors --- #
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
    
		p <- nlevels(temp.data[,match(p1, names(temp.data))])
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
		
		# --- GCA, SCA, REC, GCAxE, SCAxE, and RECxE; df, sums of squares, and mean squares --- #
		
		esteffect <- as.data.frame(cbind(M3,GCA,SCA,REC))
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
		
		REC.effect <- subset(out2, substr(rownames(out2),1,3)=="REC")
		REC.df <- sum(REC.effect$Df)
		REC.SS <- sum(REC.effect$"Sum Sq")
		REC.MS <- REC.SS/REC.df
		
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
		
		envrec <- paste(environment, ":REC",sep="")
		RECxE.effect <- subset(out2, substr(rownames(out2),1,7)==envrec)
		RECxE.df <- sum(RECxE.effect$Df)
		RECxE.SS <- sum(RECxE.effect$"Sum Sq")
		RECxE.MS <- RECxE.SS/RECxE.df
		
		# --- MAT, NONM, MATxE, and NONMCxE; df, sums of squares, and mean squares --- #
		
		esteffect <- as.data.frame(cbind(M3,GCA,SCA,MAT, NONM))
		data.all <- merge(temp.data,esteffect, by.x=c(p1,p2), by.y=c("i","j"))
		
		if (design == "CRD") {myformula3 <- paste(respvar[z], " ~ ", environment, " + ", rep, ":", environment," + . - ", p1, " - ", p2, " - ", rep, " + (. - ", p1, " - ", p2, " - ", rep, "):",environment, sep = "") }
		if (design == "RCB") {myformula3 <- paste(respvar[z], " ~ ", environment, " + ", block, ":", environment," + . - ", p1, " - ", p2, " - ", block, " + (. - ", p1, " - ", p2, " - ", block, "):",environment, sep = "") }
		model3 <- aov(formula(myformula3), data = data.all)
		out3 <- anova(model3)
		
		MAT.effect <- subset(out3, substr(rownames(out3),1,3)=="MAT")
		MAT.df <- sum(MAT.effect$Df)
		MAT.SS <- sum(MAT.effect$"Sum Sq")
		MAT.MS <- MAT.SS/MAT.df
		
		NONM.effect <- subset(out3, substr(rownames(out3),1,4)=="NONM")
		NONM.df <- sum(NONM.effect$Df)
		NONM.SS <- sum(NONM.effect$"Sum Sq")
		NONM.MS <- NONM.SS/NONM.df
		
		envmat <- paste(environment, ":MAT",sep="")
		MATxE.effect <- subset(out3, substr(rownames(out3),1,7)==envmat)
		MATxE.df <- sum(MATxE.effect$Df)
		MATxE.SS <- sum(MATxE.effect$"Sum Sq")
		MATxE.MS <- MATxE.SS/MATxE.df
		
		envnonm <- paste(environment, ":NONM",sep="")
		NONMxE.effect <- subset(out3, substr(rownames(out3),1,8)==envnonm)
		NONMxE.df <- sum(NONMxE.effect$Df)
		NONMxE.SS <- sum(NONMxE.effect$"Sum Sq")
		NONMxE.MS <- NONMxE.SS/NONMxE.df
		
		# --- computation of F and Pvalue --- #
			
		GCA.F <- GCA.MS/GCAxE.MS
		SCA.F <- SCA.MS/SCAxE.MS
		REC.F <- REC.MS/RECxE.MS
		MAT.F <- MAT.MS/MATxE.MS
		NONM.F <- NONM.MS/NONMxE.MS
		
		GCAxE.F <- GCAxE.MS/EMS
		SCAxE.F <- SCAxE.MS/EMS
		RECxE.F <- RECxE.MS/EMS
		MATxE.F <- MATxE.MS/EMS
		NONMxE.F <- NONMxE.MS/EMS
		EF <- NA
		
		GCA.P <- 1 - pf(GCA.F, GCA.df, GCAxE.df)
		SCA.P <- 1 - pf(SCA.F, SCA.df, SCAxE.df)
		REC.P <- 1 - pf(REC.F, REC.df, RECxE.df)
		MAT.P <- 1 - pf(MAT.F, MAT.df, MATxE.df)
		NONM.P <- 1 - pf(NONM.F, NONM.df, NONMxE.df)
		
		GCAxE.P <- 1 - pf(GCAxE.F, GCAxE.df, EDF)
		SCAxE.P <- 1 - pf(SCAxE.F, SCAxE.df, EDF)
		RECxE.P <- 1 - pf(RECxE.F, RECxE.df, EDF)
		MATxE.P <- 1 - pf(MATxE.F, MATxE.df, EDF)
		NONMxE.P <- 1 - pf(NONMxE.F, NONMxE.df, EDF)
		EP <- NA
		
		SV <- c("SV","GCA", "SCA", "REC", "  MAT", "  NONM", "GCAxE", "SCAxE", "RECxE", "  MATxE", "  NONMxE", "Error")
		DF <- format(rbind("Df", paste(GCA.df), paste(SCA.df), paste(REC.df), paste(MAT.df), paste(NONM.df), paste(GCAxE.df), paste(SCAxE.df), paste(RECxE.df), paste(MATxE.df), paste(NONMxE.df),paste(EDF)),justify="right")
		SSq <- format(rbind("Sum Sq",format(round(rbind(GCA.SS,SCA.SS,REC.SS,MAT.SS,NONM.SS,GCAxE.SS,SCAxE.SS,RECxE.SS,MATxE.SS,NONMxE.SS,ESS),digits=2), nsmall=2)), justify="right")
		MSq <- format(rbind("Mean Sq",format(round(rbind(GCA.MS, SCA.MS, REC.MS, MAT.MS, NONM.MS, GCAxE.MS, SCAxE.MS, RECxE.MS, MATxE.MS, NONMxE.MS, EMS),digits=2), nsmall=2)), justify="right")
		Fvalue <- format(rbind("F value",format(round(rbind(GCA.F, SCA.F, REC.F, MAT.F, NONM.F, GCAxE.F, SCAxE.F, RECxE.F, MATxE.F, NONMxE.F),digits=2),nsmall=2)," "), justify="right")
		GCA.p2<-formatC(as.numeric(format(GCA.P,scientific=FALSE)),format="f")
		SCA.p2<-formatC(as.numeric(format(SCA.P,scientific=FALSE)),format="f")
		REC.p2<-formatC(as.numeric(format(REC.P,scientific=FALSE)),format="f")
		MAT.p2<-formatC(as.numeric(format(MAT.P,scientific=FALSE)),format="f")
		NONM.p2<-formatC(as.numeric(format(NONM.P,scientific=FALSE)),format="f")
		GCAxE.p2<-formatC(as.numeric(format(GCAxE.P,scientific=FALSE)),format="f")
		SCAxE.p2<-formatC(as.numeric(format(SCAxE.P,scientific=FALSE)),format="f")
		RECxE.p2<-formatC(as.numeric(format(RECxE.P,scientific=FALSE)),format="f")
		MATxE.p2<-formatC(as.numeric(format(MATxE.P,scientific=FALSE)),format="f")
		NONMxE.p2<-formatC(as.numeric(format(NONMxE.P,scientific=FALSE)),format="f")
		P <- format(rbind("Pr(>F)",GCA.p2, SCA.p2, REC.p2, MAT.p2, NONM.p2, GCAxE.p2, SCAxE.p2, RECxE.p2, MATxE.p2, NONMxE.p2, " "), justify="right")
		
		AOV2 <- as.table(cbind (SV,DF,SSq,MSq,Fvalue,P))
		colnames(AOV2) <- c("", "", "", "", "", "")
		rownames(AOV2) <- rep("",nrow(AOV2))
		cat("\n\nANOVA TABLE")
		print(AOV2)
		result[[z]]$anova2.table <- AOV2
		
		# --- Estimates of GCA, SCA, and REC effects --- #
		
		library(doBy)
		myformula4<- paste(respvar[z], " ~ ", p1," + ",p2, sep = "")
		meandata <- summaryBy(formula(myformula4), data=temp.data)
		
		# serial to parallel of meandata
		mdata <- as.matrix(rep(0,p*p),nrow=p, ncol=p)
		dim(mdata) <- c(p,p)
		
		for (I in 1:p)  {
			for (J in 1:p)   {
				if (I != J) mdata[I,J] <- meandata[(meandata[,match(p1, names(meandata))]==I & meandata[,match(p2, names(meandata))]==J),3]
			} 
		}
		colnames(mdata) <- rownames(mdata) <- levels(temp.data[,match(p1, names(temp.data))])
		
		# printing the matrix of means
		mdata2 <- format(round(mdata, 4),nsmall=4)
		mdata3 <-noquote(format(gsub(" 0.0000", "", mdata2),justify="right"))
		cat("\n\nMATRIX OF MEANS\n")
		print(mdata3)
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
		
		colnames(G_SCA) <- rownames(G_SCA) <- levels(temp.data[,match(p1, names(temp.data))])
		G_SCA2 <-noquote(format(round(G_SCA, digits=4), nsmall=4))
		cat("\n\nGENERAL COMBINING ABILITY EFFECTS, SPECIFIC COMBINING ABILITY EFFECTS\n")
		cat("(above diagonal) AND RECIPROCAL EFFECTS (below diagonal):\n")
		cat("AVERAGED OVER ENVIRONMENTS:\n")
		print(G_SCA2)
		result[[z]]$gcasca.matrix <- G_SCA2
		
		# --- estimates of standard errors --- #
		
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
		
		rownames(VAREST) <- c(' Gi', ' Sij', ' Rij', ' Gi-Gj', ' Sij-Sik', ' Sij-Skl')
		colnames(VAREST) <- c("Std. Error", "LSD")
		cat("\n\nTABLE OF STANDARD ERRORS AND LSDs:\n")
		print(VAREST)
		result[[z]]$stderror.table <- VAREST
		
		#--- Variance Component ---#
		
		Vg <- (GCA.MS-GCAxE.MS-SCA.MS+SCAxE.MS)/(2*r*e*(p-2))
		Vs <- (SCA.MS-SCAxE.MS)/(2*r*e)
		Vr <- (REC.MS-RECxE.MS)/(2*r*e)
		Vge <- (GCAxE.MS-SCAxE.MS)/(2*r*(p-2))
		Vse <- (SCAxE.MS-EMS)/(2*r)
		Vre <- (RECxE.MS-EMS)/(2*r)
		
		VC <- round(rbind(Vg, Vs, Vr, Vge, Vse, Vre, EMS), digits=4)
		colnames(VC) <- c("Estimate")
		rownames(VC) <- c(" GCA", " SCA", " REC"," GCAxE", " SCAxE", " RECxE", " Error")
		TABLE <- as.table(VC)
		cat("\n\nESTIMATES OF VARIANCE COMPONENTS:\n")
		print(TABLE)
		result[[z]]$var.components <-TABLE
		
		if (Vg<0) Vg <- 0
		if (Vs<0) Vs <- 0
		if (Vr<0) Vr <- 0
		if (Vge<0) Vge <- 0
		if (Vse<0) Vse <- 0
		if (Vre<0) Vre <- 0
		
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
		h2N <- VA / VP                 	    # individual based
		Dominance.ratio <- sqrt(2*VD/VA)   
		
		Estimate <- round(rbind(VA, VD, VAE, VDE, h2N, h2B, Dominance.ratio), digits=4)
		colnames(Estimate) <- c("Estimate")
		rownames(Estimate) <- c(" VA", " VD", " VAE", " VDE", " h2-narrow sense", " H2-broad sense", " Dominance Ratio")
		title <- "Estimates of genetic variance components"
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

