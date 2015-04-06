nc3TestME <-
function(design = c("CRD", "RCB"), data, respvar, tester, f2lines, rep=NULL, block=NULL, individual=NULL, environment=NULL){
	
  options(show.signif.stars=FALSE)
  data <- eval(parse(text = data))
	
	#trim the strings
	respvar <- trim.strings(respvar)
	tester <- trim.strings(tester)
	f2lines <- trim.strings(f2lines)
	block <- trim.strings(block)
	rep <- trim.strings(rep)
	individual <- trim.strings(individual)
	environment <-trim.strings(environment)
  
  # --- create titles --- #
  cat("\nMULTIPLE ENVIRONMENT ANALYSIS\n")
  cat("\nDESIGN: NORTH CAROLINA EXPERIMENT III IN ",design, "\n", sep="")
	
	#define factors
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
		
		# --- LMER for the design --- #
		if (design == "CRD") {myformula1 <- paste(respvar[i], " ~ 1 + (", tester, ") + (1|", environment, ") + (1|", f2lines, ") + (1|", tester, ":", f2lines, ") + (1|", environment, ":", tester, ") + (1|", environment, ":", f2lines, ") + (1|", environment, ":", f2lines, ":", tester,")", sep = "") }
		if (design == "RCB") {myformula1 <- paste(respvar[i], " ~ 1 + (", tester, ") + (1|", environment, ") + (1|", block, ":", environment, ") + (1|", f2lines, ") + (1|", tester, ":", f2lines, ") + (1|", environment, ":", tester, ") + (1|", environment, ":", f2lines, ") + (1|", environment, ":", f2lines, ":", tester,")", sep = "") }
		library(lme4)
		model <- lmer(formula(myformula1), data = temp.data)
		result[[i]]$lmer.result <- summary(model)
		
		# --- edit format of lmer output before printing --- #
		fixeff.table<-round(summary(model)@coefs, digits=4)
		rownames(fixeff.table)<-c(" (Intercept)", paste("",tester))
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
		cat("\n Fixed Effects:\n")
		print(fixeff.table)
		cat("\n Random Effects:")
		print(rematNew)
		
		#--- Estimates of genetic variance components ---#
		varcomp <- summary(model)@REmat
		Ve <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
		Vf_m <- as.numeric(varcomp[3,3])
		Vm <- as.numeric(varcomp[4,3])
		Vef_m <- as.numeric(varcomp[1,3])
		Vem <- as.numeric(varcomp[2,3])
		Vef <- as.numeric(varcomp[5,3])
		
		m <- nlevels(temp.data[,match(f2lines, names(temp.data))])
		e <- nlevels(temp.data[,match(environment, names(temp.data))])	
		if (design == "CRD") {r <- nlevels(temp.data[,match(rep, names(temp.data))])}
		if (design == "RCB") {r <- nlevels(temp.data[,match(block, names(temp.data))]) }
		F<-0
						
		VA <- 4*Vm
		VD <- 2*Vf_m
		if (VD < 0) VD <- 0
		VAE <- 4*Vem
		VDE <- 2*Vef_m
		if (VDE < 0) VDE <- 0
		VE <- Ve
		VP <- VA + VD + VAE + VDE + VE
		h2B <- (VA + VD) / VP            	    # individual based
		h2N <- VA / VP             		        # individual based
		h2f <- Vm /(Vm + VE/(r*e))   		      # family heritability
		Dominance.ratio <- sqrt(2*VD/VA)
		
		Estimate <- formatC(rbind(VA, VAE,  VD, VDE, h2N, h2B, h2f, Dominance.ratio), format="f")
		with.colheader<-format(rbind("Estimate", Estimate), justify="right")
		colnames(with.colheader) <- c("")
		rownames(with.colheader) <- c("", " VA", " VAxE", " VD", " VDxE", " h2-narrow sense", " H2-broad sense", " h2-family", " Dominance Ratio")
		TABLE <- as.table(with.colheader)
		cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:")
		print(TABLE)
		result[[i]]$genvar.components <- TABLE
		cat("\n")
		cat("\n==============================================================\n")
	}## end of loop (i)
	detach("package:lme4")
  	return(list(output = result))
}

