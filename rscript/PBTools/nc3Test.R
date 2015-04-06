nc3Test <-
function(design = c("CRD", "RCB"), data, respvar, tester, f2lines, rep=NULL, block=NULL, individual=NULL, environment=NULL) {
	
	options(show.signif.stars=FALSE)
	data <- eval(parse(text = data))
	
	#trim the strings
	respvar <- trim.strings(respvar)
	tester <- trim.strings(tester)
	f2lines <- trim.strings(f2lines)
	block <- trim.strings(block)
	rep <- trim.strings(rep)
	individual <- trim.strings(individual)
	if (!is.null(environment)) {environment <-trim.strings(environment)}
	
	# --- create titles --- #
	cat("\nDESIGN: NORTH CAROLINA EXPERIMENT III IN ",design, "\n", sep="")
	
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
	      crdVars<-c(respvar[i], tester, f2lines, rep, environment)
	      rcbVars<-c(respvar[i], tester, f2lines, block, environment)
	    } else {
	      crdVars<-c(respvar[i], tester, f2lines, rep)
	      rcbVars<-c(respvar[i], tester, f2lines, block)
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
			temp.data[,match(tester, names(temp.data))] <- factor(trim.strings(temp.data[,match(tester, names(temp.data))]))
			temp.data[,match(f2lines, names(temp.data))] <- factor(trim.strings(temp.data[,match(f2lines, names(temp.data))]))
			if (design == "CRD") {temp.data[,match(rep, names(temp.data))] <- factor(trim.strings(temp.data[,match(rep, names(temp.data))])) }
			if (design == "RCB") {temp.data[,match(block, names(temp.data))] <- factor(trim.strings(temp.data[,match(block, names(temp.data))]))	}
			
			# --- check if raw data is balanced. If not, generate estimates for missing values --- #
			temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[i], names(temp.data))]) == FALSE))
			if (design == "CRD") {
				tempDataForAnova<-temp.data[,c(f2lines, tester, rep, respvar[i])]
				balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], f2lines, tester, rep)
			}
			if (design == "RCB") {
				tempDataForAnova<-temp.data[,c(f2lines, tester, block, respvar[i])]
				balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], f2lines, tester, block)
			}
      
	    # --- data summary --- #
	    funcTrialSum <- class.information2(names(temp.data),temp.data)
	    cat("\nDATA SUMMARY: ","\n", sep="")
	    print(funcTrialSum)
	    cat("\nNumber of observations read: ",nrow(temp.data), sep="")
	    cat("\nNumber of missing observations: ",nrow(balancedData)-nrow(temp.data), sep="")
	    
	    # --- ANOVA for NC3 experiment --- #
	    cat("\n\n\nANOVA TABLE FOR THE EXPERIMENT\n")
			if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
				if (nrow(temp.data) == nrow(balancedData)) {
				  anovaRemark <- "REMARK: Raw dataset is balanced."
					dataForAnova<-tempDataForAnova  
				} else {
					if (design == "CRD") {dataForAnova<-estimateMissingData(design="CRD", data=balancedData, respvar[i], f2lines, tester, rep)  }
					if (design == "RCB") {dataForAnova<-estimateMissingData(design="RCB", data=balancedData, respvar[i], f2lines, tester, block)  }
					anovaRemark  <- "REMARK: Raw data and estimates of the missing values are used."
				}  

				if (design == "CRD") {myformula <- paste(respvar[i], " ~ ", tester, "*", f2lines, sep = "")  }
				if (design == "RCB") {myformula <- paste(respvar[i], " ~ ", block, " + ", tester, "*", f2lines, sep = "")  }
				anova.factorial<-summary(aov(formula(myformula), data=dataForAnova))		
				print(anova.factorial)
				cat("-------\n")
				cat(anovaRemark)
				result[[i]]$site[[j]]$nc3.anova <- anova.factorial

			} else {anovaRemark <- "ERROR: Too many missing values. Cannot perform ANOVA for balanced data." 
			        cat("\n",anovaRemark)
			}
		
			# --- LMER for the design --- #
			if (design == "CRD") {myformula1 <- paste(respvar[i], " ~ 1 + ", tester, " + (1|", f2lines, ") + (1|",tester,":",f2lines,")", sep = "") }
			if (design == "RCB") {myformula1 <- paste(respvar[i], " ~ 1 + ", tester, " + (1|", block, ") + (1|", f2lines, ") + (1|",tester,":",f2lines,")", sep = "") }
			library(lme4)
			model <- lmer(formula(myformula1), data = temp.data)
	    result[[i]]$site[[j]]$lmer.result <- summary(model)
		
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
			Vf_m <- as.numeric(varcomp[1,3])
			Vm <- as.numeric(varcomp[2,3])
		
			m <- nlevels(temp.data[,match(f2lines, names(temp.data))])
			if (design == "CRD") {r <- nlevels(temp.data[,match(rep, names(temp.data))])}
			if (design == "RCB") {r <- nlevels(temp.data[,match(block, names(temp.data))]) }
			F<-0
		
			VA <- 4*Vm
			VD <- 2*Vf_m
			if (VD < 0) VD <- 0
			VE <- Ve
			# VE <- Ve -(1/4)*VA -(1/2)*VD      # formula for individual; taken form Kearsey and Pooni
			VP <- VA + VD + VE
			h2N <- VA / VP                      # individual based
			h2B <- (VA + VD) / VP               # individual based
			h2f <- Vm/(Vm+VE/r)                 # family heritability
		
			Dominance.ratio <- sqrt(2*VD/VA)    # will be undefined if VD is negative 
			Estimate <- formatC(rbind(VA, VD, h2N, h2B, h2f, Dominance.ratio), format="f")
			with.colheader<-format(rbind("Estimate", Estimate), justify="right")
			colnames(with.colheader) <- c("")
			rownames(with.colheader) <- c("", " VA", " VD", " h2-narrow sense", " H2-broad sense", " h2-family", " Dominance Ratio")
	    TABLE <- as.table(with.colheader)
	    cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:")
	    print(TABLE)
	    cat("\n")
	    result[[i]]$site[[j]]$genvar.components <- TABLE
		} ## end of for loop (j)
	  cat("\n==============================================================\n")
	}## end of loop (i)
	detach("package:lme4")
	return(list(output = result))
}

