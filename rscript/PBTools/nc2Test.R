nc2Test <-
function(design = c("CRD", "RCB"), data, respvar, female, male, rep=NULL, block=NULL, inbred=TRUE, individual=NULL, environment=NULL) {
	
	options(show.signif.stars=FALSE)
	data <- eval(parse(text = data))
	
	# --- trim the strings --- #
	respvar <- trim.strings(respvar)
	female <- trim.strings(female)
	male <- trim.strings(male)
	block <- trim.strings(block)
	rep <- trim.strings(rep)
	individual <- trim.strings(individual)
	if (!is.null(environment)) {environment <-trim.strings(environment)}
	
	# --- create titles --- #
	if (inbred) {parentsType<-"INBRED"
	} else {parentsType<-"CROSS"}
	cat("\nDESIGN: NORTH CAROLINA EXPERIMENT II IN ",design, " (", parentsType, ")\n", sep="")

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
	      crdVars<-c(respvar[i], female, male, rep, environment)
	      rcbVars<-c(respvar[i], female, male, block, environment)
	    } else {
	      crdVars<-c(respvar[i], female, male, rep)
	      rcbVars<-c(respvar[i], female, male, block)
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
			temp.data[,match(female, names(temp.data))] <- factor(trim.strings(temp.data[,match(female, names(temp.data))]))
			temp.data[,match(male, names(temp.data))] <- factor(trim.strings(temp.data[,match(male, names(temp.data))]))
			if (design == "CRD") {temp.data[,match(rep, names(temp.data))] <- factor(trim.strings(temp.data[,match(rep, names(temp.data))])) }
			if (design == "RCB") {temp.data[,match(block, names(temp.data))] <- factor(trim.strings(temp.data[,match(block, names(temp.data))]))  }
			
			# --- check if raw data is balanced. If not, generate estimates for missing values --- #
			temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[i], names(temp.data))]) == FALSE))
			if (design == "CRD") {
				tempDataForAnova<-temp.data[,c(male, female, rep, respvar[i])]
				balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], male, female, rep)
			}
			if (design == "RCB") {
				tempDataForAnova<-temp.data[,c(male, female, block, respvar[i])]
				balancedData<-generateBalancedData(design="FACTORIAL", data=tempDataForAnova, respvar[i], male, female, block)
			}

			# --- data summary --- #
			funcTrialSum <- class.information2(names(temp.data),temp.data)
			cat("\nDATA SUMMARY: ","\n", sep="")
			print(funcTrialSum)
			cat("\nNumber of observations read: ",nrow(temp.data), sep="")
			cat("\nNumber of missing observations: ",nrow(balancedData)-nrow(temp.data), sep="")

			# --- ANOVA for NC2 experiment --- #
	    cat("\n\n\nANOVA TABLE FOR THE EXPERIMENT\n")
			if ((nrow(temp.data)/nrow(balancedData)) >= 0.90) {
				if (nrow(temp.data) == nrow(balancedData)) {
				  anovaRemark <- "REMARK: Raw dataset is balanced."
					dataForAnova<-tempDataForAnova  
				} else {
					if (design == "CRD") {dataForAnova<-estimateMissingData(design="CRD", data=balancedData, respvar[i], male, female, rep)  }
					if (design == "RCB") {dataForAnova<-estimateMissingData(design="RCB", data=balancedData, respvar[i], male, female, block)  }
					anovaRemark  <- "REMARK: Raw data and estimates of the missing values are used."
				}  

				if (design == "CRD") {myformula <- paste(respvar[i], " ~ ", male, "*", female, sep = "")  }
				if (design == "RCB") {myformula <- paste(respvar[i], " ~ ", block, " + ", male, "*", female, sep = "")  }
				anova.factorial<-summary(aov(formula(myformula), data=dataForAnova))		
				print(anova.factorial)
				cat("-------\n")
				cat(anovaRemark)
				result[[i]]$site[[j]]$nc2.anova <- anova.factorial

			} else {anovaRemark <- "ERROR: Too many missing values. Cannot perform ANOVA for balanced data." 
			        cat("\n",anovaRemark)
              }
      
			# --- LMER for factorial ---#
			if (design == "CRD") {myformula1 <- paste(respvar[i], " ~ 1 + (1|", male, ") + (1|", female, ") + (1|",male,":",female,")", sep = "") }
			if (design == "RCB") {myformula1 <- paste(respvar[i], " ~ 1 + (1|", block, ") + (1|", male, ") + (1|", female, ") + (1|",male,":",female,")", sep = "") }
			library(lme4)
			model <- lmer(formula(myformula1), data = temp.data)
	    result[[i]]$site[[j]]$lmer.result <- summary(model)
		
	    # --- edit format of lmer output before printing --- #
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
	    print(round(summary(model)@coefs, digits=4))
	    cat("\n Random Effects:")
	    print(rematNew)
		
			# --- Estimates of genetic variance components --- #
			varcomp <- summary(model)@REmat
			Ve <- as.numeric(varcomp[varcomp[,1] == "Residual", "Variance"])
			Vm_f <- as.numeric(varcomp[1,3])
			Vf <- as.numeric(varcomp[2,3])
			Vm <- as.numeric(varcomp[3,3])
		
			f <- nlevels(temp.data[,match(female, names(temp.data))])
			m <- nlevels(temp.data[,match(male, names(temp.data))])
			N <- NROW(temp.data)
			if (design == "CRD") {r <- nlevels(temp.data[,match(rep, names(temp.data))])}
			if (design == "RCB") {r <- nlevels(temp.data[,match(block, names(temp.data))]) }
	    
      if (inbred) F<-1
	    else F<-0
		
			VA <- (2/(1+F))*(Vm+Vf)
			VD <- (4/(1+F)^2)*Vm_f
			if (VD < 0) VD <- 0
			VE <- Ve
			VP <- VA + VD + VE
			h2N <- VA / VP                         # plot-mean based, narrow sense
			h2B <- (VA + VD) / VP     	    	   # plot-mean based, broad sense
			# VE <- Ve - (1/2)*VA - (3/4)*(VD)     # formula for individual
			# formula from Kearsey and Pooni
			Dominance.ratio <- sqrt(2*VD/VA)
      
			Estimate <- formatC(rbind(VA, VD, h2N, h2B, Dominance.ratio), format="f")
			with.colheader<-format(rbind("Estimate", Estimate), justify="right")
			colnames(with.colheader) <- c("")
			rownames(with.colheader) <- c("", " VA", " VD", " Narrow sense heritability(plot-mean based)", " Broad sense heritability(plot-mean based)", " Dominance Ratio")
	    TABLE <- as.table(with.colheader)
	    cat("\n\nESTIMATES OF GENETIC VARIANCE COMPONENTS:")
	    print(TABLE)
	    result[[i]]$site[[j]]$genvar.components <- TABLE
		
			# --- Estimates of heritability values --- #
			# --- Family Selection ---#
			H2fm <- Vm/(Vm + Vm_f/f + Ve/(r*f))
			H2ff <- Vf/(Vf + Vm_f/m + Ve/(r*m))
			H2ffs <- (f*Vm + m*Vf + Vm_f)/(f*Vm + m*Vf + Vm_f + Ve/r)
		
			# --- For individual selection --- #
			h2m <- (4/(1+F))*Vm/(Vm + Vm_f + Ve)
			H2m <- ((4/(1+F))*Vm + (4/(1+F)^2)*Vm_f)/(Vm + Vm_f + Ve)
		
			h2f <- (4/(1+F))*Vf/(Vf + Vm_f + Ve)
			H2f <- ((4/(1+F))*Vf + (4/(1+F)^2)*Vm_f)/(Vf + Vm_f + Ve)
			
			h2fs <- (2/(1+F))*(Vm+Vf)/(Vm + Vf + Vm_f + Ve)
			H2fs <- ((2/(1+F))*(Vm+Vf) + (4/(1+F)^2)*Vm_f)/(Vm + Vf + Vm_f + Ve)
			
			rowMale2<-paste("",male)
			rowFemale2<-paste("",female)
			family <- round(rbind(H2fm, H2ff, H2ffs), digits=2)
			narrowsense <- round(rbind(h2m, h2f, h2fs), digits=2)
			broadsense <- round(rbind(H2m, H2f, H2fs), digits=2)
		
			TABLE2 <- cbind(family, narrowsense, broadsense)
			colnames(TABLE2) <- c("Family Selection", "Narrow Sense", "Broad sense")
			rownames(TABLE2) <- c(rowMale2, rowFemale2, " Full-sib")
	    TABLE2_final <- as.table(TABLE2)
	    result[[i]]$site[[j]]$heritability <- TABLE2_final
	    cat("\n\nESTIMATES OF HERITABILITY:\n")
	    print(TABLE2_final) 
	    cat("\n")
		} ## end of for loop (j)
	cat("\n==============================================================\n")  
	}## end of loop (i)
	detach("package:lme4")
	return(list(output = result))
}

