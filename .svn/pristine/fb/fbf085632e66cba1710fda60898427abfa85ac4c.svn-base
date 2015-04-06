# -------------------------------------------------------
# SINGLE SITE ANALYSIS
# File Created by: Alaine A. Gulles 07.01.2011
# File Modified by: Alaine A. Gulles 07.01.2011
# Script Created by: Violeta Bartolome
# Script Modified by: Violeta Bartolome
#                     & Alaine A. Gulles 07.01.2011
# --------------------------------------------------------

# --------------------------------------------------------
# ARGUMENTS:
# data - a dataframe
# exptl.design = 1, RCB
#	   = 2, Aug RCB
#	   = 3, Aug LS
#        = 4, Alpha Lattice
#	   = 5, Row-Column 
# respvar - a string; variable name of the response variable
# geno - a string; variable name of the treatment/genotype variable
# row - a string; variable name of the blocking variable or row variable
# column - a string; variable name of the column variable
#        - NULL, if design is RCB, Aug RCB, Alpha Lattice
# rep - a string; variable name of the replication variable
#       - NULL, if design is RCB, Aug RCB, Aug LS, 
# env - a string; variable name of the environment variable
# is.random - logical; indicating wheather genotype/treatment is random or not; default value is FALSE (FIXED factor)
# --------------------------------------------------------


ssa.test <- function(exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol"), 
		data, respvar, geno, row, column = NULL, rep = NULL, env = NULL, is.random = FALSE, exclude = NULL) {
	
	exptl.design <- match.arg(exptl.design)
	design <- c("Randomized Complete Block (RCB)", "Augmented RCB",
			"Augmented Latin Square", "Alpha Lattice", "Row-Column")
	
	if (exptl.design == "RCB" || exptl.design == "AugRCB") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }} 
	if (exptl.design == "AugLS") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "Alpha") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(rep, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "RowCol") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data))) || is.na(match(rep, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	
	if (is.null(env)) {
		env = "EnvLevel"
		data <- cbind(data, EnvLevel=1)
	}
	
	#define all factors #added by RIZAM 08/01/11
	data[,match(geno, names(data))] <- factor(trim.strings(data[,match(geno, names(data))]))
	data[,match(env, names(data))] <- factor(trim.strings(data[,match(env, names(data))]))
	data[,match(row, names(data))] <- factor(trim.strings(data[,match(row, names(data))]))	
	if (exptl.design == "AugLS") { 
#			data[,match(row, names(data))] <- factor(trim.strings(data[,match(row, names(data))]))	
		data[,match(column, names(data))] <- factor(trim.strings(data[,match(column, names(data))]))	
	}
	if (exptl.design == "Alpha") { 
#			data[,match(row, names(data))] <- factor(trim.strings(data[,match(row, names(data))]))	
		data[,match(rep, names(data))] <- factor(trim.strings(data[,match(rep, names(data))]))	
	}
	
	if (exptl.design == "RowCol") { 
#			data[,match(row, names(data))] <- factor(trim.strings(data[,match(row, names(data))]))	
		data[,match(column, names(data))] <- factor(trim.strings(data[,match(column, names(data))]))	
		data[,match(rep, names(data))] <- factor(trim.strings(data[,match(rep, names(data))]))	
	}
	
	
	result <- list()
	for (i in (1:length(respvar))) {
		result[[i]] <- list()
		result[[i]]$respvar <- respvar[i]
		for (j in (1:nlevels(data[,match(env, names(data))]))) {
			temp.data <- data[sort(match(c(respvar[i], geno, row, column, rep, env), names(data)))]
			result[[i]]$site[[j]] <- list()
			result[[i]]$site[[j]]$env <- levels(temp.data[,match(env, names(temp.data))])[j]
			temp.data <- subset(temp.data, temp.data[,match(env, names(temp.data))] == levels(temp.data[,match(env, names(temp.data))])[j])
			result[[i]]$site[[j]]$obsread <- nrow(temp.data)
			temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[i], names(temp.data))]) == F))
			result[[i]]$site[[j]]$obsused <- nrow(temp.data)
			no.reps <- data.frame(n = tapply(eval(parse(text = paste("temp.data$", respvar[i], sep = ""))), eval(parse(text = paste("temp.data$", geno, sep = ""))), FUN = length))
			no.reps <- as.numeric(1/mean(1/no.reps)) 
			result[[i]]$site[[j]]$numreps <- no.reps
			
			# --- CONSTRUCT THE MODEL --- #
			
			# --- CONSTRUCT THE MODEL: RANDOM FACTOR --- #
			if (is.random) {
				if (is.null(exclude)) {
					if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula1 <- paste(respvar[i], " ~ 1 + (1|", row, ") + (1|", geno,")", sep = "") }
					if (exptl.design == "AugLS") { myformula1 <- paste(respvar[i], " ~ 1 + (1|", row, ") + (1|", column, ") + (1|", geno,")", sep = "") }
#					if (exptl.design == "Alpha") { myformula1 <- paste(y, " ~ 1 + (1|", rep,"/", row,") + (1|", geno,")", sep = "") }
#					if (exptl.design == "RowCol") { myformula1 <- paste(y, " ~ 1 + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|", geno,")", sep = "") }
					if (exptl.design == "Alpha") { myformula1 <- paste(respvar[i], " ~ 1 + (1|", rep,"/", row,") + (1|", geno,")", sep = "") } #by RIZAM
					if (exptl.design == "RowCol") { myformula1 <- paste(respvar[i], " ~ 1 + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|", geno,")", sep = "") } #by RIZAM
				} else {
					for (k in (1:nrow(temp.data))) {
						if (is.na(match(temp.data[k,match(geno, names(temp.data))], exclude))) {
							temp.data$Test[k] <- levels(temp.data[,match(geno, names(temp.data))])[temp.data[k,match(geno, names(temp.data))]]
							temp.data$Check[k] <- "0"
						}else {
							temp.data$Test[k] <- "0"
							temp.data$Check[k] <- exclude[match(temp.data[k,match(geno, names(temp.data))], exclude)]
						}
					}
					temp.data$Check <- factor(temp.data$Check)
					temp.data$Test <- factor(temp.data$Test)
					
					if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula1 <- paste(respvar[i], " ~ 1 + Check + (1|", row, ") + (1|Test:Check)", sep = "") }
					if (exptl.design == "AugLS") { myformula1 <- paste(respvar[i], " ~ 1 + Check + (1|", row, ") + (1|", column, ") + (1|Test:Check)", sep = "") }
					if (exptl.design == "Alpha") { myformula1 <- paste(respvar[i], " ~ 1 + Check + (1|", rep,"/", row,") + (1|Test:Check)", sep = "") }
					if (exptl.design == "RowCol") { myformula1 <- paste(respvar[i], " ~ 1 + Check + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|Test:Check)", sep = "") }
				}
			}
			
			# --- CONSTRUCT THE MODEL: FIXED FACTOR --- #
			if (!is.random) {
				if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula1 <- paste(respvar[i], " ~ 1 + ", geno," + (1|", row, ")", sep = "") }
				if (exptl.design == "AugLS") { myformula1 <- paste(respvar[i], " ~ 1 + ", geno," + (1|", row, ") + (1|", column, ")", sep = "") }
				if (exptl.design == "Alpha") { myformula1 <- paste(respvar[i], " ~ 1 + ", geno," + (1|", rep,"/", row,")", sep = "") }
				if (exptl.design == "RowCol") { myformula1 <- paste(respvar[i], " ~ 1 + ", geno," + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "") }
			} ## -- end stmt -- if (!is.random) 
			
			library(lme4) #added by RIZAM 07/29/11
			model <- lmer(formula(myformula1), data = temp.data)
			result[[i]]$site[[j]]$formula <- myformula1
			result[[i]]$site[[j]]$model <- model
			
			# --- VARIANCE COMPONENTS --- #
			varcomp <- NULL
			for (k in (1:length(VarCorr(model)))) { varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[k], Variance = VarCorr(model)[[k]][1], Std.Dev. = attr(VarCorr(model)[[k]], "stddev")[[1]])) }
			varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model), "sc")**2, Std.Dev. = attr(VarCorr(model), "sc")))
			attr(varcomp, "heading") <- "Variance Components for Random Effects\n"
			result[[i]]$site[[j]]$varcomp.table <- varcomp
			
			#for saving variance and num of reps by RIZAM 083111
			result[[i]]$site[[j]]$varcompnRep <- as.data.frame(attr(VarCorr(model), "sc")**2)
			result[[i]]$site[[j]]$varcompnRep$numRep <- result[[i]]$site[[j]]$numreps
			result[[i]]$site[[j]]$varcompnRep$env <- result[[i]]$site[[j]]$env[[1]]
			colnames(result[[i]]$site[[j]]$varcompnRep) <- c(paste(respvar[i],"sigma2",sep="_"),paste(respvar[i],"No.rep",sep="_"),env)
			if (j == 1) {result[[i]]$out.sigma2 <- result[[i]]$site[[j]]$varcompnRep 
			} else {result[[i]]$out.sigma2 <- rbind(result[[i]]$out.sigma2, result[[i]]$site[[j]]$varcompnRep)}
			
			if (is.random) {
				# --- TEST SIGNIFICANCE OF TREATMENT EFFECT USING LRT --- #
				if (is.null(exclude)) { myformula2 <- gsub(paste(" + (1|", geno,")", sep = ""), "", myformula1, fixed = TRUE) 
				} else { myformula2 <- gsub("+ (1|Test:Check)", "", myformula1, fixed = TRUE)  }
				model1 <- lmer(formula(myformula1), data = temp.data, REML = F)
				model2 <- lmer(formula(myformula2), data = temp.data, REML = F)
				anova.table1 <- anova(model2, model1)
				attr(anova.table1, "heading")[3] <- paste("model2: ", myformula2, sep = "")
				attr(anova.table1, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
				attr(anova.table1, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
				
				anova.table2 <- anova.table1[2,c(5,ncol(anova.table1))]
				rownames(anova.table2) <- geno
				attr(anova.table2, "heading")[1] <- attr(anova.table1, "heading")[1]
				attr(anova.table2, "heading")[2] <- paste("Environment Variable: ", env, " = ", levels(temp.data[,match(env, names(temp.data))])[j], sep = "")
				attr(anova.table2, "heading")[3] <- paste("Response Variable: ", respvar[i], "\n", sep = "")
				result[[i]]$site[[j]]$anova.table <- anova.table2
				result[[i]]$site[[j]]$model1 <- model1
				
				# --- ESTIMATE HERITABILITY --- #
#				if (is.null(exclude)) {result[[i]]$site[[j]]$heritability <- varcomp[varcomp[,1] == geno, "Variance"]/(varcomp[varcomp[,1] == geno, "Variance"] + (varcomp[varcomp[,1] == "Residual", "Variance"]/no.reps))
#				} else {result[[i]]$site[[j]]$heritability <- varcomp[varcomp[,1] == "Test:Check", "Variance"]/(varcomp[varcomp[,1] == "Test:Check", "Variance"] + (varcomp[varcomp[,1] == "Residual", "Variance"]/no.reps))}
				if (is.null(exclude)) {herit <- varcomp[varcomp[,1] == geno, "Variance"]/(varcomp[varcomp[,1] == geno, "Variance"] + (varcomp[varcomp[,1] == "Residual", "Variance"]/no.reps))
				} else {herit <- varcomp[varcomp[,1] == "Test:Check", "Variance"]/(varcomp[varcomp[,1] == "Test:Check", "Variance"] + (varcomp[varcomp[,1] == "Residual", "Variance"]/no.reps))}
				#herit <- as.matrix(herit)									#RIZAM 092111
				herit <- as.matrix(round(herit,digits = 2))					#RIZAM 092111
				#herit <- as.table(round(herit,digits = 2))									#RIZAM 122111
				rownames(herit) <- ""                    			#RIZAM 092111
				result[[i]]$site[[j]]$heritability <- herit[1,1]  	#RIZAM 092111
				#result[[i]]$site[[j]]$heritability <- herit[1]  	#RIZAM 122111
				
				# --- PREDICTED MEANS OF GENOTYPES -- #
				if (is.null(exclude)) { sumStat.table <- eval(parse(text = paste("coef(model)$", geno, sep = ""))); sumStat.table <- cbind(rownames(sumStat.table), sumStat.table) 
				} else {
					sumStat.table <- coef(model)$"Test:Check"
					temp.names <- t(as.data.frame(strsplit(rownames(sumStat.table), ":")))
					sumStat.table$Test <- temp.names[,1]
					sumStat.table$Check <- temp.names[,2]
					sumStat.table <- subset(sumStat.table, Test != "0", select = c("Test", "(Intercept)"))
				}
				rownames(sumStat.table) <- NULL
				colnames(sumStat.table) <- c(geno, "Means")
				attr(sumStat.table, "heading") <- paste("Predicted Means of ", geno, "\n", sep = "")
				result[[i]]$site[[j]]$summary.statistic <- sumStat.table
				
				#For saving to file
				result[[i]]$site[[j]]$sum.out <- sumStat.table 												#added by RIZAM 083111
				#			colnames(result[[i]]$site[[j]]$sum.out) <- c(geno,paste(result[[i]]$respvar,"PredMean",sep="_"))
				result[[i]]$site[[j]]$sum.out$Env <- result[[i]]$site[[j]]$env[[1]] 						#added by RIZAM 083111
				colnames(result[[i]]$site[[j]]$sum.out) <- c(geno,paste(result[[i]]$respvar,"PredMean",sep="_"),env)
				
				if (j==1) {result[[i]]$means.out <- result[[i]]$site[[j]]$sum.out							#added by RIZAM 083111
				} else {result[[i]]$means.out <- rbind(result[[i]]$means.out,result[[i]]$site[[j]]$sum.out)}	#added by RIZAM 083111
				
#				# Prediction intervals of Genotypes means #added 090411
#				if (length(levels(temp.data$Geno)) <= 50) {
#				print(dotplot(ranef(model1, postVar=TRUE))$Geno) } else 
#					{ (qqmath(ranef(model1, postVar=TRUE))$Geno) };
				
			}
			
			if (!is.random) {
				# --- DISPLAY ANOVA TABLE --- #
				anova.table <- anova(model)
				attr(anova.table, "heading")[1] <- paste("Analysis of Variance for ", env, " = ", levels(temp.data[,match(env, names(temp.data))])[j], sep = "")
				attr(anova.table, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")
				
				# --- TEST SIGNIFICANCE OF TREATMENT EFFECT USING MAXIMUM LIKELIHOOD RATIO TEST --- #
				myformula2 <- gsub(paste(" + ", geno, sep = ""), "", myformula1, fixed = TRUE)		
				model1 <- lmer(formula(myformula1), data = temp.data, REML = F)
				model2 <- lmer(formula(myformula2), data = temp.data, REML = F)
				anova.table1 <- anova(model2, model1)
				attr(anova.table1, "heading")[3] <- paste("model2: ", myformula2, sep = "")
				attr(anova.table1, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
				attr(anova.table1, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
				
				anova.table2 <- cbind(anova.table, anova.table1[2,5:ncol(anova.table1)])
				attr(anova.table2, "class") <- c("anova", "data.frame")
				attr(anova.table2, "heading")[1] <- attr(anova.table1, "heading")[1]
				attr(anova.table2, "heading")[2] <- paste("Environment Variable: ", env, " = ", levels(temp.data[,match(env, names(temp.data))])[j], sep = "")
				attr(anova.table2, "heading")[3] <- attr(anova.table, "heading")[2]
				result[[i]]$site[[j]]$anova.table <- anova.table2
				
				# --- COMPUTE TREATMENT MEANS --- #
				myformula3 <- gsub("~ 1", "~ 0", myformula1)
				model3 <- lmer(formula(myformula3), data = temp.data)
				sumStat.table <- data.frame(summary(model3)@coefs)[,1:2]
				rownames(sumStat.table) <- gsub(geno,"",rownames(sumStat.table))
				sumStat.table <- cbind(rownames(sumStat.table), sumStat.table)
				rownames(sumStat.table) <- NULL
				colnames(sumStat.table) <- c(geno, "Mean", "StdErrMean")
				result[[i]]$site[[j]]$summary.statistic <- sumStat.table
				result[[i]]$site[[j]]$summary.statistic <- sumStat.table
				
				#For saving to file
				result[[i]]$site[[j]]$sum.out <- sumStat.table 												#added by RIZAM 08/30/11
#				colnames(result[[i]]$site[[j]]$sum.out) <- c(geno, paste(result[[i]]$respvar,"Mean",sep="_"), paste(result[[i]]$respvar,"StdErrMean",sep="_"))
				result[[i]]$site[[j]]$sum.out$Env <- result[[i]]$site[[j]]$env[[1]] 						#added by RIZAM 08/30/11
				colnames(result[[i]]$site[[j]]$sum.out) <- c(geno, paste(result[[i]]$respvar,"Mean",sep="_"), paste(result[[i]]$respvar,"StdErrMean",sep="_"),env)
#				colnames(result[[i]]$site[[j]]$sum.out) <- c(names(result[[i]]$site[[j]]$summary.statistic[1]), paste(result[[i]]$respvar,"Mean",sep="_"), paste(result[[i]]$respvar,"StdErrMean",sep="_"))

				if (j==1) {result[[i]]$meansse.out <- result[[i]]$site[[j]]$sum.out							#added by RIZAM 08/30/11 #rev 031612
				} else {result[[i]]$meansse.out <- rbind(result[[i]]$meansse.out,result[[i]]$site[[j]]$sum.out)}	#added by RIZAM 08/30/11 #rev 031612

				if (j==1) {result[[i]]$means.out <- result[[i]]$site[[j]]$sum.out[-3]							#added by RIZAM 08/30/11 #rev 031612
				} else {result[[i]]$means.out <- rbind(result[[i]]$means.out,result[[i]]$site[[j]]$sum.out[-3])}	#added by RIZAM 08/30/11 #rev 031612
				
			} ### end stmt -- if (!is.random)
			
			result[[i]]$site[[j]]$residuals <- resid(model1)
			result[[i]]$site[[j]]$fitted.values <- fitted(model1)
			result[[i]]$site[[j]]$data <- temp.data
			
			if (is.random) {names.resid <- paste(respvar[i],"resid_random",sep="_")
			} else {names.resid <- paste(respvar[i],"resid_fixed",sep="_")}
			
			#to create data frame for residuals
			if (j==1) {
				result[[i]]$residuals <- as.data.frame(result[[i]]$site[[j]]$residuals)
				names(result[[i]]$residuals) <- names.resid
				result[[i]]$residuals <- cbind(temp.data,result[[i]]$residuals)
			} else {
				resid2 <- as.data.frame(result[[i]]$site[[j]]$residuals)
				names(resid2) <- names.resid
				resid2 <- cbind(temp.data,resid2)
				result[[i]]$residuals <- rbind(result[[i]]$residuals,resid2)			
			} 
		} ## -- end stmt -- for (j in (1:nlevels(data[,match(env, names(data))])))
		
		#to consolidate means and variances for fixed and random?
		
#		if (!is.random) { 
		if (i==1) {means.out.all <- result[[i]]$means.out									#added by RIZAM 083011
		} else {means.out.all <- merge(means.out.all,result[[i]]$means.out, by=c(env,geno))}	#added by RIZAM 083011

		if (!is.random) { #060112
			if (i==1) {meansse.out.all <- result[[i]]$meansse.out										#added by RIZAM 083011 #rev 031612
			} else {meansse.out.all <- merge(meansse.out.all,result[[i]]$meansse.out, by=c(env,geno))}	#added by RIZAM 083011 #rev031612
		} else {
			meansse.out.all <- NULL
		}
		
		if (i==1) {varrep.out.all <- result[[i]]$out.sigma2									#added by RIZAM 083111
		} else {varrep.out.all <- merge(varrep.out.all,result[[i]]$out.sigma2,by=c(env))}		#added by RIZAM 083111
		
		if (exptl.design == "RCB" | exptl.design == "AugRCB") { byVariables <- c(env,geno,row)}
		if (exptl.design == "AugLS") { byVariables <- c(env,geno,row,col)}
		if (exptl.design == "Alpha") { byVariables <- c(env,geno,row,rep)}
		if (exptl.design == "RowCol") { byVariables <- c(env,geno,row,col,rep)}
		
		if (i==1) {resid.out.all <- result[[i]]$residuals									#added by RIZAM 083111
		} else {
			resid.out.all <- merge(resid.out.all,result[[i]]$residuals,by= byVariables,sort = TRUE)
		}
	} ## -- end stmt -- for (i in (1:length(respvar)))
	
	#create resid.out.all for second data set -random
	residuals.data <- resid.out.all
	for (i in (1:length(respvar))) {
		residuals.data[,match(respvar[i], names(residuals.data))] <- NULL
	}
	detach("package:lme4")
	return(list(design = design[match(exptl.design, c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol"))], 
					output = result,
					means = means.out.all,
					meansse = meansse.out.all,
					varrep = varrep.out.all,
					residuals = resid.out.all,
					residuals.data = residuals.data,
					byVars = byVariables))
}
