
GEOneStage.test <- function(exptl.design = c("RCB", "Alpha", "RowCol"),  data, respvar, geno, row, column = NULL, rep = NULL, env, is.genoRandom = FALSE) {  
	library(lme4) #RIZAM 091911
	
	prev.opt <- options()$warn
	options(warn = -1)
	exptl.design <- match.arg(exptl.design)
	design <- c("Randomized Complete Block (RCB)", "Alpha Lattice", "Row-Column")

	if (exptl.design == "RCB") 	{ if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }} 
	if (exptl.design == "Alpha") 	{ if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(rep, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "RowCol") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data))) || is.na(match(rep, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}

	data[,match(geno, names(data))] <- factor(data[,match(geno, names(data))])
	data[,match(row, names(data))]  <- factor(data[,match(row, names(data))])
	data[,match(env, names(data))]  <- factor(data[,match(env, names(data))])
	if (!is.null(rep)) data[,match(rep, names(data))] <- factor(data[,match(rep, names(data))])
	if (!is.null(column)) data[,match(column, names(data))] <- factor(data[,match(column, names(data))])

	result <- list()
	for (i in (1:length(respvar))) {
		result[[i]] <- list()
		result[[i]]$respvar <- respvar[i]
		temp.data <- subset(data, subset = (is.na(data[,match(respvar[i], names(data))]) == F))
		sumStat.Env <- summaryStat(temp.data[match(respvar[i], names(temp.data))], temp.data[match(env, names(temp.data))], c("min", "mean", "max", "var", "sd"))
		sumStat.Env <- sumStat.Env[,c(2:ncol(sumStat.Env))]

		result[[i]]$obsread <- nrow(data)
		result[[i]]$obsused <- nrow(temp.data)

		# --- CONSTRUCT THE MODEL --- #
		if (is.genoRandom) trt.stmt <- paste("(1|", geno,")", sep = "") else trt.stmt <- paste(geno, sep = "")		
		if (exptl.design == "RCB")    myformula1 <- paste(respvar[i], " ~ 1 + ", trt.stmt, " + (1|", env,") + (1|", row,":", env,") + (1|", geno,":", env,")", sep = "") 
		if (exptl.design == "Alpha")  myformula1 <- paste(respvar[i], " ~ 1 + ", trt.stmt, " + (1|", env,") + (1|", rep,":", env,") + (1|", rep,":", row,":", env,") + (1|", geno,":", env,")", sep = "")
		if (exptl.design == "RowCol") myformula1 <- paste(respvar[i], " ~ 1 + ", trt.stmt, " + (1|", env,") + (1|", rep,":", env,") + (1|", rep,":", row,":", env,") + (1|", rep,":", column,":", env,") + (1|", geno,":", env,")", sep = "")
		model <- lmer(formula(myformula1), data = temp.data)
		result[[i]]$formula <- myformula1
		result[[i]]$model <- model

		# --- VARIANCE COMPONENTS --- #

		varcomp <- NULL
		for (j in (1:length(VarCorr(model)))) { varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[j], Variance = VarCorr(model)[[j]][1], Std.Dev. = attr(VarCorr(model)[[j]], "stddev")[[1]])) }
		varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model), "sc")**2, Std.Dev. = attr(VarCorr(model), "sc")))
		result[[i]]$varcomp.table <- varcomp

# 		# --- TEST OF SIGNIFICANCE OF GENO EFFECT USING LRT && ANOVA TABLE IF GENO IS FIXED --- #
# 		myformula2 <- gsub(paste(" + ", trt.stmt, sep = ""), "", myformula1, fixed = TRUE)
# 		model1 <- lmer(formula(myformula1), data = temp.data, REML = F)
# 		model2 <- lmer(formula(myformula2), data = temp.data, REML = F)
# 		temp.anova <- anova(model2, model1)
# 		attr(temp.anova, "heading")[3] <- paste("model2: ", myformula2, sep = "")
# 		attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
# 		attr(temp.anova, "heading")[1] <- paste("Test for Genotype Effect using Maximum likelihood ratio test", sep = "")
# 		if (is.genoRandom) { anova.table2 <- temp.anova[2,c(5,ncol(temp.anova))]; rownames(anova.table2) <- geno }
# 		else { anova.table1 <- anova(model1); anova.table2 <- cbind(anova.table1, temp.anova[2,5:ncol(temp.anova)]); attr(anova.table2, "class") <- c("anova", "data.frame") }
# 		attr(anova.table2, "heading")[1] <- paste("Test for Genotype Effect using Maximum likelihood ratio test", sep = "")
# 		attr(anova.table2, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")

    #--- SEPARATED ANOVA AND TEST IF GENO IS FIXED ---# RIZAM 091811
    # --- TEST OF SIGNIFICANCE OF GENO EFFECT USING LRT && ANOVA TABLE IF GENO IS FIXED --- #
  	myformula2 <- gsub(paste(" + ", trt.stmt, sep = ""), "", myformula1, fixed = TRUE)
		model1 <- lmer(formula(myformula1), data = temp.data, REML = F)
		model2 <- lmer(formula(myformula2), data = temp.data, REML = F)
		temp.anova <- anova(model2, model1)
		attr(temp.anova, "heading")[3] <- paste("model2: ", myformula2, sep = "")
		attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
		attr(temp.anova, "heading")[1] <- paste("Test for Genotype Effect using Maximum likelihood ratio test", sep = "")
#     test of gen
		anova.table2 <- temp.anova[2,c(5,ncol(temp.anova))]; rownames(anova.table2) <- geno
  	attr(anova.table2, "heading")[1] <- paste("Test for Genotype Effect using Maximum likelihood ratio test", sep = "")
		attr(anova.table2, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")
    result[[i]]$testsig.Geno <- anova.table2  

    if (!is.genoRandom){
      #ANOVA Table
      anova.table1 <- anova(model1); 
      anova.table2 <- cbind(anova.table1, temp.anova[2,c(5,ncol(temp.anova))]); 
      attr(anova.table2, "class") <- c("anova", "data.frame")
  	  result[[i]]$anova.table1 <- anova.table1
      result[[i]]$anova.table.test <- anova.table2
    }

    # --- TEST OF SIGNIFICANCE OF ENVIRONMENT EFFECT USING LRT --- #
		myformula3 <- gsub(paste(" + (1|", env, ")", sep = ""), "", myformula1, fixed = TRUE)
		model3 <- lmer(formula(myformula3), data = temp.data, REML = F)
		temp.anova <- anova(model3, model1)
		attr(temp.anova, "heading")[3] <- paste("model3: ", myformula3, sep = "")
		attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
		attr(temp.anova, "heading")[1] <- paste("Test for Environment Effect using Maximum likelihood ratio test", sep = "")
		
		anova.table3 <- temp.anova[2, c(5, ncol(temp.anova))]
		rownames(anova.table3) <- env
		attr(anova.table3, "heading")[1] <- paste("Test for Environment Effect by Maximum likelihood ratio test", sep = "")
		attr(anova.table3, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")

		# --- TEST OF SIGNIFICANCE OF GENOTYPE X ENVIRONMENT EFFECT USING LRT --- #
		myformula4 <- gsub(paste(" + (1|", geno, ":", env, ")", sep = ""), "", myformula1, fixed = TRUE)
		model4 <- lmer(formula(myformula4), data = temp.data, REML = F)
		temp.anova <- anova(model4, model1)
		attr(temp.anova, "heading")[3] <- paste("model4: ", myformula4, sep = "")
		attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
		attr(temp.anova, "heading")[1] <- paste("Test for Genotype X Environment Effect using Maximum likelihood ratio test", sep = "")
		                                                  # G X E Effect
		anova.table4 <- temp.anova[2, c(5, ncol(temp.anova))]
		rownames(anova.table4) <- paste(geno, ":", env, sep = "")
		attr(anova.table4, "heading")[1] <- paste("Test for Genotype X Environment Effect using Maximum likelihood ratio test", sep = "")
		attr(anova.table4, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")

		#result[[i]]$testsig.Geno    <- anova.table2
		result[[i]]$testsig.Env     <- anova.table3
		result[[i]]$testsig.GenoEnv <- anova.table4

		# --- PREDICTED MEANS/MEANS OF GENOTYPE --- #
		if (is.genoRandom) { 
			# --- ESTIMATE HERITABILITY --- #

			no.reps <- data.frame(n = tapply(eval(parse(text = paste("temp.data$", respvar[i], sep = ""))), eval(parse(text = paste("temp.data$", geno, sep = ""))), FUN = length))
			no.reps <- as.numeric(1/mean(1/no.reps)) 
			genetic.var <- varcomp[varcomp[,1] == geno, "Variance"]
			ge.var <- varcomp[varcomp[,1] == paste(geno, ":", env, sep = ""), "Variance"]
			resid.var <- varcomp[varcomp[,1] == "Residual", "Variance"]
			heritability <- genetic.var/(genetic.var + (ge.var/nlevels(temp.data[,env])) + (resid.var/(no.reps*nlevels(temp.data[,env]))))
			#heritability <- as.matrix(heritability)         		#RIZAM 091811
			heritability <- as.matrix(round(heritability,digits = 2))       #RIZAM 011112
			rownames(heritability) <- ""                    		#RIZAM 091811
			result[[i]]$heritability <- heritability[1,1]   		#RIZAM 091811
			#result[[i]]$heritability <- heritability

      # --- PREDICTED MEANS --- #
			sumStat.Geno <- eval(parse(text = paste("coef(model)$", geno, sep = ""))); 
			sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno) 
			colnames(sumStat.Geno) <- c(geno, "Mean")
		} 
		else {
			myformula5 <- gsub("~ 1", "~ 0", myformula1, fixed = TRUE)
			model.noint <- lmer(formula(myformula5), data = temp.data)
			sumStat.Geno <- data.frame(summary(model.noint)@coefs)[,1:2]
			rownames(sumStat.Geno) <- gsub(geno,"", rownames(sumStat.Geno))
			sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno)
			colnames(sumStat.Geno) <- c(geno, "Mean", "StdErrMean")
		}
		rownames(sumStat.Geno) <- NULL

		# --- GENOTYPE X ENVIRONMENT MEANS --- #
		# -- GENOTYPE EFFECT:
		intercept <- fixef(model)[[1]]
		if (is.genoRandom) { geno.effect <- eval(parse(text = paste("ranef(model)$", geno, sep = ""))); geno.effect <- data.frame(rownames(geno.effect), geno.effect) }
		else { geno.effect <- as.data.frame(fixef(model)[-1]); geno.effect <- data.frame(gsub(geno, "", rownames(geno.effect)), geno.effect) }
		colnames(geno.effect) <- c(geno, "geno_effect")
		rownames(geno.effect) <- NULL
		# -- ENVIRONMENT EFFECT:
		env.effect <- eval(parse(text = paste("ranef(model)$", env, sep = "")))
		env.effect <- data.frame(gsub(env, "", rownames(env.effect)), env.effect)
		colnames(env.effect) <- c(env, "env_effect")
		rownames(env.effect) <- NULL
		# -- G X E EFFECT:
		GXE.effect <- as.data.frame(eval(parse(text = paste("ranef(model)$'", geno, ":", env, "'",sep = ""))))
		names <- t(as.data.frame(strsplit(rownames(GXE.effect), ":")))
		GXE.effect <- data.frame(names[,1], names[,2], GXE.effect+intercept)
		colnames(GXE.effect) <- c(geno, env, "ge_effect")
		rownames(GXE.effect) <- NULL
		# -- G X E MEANS:
		sumStat.GenoEnv <- merge(GXE.effect, env.effect, by = env, all = TRUE)
		sumStat.GenoEnv <- merge(sumStat.GenoEnv, geno.effect, by = geno, all = TRUE)
		sumStat.GenoEnv <- data.frame(sumStat.GenoEnv[,match(geno, names(sumStat.GenoEnv))], sumStat.GenoEnv[,match(env,names(sumStat.GenoEnv))], rowSums(subset(sumStat.GenoEnv, select = c(ge_effect, env_effect, geno_effect)), na.rm = TRUE))
		colnames(sumStat.GenoEnv) <- c(geno, env, paste(respvar[i], "means", sep = "_"))
		result[[i]]$means.Geno    <- sumStat.Geno
		result[[i]]$means.Env     <- sumStat.Env
		result[[i]]$means.GenoEnv <- sumStat.GenoEnv
		result[[i]]$residuals <- resid(model1)
		result[[i]]$fitted.values <- fitted(model1)
		result[[i]]$data <- temp.data
	} ### end stmt -- for (i in (1:length(respvar))) 
	options(warn = prev.opt)
	detach("package:lme4")
	return(list(design = design[match(exptl.design, c("RCB", "Alpha", "RowCol"))], output = result))
}

