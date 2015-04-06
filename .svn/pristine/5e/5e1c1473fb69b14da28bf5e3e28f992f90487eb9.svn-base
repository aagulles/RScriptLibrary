GETwoStage.test <- function(data, respvar, stderr = NULL, sigma2, numrep, geno, env, weight = c("none", "stderr", "stdmean"), is.genoRandom = FALSE) {
	                                     
 library(lme4) 
	# --- CHECK INPUT --- #
	if (is.na(match(respvar, names(data))) ||  
# 	    is.na(match(stderr , names(data))) || 
	    is.na(match(sigma2 , names(data))) || 
	    is.na(match(numrep, names(data))) || 
	    is.na(match(geno, names(data))) || 
	    is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }
	
	weight <- match.arg(weight)

	data[,match(geno, names(data))] <- factor(data[,match(geno, names(data))])
	data[,match(env, names(data))] <- factor(data[,match(env, names(data))])
	result <- list()

	for (i in (1:length(respvar))) {
		temp.data <- subset(data, subset = (is.na(data[,match(respvar[i], names(data))]) == F))
		if (weight == "stderr")  temp.data$weight <- 1/(temp.data[, match(stderr[i], names(temp.data))]^2) else temp.data$weight <- NULL
		if (weight == "stdmean") temp.data$Y <- temp.data[,match(respvar[i], names(temp.data))]/temp.data[,match(stderr[i], names(temp.data))] else temp.data$Y <- temp.data[,match(respvar[i], names(temp.data))]
		result[[i]] <- list()
		result[[i]]$respvar <- respvar[i]
		result[[i]]$obsread <- nrow(data)
		result[[i]]$obsused <- nrow(temp.data)
		
		# --- CONSTRUCT THE MODEL --- #
		if (is.genoRandom) trt.stmt <- paste("(1|", geno,")", sep = "") else trt.stmt <- paste(geno, sep = "")		
		myformula1 <- paste("Y ~ 1 + ", trt.stmt, " + (1|", env,")", sep = "") 
		model <- lmer(formula(myformula1), weights = temp.data$weight, data = temp.data)
		result[[i]]$formula <- myformula1
		result[[i]]$model <- model

		# --- SUMMARY STATISTICS PER SITE --- #
		sumStat.Env <- data.frame(tapply(temp.data[, match(respvar[i], names(temp.data))], temp.data[, match(env, names(temp.data))], mean))
		sumStat.Env <- data.frame(rownames(sumStat.Env), sumStat.Env)
		colnames(sumStat.Env) <- c(env, "means")
		rownames(sumStat.Env) <- NULL

		numRep <- 1/mean(1/tapply(temp.data[, match(numrep[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean))
		numEnv <- nlevels(temp.data[, match(env, names(temp.data))])
		error.variance <- sum(tapply(temp.data[, match(sigma2[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean) * tapply(temp.data[, match(numrep[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean))/sum(tapply(temp.data[, match(numrep[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean))
		if (weight == "none")    sigma2.ge <- attr(VarCorr(model), "sc")**2 - error.variance
		if (weight == "stderr")  sigma2.ge <- 1
		if (weight == "stdmean") sigma2.ge <- attr(VarCorr(model), "sc")**2 * error.variance

		# --- VARIANCE COMPONENTS --- #
		varcomp <- NULL
		for (j in (1:length(VarCorr(model)))) { varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[j], Variance = VarCorr(model)[[j]][1], Std.Dev. = attr(VarCorr(model)[[j]], "stddev")[[1]])) }
		varcomp <- rbind(varcomp, data.frame(Groups = paste(geno, ":", env, sep = ""), Variance = sigma2.ge, Std.Dev. = sqrt(sigma2.ge)))
		varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model), "sc")**2, Std.Dev. = attr(VarCorr(model), "sc")))
		result[[i]]$varcomp.table <- varcomp

    #--- COMBINED ANOVA AND TEST ---# RIZAM 091811
# 		# --- TEST OF SIGNIFICANCE OF GENO EFFECT USING LRT && ANOVA TABLE IF GENO IS FIXED --- #
# 		myformula2 <- gsub(paste(" + ", trt.stmt, sep = ""), "", myformula1, fixed = TRUE)
# 		model1 <- lmer(formula(myformula1), weights = temp.data$weight, data = temp.data, REML = F)
# 		model2 <- lmer(formula(myformula2), weights = temp.data$weight, data = temp.data, REML = F)
# 		temp.anova <- anova(model2, model1)
# 		attr(temp.anova, "heading")[3] <- paste("model2: ", myformula2, sep = "")
# 		attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
# 		attr(temp.anova, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
# 		if (is.genoRandom) { anova.table2 <- temp.anova[2,c(5,ncol(temp.anova))]; rownames(anova.table2) <- geno } else { anova.table1 <- anova(model1); anova.table2 <- cbind(anova.table1, temp.anova[2,c(5,ncol(temp.anova))]); attr(anova.table2, "class") <- c("anova", "data.frame") }
# 		attr(anova.table2, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
# 		attr(anova.table2, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")
# 		result[[i]]$testsig.Geno <- anova.table2	

    #--- SEPARATED ANOVA AND TEST IF GENO IS FIXED ---# RIZAM 091811
    # --- TEST OF SIGNIFICANCE OF GENO EFFECT USING LRT && ANOVA TABLE IF GENO IS FIXED --- #
		myformula2 <- gsub(paste(" + ", trt.stmt, sep = ""), "", myformula1, fixed = TRUE)
		model1 <- lmer(formula(myformula1), weights = temp.data$weight, data = temp.data, REML = F)
		model2 <- lmer(formula(myformula2), weights = temp.data$weight, data = temp.data, REML = F)
		temp.anova <- anova(model2, model1)
		attr(temp.anova, "heading")[3] <- paste("model2: ", myformula2, sep = "")
		attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
		attr(temp.anova, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
#     test of gen
		anova.table2 <- temp.anova[2,c(5,ncol(temp.anova))]; rownames(anova.table2) <- geno
  	attr(anova.table2, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
		attr(anova.table2, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")
    result[[i]]$testsig.Geno <- anova.table2  

    if (!is.genoRandom){
      #ANOVA Table
      anova.table1 <- anova(model1)
      anova.table2 <- cbind(anova.table1, temp.anova[2,c(5,ncol(temp.anova))]); attr(anova.table2, "class") <- c("anova", "data.frame")
  	  result[[i]]$anova.table1 <- anova.table1
      result[[i]]$anova.table.test <- anova.table2
    }
    #######
    
    # --- TEST OF SIGNIFICANCE OF ENVIRONMENT EFFECT USING LRT --- #
		if (is.genoRandom) {
			myformula3 <- gsub(paste(" + (1|", env, ")", sep = ""), "", myformula1, fixed = TRUE)
			model3 <- lmer(formula(myformula3), weights = temp.data$weight, data = temp.data, REML = F)
			temp.anova <- anova(model3, model1)
			attr(temp.anova, "heading")[3] <- paste("model3: ", myformula3, sep = "")
			attr(temp.anova, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
			attr(temp.anova, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
			anova.table3 <- temp.anova[2, c(5, ncol(temp.anova))]
			rownames(anova.table3) <- env
			attr(anova.table3, "heading")[1] <- paste("Linear Mixed Model fit by Maximum likelihood ratio test", sep = "")
			attr(anova.table3, "heading")[2] <- paste("Response Variable: ", respvar[i], "\n", sep = "")
			result[[i]]$testsig.Env <- anova.table3

			# --- ESTIMATE HERITABILITY --- #
			genetic.var <- varcomp[varcomp[,1] == geno, "Variance"]
			heritability <- genetic.var/(genetic.var + (sigma2.ge/numEnv) + (error.variance/(numRep*numEnv)))
	      		heritability <- as.matrix(round(heritability,digits = 2))         #RIZAM 091811
      			rownames(heritability) <- ""                    #RIZAM 091811
			#result[[i]]$heritability <- heritability
	      		result[[i]]$heritability <- heritability[1,1]   #RIZAM 091811

			# --- PREDICTED MEANS/MEANS OF GENOTYPE --- #
			sumStat.Geno <- eval(parse(text = paste("coef(model)$", geno, sep = ""))); 
			sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno) 
			colnames(sumStat.Geno) <- c(geno, "Mean")
		} 
		else {
			myformula4 <- gsub("~ 1", "~ 0", myformula1, fixed = TRUE)
			model.noint <- lmer(formula(myformula4), weights = temp.data$weight, data = temp.data)
			sumStat.Geno <- data.frame(summary(model.noint)@coefs)[,1:2]
			rownames(sumStat.Geno) <- gsub(geno,"", rownames(sumStat.Geno))
			sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno)
			colnames(sumStat.Geno) <- c(geno, "Mean", "StdErrMean")
		}
		rownames(sumStat.Geno) <- NULL
		result[[i]]$means.Geno <- sumStat.Geno
		result[[i]]$means.Env  <- sumStat.Env
		result[[i]]$residuals  <- resid(model1)
		result[[i]]$fitted.values  <- fitted(model1)
		result[[i]]$data  <- subset(temp.data, select = c(env, geno, respvar[i], stderr[i], sigma2[i], numrep[i]))
		
	}
	detach("package:lme4")
	return(result)
}