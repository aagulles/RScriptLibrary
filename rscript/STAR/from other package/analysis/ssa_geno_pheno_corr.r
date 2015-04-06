genoNpheno.corr <- function(exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol"), 
			    data, respvar, geno, row, column = NULL, rep = NULL, env, exclude = NULL){
	library(lme4) #092711
	exptl.design <- match.arg(exptl.design)
	if (length(respvar) <= 1) stop("Correlation cannot be performed. At least two response variable is required.")
	if (exptl.design == "RCB" || exptl.design == "AugRBC") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }} 
	if (exptl.design == "AugLS") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "Alpha") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(rep, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "RowCol") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data))) || is.na(match(rep, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}

	data[,match(geno, names(data))] <- factor(trim.strings(data[,match(geno, names(data))]))
	data[,match(env, names(data))] <- factor(trim.strings(data[,match(env, names(data))]))
	data[,match(row, names(data))] <- factor(trim.strings(data[,match(row, names(data))]))
	if (!is.null(column)) data[,match(column, names(data))] <- factor(trim.strings(data[,match(column, names(data))]))
	if (!is.null(rep)) data[,match(rep, names(data))] <- factor(trim.strings(data[,match(rep, names(data))]))
	
	design <- c("Randomized Complete Block (RCB)", "Augmented RCB",
			"Augmented Latin Square", "Alpha Lattice", "Row-Column")

	GenoCorr <- list()
	PhenoCorr <- list()
	NumObs <- list()

	if (is.null(exclude)) {
		if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula3 <- paste("Y ~ 1 + (1|", row, ") + (1|", geno,")", sep = "") }
		if (exptl.design == "AugLS") { myformula3 <- paste("Y ~ 1 + (1|", row, ") + (1|", column, ") + (1|", geno,")", sep = "") }
		if (exptl.design == "Alpha") { myformula3 <- paste("Y ~ 1 + (1|", rep,"/", row,") + (1|", geno,")", sep = "") }
		if (exptl.design == "RowCol") { myformula3 <- paste("Y ~ 1 + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|", geno,")", sep = "") }
	}else {
		if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula3 <- paste("Y ~ 1 + Check + (1|", row, ") + (1|Test:Check)", sep = "") }
		if (exptl.design == "AugLS") { myformula3 <- paste("Y ~ 1 + Check + (1|", row, ") + (1|", column, ") + (1|Test:Check)", sep = "") }
		if (exptl.design == "Alpha") { myformula3 <- paste("Y ~ 1 + Check + (1|", rep,"/", row,") + (1|Test:Check)", sep = "") }
		if (exptl.design == "RowCol") { myformula3 <- paste("Y ~ 1 + Check + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|Test:Check)", sep = "") }
	}

	for (i in (1:nlevels(data[,match(env, names(data))]))) {
		GenoCorr[[i]] <- matrix(NA, nrow = length(respvar), ncol = length(respvar))
		PhenoCorr[[i]] <- matrix(NA, nrow = length(respvar), ncol = length(respvar))
		NumObs[[i]] <- matrix(NA, nrow = length(respvar), ncol = length(respvar))
		temp.data <- subset(data, data[,match(env, names(data))] == levels(data[,match(env, names(data))])[i])
		for (j in (1:(length(respvar) - 1))) {
			for (k in (1:length(respvar))) {
				if (k > j) {
					temp.data$Y <- temp.data[,respvar[[j]]] + temp.data[,respvar[[k]]]
					mydata <- subset(temp.data, subset = (is.na(Y) == F)) 
					if (!is.null(exclude)) {
						trmt.label <- "Test:Check"
						for (l in (1:nrow(mydata))) {
							if (is.na(match(mydata[l,match(geno, names(mydata))], exclude))) {
								mydata$Test[l] <- levels(mydata[,match(geno, names(mydata))])[mydata[l,match(geno, names(mydata))]]
								mydata$Check[l] <- "0"
							}
							else {
								mydata$Test[l] <- "0"
								mydata$Check[l] <- exclude[match(mydata[l,match(geno, names(mydata))], exclude)]
							}
						}
						mydata$Test <- factor(mydata$Test)
						mydata$Check <- factor(mydata$Check)
					} else { trmt.label <- geno }
					myformula1 <- gsub("Y", respvar[j], myformula3, fixed = TRUE)
					myformula2 <- gsub("Y", respvar[k], myformula3, fixed = TRUE)
					
					# --- DEFINE MODEL --- #
					model1 <- lmer(formula(myformula1), data = mydata)
					model2 <- lmer(formula(myformula2), data = mydata)
					model3 <- lmer(formula(myformula3), data = mydata)
					# --- COMPUTE HARMONIC MEANS --- #
					no.reps <- data.frame(n = tapply(eval(parse(text = paste("mydata$Y"))), eval(parse(text = paste("mydata$", geno, sep = ""))), FUN = length))
					no.reps <- as.numeric(1/mean(1/no.reps)) 
					# --- COMPUTE GENOTYPIC AND PHENOTYPIC VARIANCES --- #
					varcomp1 <- NULL
					for (l in (1:length(VarCorr(model1)))) { varcomp1 <- rbind(varcomp1, data.frame(Groups = names(VarCorr(model1))[l], Variance = VarCorr(model1)[[l]][1], Std.Dev. = attr(VarCorr(model1)[[l]], "stddev")[[1]])) }
					varcomp1 <- rbind(varcomp1, data.frame(Groups = "Residual", Variance = attr(VarCorr(model1), "sc")**2, Std.Dev. = attr(VarCorr(model1), "sc")))
					genetic.var1 <- varcomp1[varcomp1[,1] == trmt.label, "Variance"]
					pheno.var1 <- varcomp1[varcomp1[,1] == trmt.label, "Variance"] + (varcomp1[varcomp1[,1] == "Residual", "Variance"]/no.reps)

					varcomp2 <- NULL
					for (l in (1:length(VarCorr(model2)))) { varcomp2 <- rbind(varcomp2, data.frame(Groups = names(VarCorr(model2))[l], Variance = VarCorr(model2)[[l]][1], Std.Dev. = attr(VarCorr(model2)[[l]], "stddev")[[1]])) }
					varcomp2 <- rbind(varcomp2, data.frame(Groups = "Residual", Variance = attr(VarCorr(model2), "sc")**2, Std.Dev. = attr(VarCorr(model2), "sc")))
					genetic.var2 <- varcomp2[varcomp2[,1] == trmt.label, "Variance"]
					pheno.var2 <- varcomp2[varcomp2[,1] == trmt.label, "Variance"] + (varcomp2[varcomp2[,1] == "Residual", "Variance"]/no.reps)

					varcomp3 <- NULL
					for (l in (1:length(VarCorr(model3)))) { varcomp3 <- rbind(varcomp3, data.frame(Groups = names(VarCorr(model3))[l], Variance = VarCorr(model3)[[l]][1], Std.Dev. = attr(VarCorr(model3)[[l]], "stddev")[[1]])) }
					varcomp3 <- rbind(varcomp3, data.frame(Groups = "Residual", Variance = attr(VarCorr(model3), "sc")**2, Std.Dev. = attr(VarCorr(model3), "sc")))
					genetic.var <- (varcomp3[varcomp3[,1] == trmt.label, "Variance"] - genetic.var1 - genetic.var2)/2
					pheno.var <- genetic.var + (((varcomp3[varcomp3[,1] == "Residual", "Variance"] - varcomp1[varcomp1[,1] == "Residual", "Variance"] - varcomp2[varcomp2[,1] == "Residual", "Variance"])/2)/no.reps)

					GenoCorr[[i]][j,k] <- genetic.var/sqrt(genetic.var1*genetic.var2)
					PhenoCorr[[i]][j,k] <- pheno.var/sqrt(pheno.var1*pheno.var2)
					NumObs[[i]][j,k] <- length(mydata$Y)
					GenoCorr[[i]][k,j] <- GenoCorr[[i]][[j,k]]
					PhenoCorr[[i]][k,j] <- PhenoCorr[[i]][j,k]
					NumObs[[i]][k,j] <- NumObs[[i]][j,k]
				}  ## end stmt -- if (k > j)
			} ## end stmt -- for (k in (1:length(respvar)))
			GenoCorr[[i]] <- as.table(GenoCorr[[i]])
			colnames(GenoCorr[[i]]) <- rownames(GenoCorr[[i]]) <- respvar
			PhenoCorr[[i]] <- as.table(PhenoCorr[[i]])
			colnames(PhenoCorr[[i]]) <- rownames(PhenoCorr[[i]]) <- respvar
			NumObs[[i]] <- as.table(NumObs[[i]])
			colnames(NumObs[[i]]) <- rownames(NumObs[[i]]) <- respvar
		} ## end stmt -- for (j in (1:(length(respvar) - 1)))
	} ## end stmt -- for (i in (1:nlevels(data[,match(env, names(data))])))
	detach("package:lme4")
	return(list(GenoCorr = GenoCorr, PhenoCorr = PhenoCorr, NumObs = NumObs))
}
