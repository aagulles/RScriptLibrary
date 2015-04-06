ssa.pairwise <- function(trmt.levels, model, type = "Tukey", alpha = 0.05, control.level = NULL) {
	library(lme4)
	library(multcomp) #added by RIZAM 082311
	temp <- NULL
	n <- c(1:length(trmt.levels))
	names(n) <- trmt.levels
	if (type == "Dunnett") { contrast1 <- contrMat(n, type = "Dunnett", base = as.numeric(match(control.level, names(n)))) }
	if (type == "Tukey") contrast1 <- contrMat(n, type = "Tukey")
	mc1 <- glht(model, linfct = contrast1)
	interval <- confint(mc1, level = 1 - alpha)
	interval.confint <- as.data.frame(interval$confint)
	signif2 <- subset(interval.confint, as.logical(lwr <= 0 & 0 <= upr) == F)
	if (nrow(signif2) != 0) {
		trmt.comp <- strsplit(rownames(signif2), split = " - ")
		for (i in (1:nrow(signif2))) temp <- rbind(temp, trmt.comp[[i]])
		temp <- data.frame(temp, signif2, row.names = NULL)
		colnames(temp)[1:2] <- c("Trmt[i]", "Trmt[j]")
		return(temp)
#	} else {return(temp <- signif2)}
	} else {
		temp <- as.matrix("(No significant pairwise comparisons.)")
		row.names(temp) <- ""
		return(noquote(temp[1,1]))
	}
	detach("package:multcomp")
	detach("package:lme4")	
}
