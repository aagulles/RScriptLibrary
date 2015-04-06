# for STAR
# 4 april 2013
ssa.pairwise <- function(model, type = "Tukey", alpha = 0.05, control.level = NULL) {
	#library(lme4)
	#library(multcomp) #added by RIZAM 082311
	temp <- NULL
     trmt.levels <- levels(model@frame[,2])
	n <- c(1:length(trmt.levels))
	names(n) <- trmt.levels
	if (type == "Dunnett") { contrast1 <- contrMat(n, type = "Dunnett", base = as.numeric(match(control.level, names(n)))) }
	if (type == "Tukey") contrast1 <- contrMat(n, type = "Tukey")
     command <- paste("glht(model, linfct = mcp(", names(model@frame)[2]," = contrast1))", sep = "")
	#mc1 <- glht(model, linfct = mcp(ENTRIES = contrast1))
     mc1 <- eval(parse(text = paste("glht(model, linfct = mcp(", names(model@frame)[2]," = contrast1))", sep = "")))
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
#	detach("package:multcomp")
#	detach("package:lme4")	
}
