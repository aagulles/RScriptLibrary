# --------------------------------------------------------------------------------
# ConstructVarCorrTable
# Created by: Alaine Gulles for Internation Rice Research Institute 09.10.2014
# Description: Construction the Variance Table
# --------------------------------------------------------------------------------

ConstructVarCorrTable <- function(model) UseMethod("ConstructVarCorrTable")

ConstructVarCorrTable.default <- function(model) {
     if (class(model) != "lmerMod") { stop("The argument 'model' should be of class lmerMod.") }
     varcomp <- NULL
     for (k in (1:length(VarCorr(model)))) { varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[k], Variance = VarCorr(model)[[k]][1], Std.Dev. = attr(VarCorr(model)[[k]], "stddev")[[1]])) }
     varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model), "sc")**2, Std.Dev. = attr(VarCorr(model), "sc")))
     names(varcomp)[3] <- "Std. Deviation"
     attr(varcomp, "heading") <- "Variance Components for Random Effects\n"
     return(varcomp)
}