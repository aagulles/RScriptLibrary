# -----------------------------------------------------------------------------------
# DataTransformation: Performs logarithm, natural logarithm and square root 
#                     transformation
# Created by: Alaine A. Gulles 06.07.2011 for International Rice Research Institute
# Modified by: Alaine A. Gulles 08.31.2012
# Note: Modified function of CreateNewVariable function of RCropStat
#     : Need the functions save in powerTransform.R
# -----------------------------------------------------------------------------------
# Arguments:
# data = data frame or a character string which indicate the name of the data frame
# var = a character string which indicate the name of the variable to be transformed
# transformation = the type of transformation to be performed
# targetName = a chara cter string which indicate the variable name to be created
# -----------------------------------------------------------------------------------

DataTransformation <- function(data, var, transformation = c("log", "ln", "sqrt", "exp", "power", "standardize"), targetName = "newVar") UseMethod("DataTransformation")

DataTransformation.default <- function(data, var, transformation = c("log", "ln", "sqrt", "exp", "power", "standardize"), targetName = "newVar") {
	if (is.character(data)) { data <- eval(parse(text = data)) }

	transformation <- match.arg(transformation)
	if (min(data[,var]) < 0) { stop("The variable contains negative observation. Transformation cannot be performed on the variable.") }

	if (min(data[,var]) == 0) { tempvar <- data[var] + 0.5
	} else { tempvar <- data[var] }
	
	if (transformation == "log") 	{ tempvar <- log10(tempvar) 	}
	if (transformation == "ln")  	{ tempvar <- log(tempvar)   	}
	if (transformation == "sqrt") { tempvar <- sqrt(tempvar) 	}
	if (transformation == "exp") 	{ tempvar <- exp(tempvar) 	}
	if (transformation == "standardize"){ tempvar <- (tempvar - mean(tempvar))/sd(tempvar) }	
	if (transformation == "power"){ 
		tempResult <- eval(parse(text = paste("powerTransform(",var," ~ 1, tempvar)", sep = "")))
		tempvar <- basicPower(tempvar, tempResult$lambda)
		remove("tempResult")
	}
	tempName <- make.unique(c(names(data), targetName), sep = "")
	data <- data.frame(data, tempvar)
	colnames(data)[ncol(data)] <- tempName[ncol(data)]
	remove(list = c("transformation", "tempvar", "tempName"))
	return(data)
}