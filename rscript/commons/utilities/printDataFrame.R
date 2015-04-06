# RCropStat Utilities

# -----------------------------------------------------------------------------------
# printDataFrame: Print the dataframe
# Created by: Alaine A. Gulles 04.30.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.10.2012
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# ARGUMENTS:
# dataFrame = the data frame to be printed
# -----------------------------------------------------------------------------------

printDataFrame <- function(dataFrame, border = TRUE) {

	if (!is.data.frame(dataFrame)) { stop("The argument 'dataFrame' should be of class data.frame.") }
	dataChar <- DataAttribute(dataFrame)[,1:2]
	dataChar[,1] <- as.character(dataChar[,1])

	colWidth <- NULL
	onlyDecimal <- NULL
	for (i in (1:nrow(dataChar))) {
		if (dataChar[i,2] == "numeric") {
			if (all(is.na(dataFrame[,i]))) { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 1)
			} else {
				if (all(dataFrame[,i] < 1) && all(dataFrame[,i] > -1)) { 
					colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 5) + 1)
					onlyDecimal <- c(onlyDecimal, i)
				} else {
					if (any(is.na(is.wholenumber(dataFrame[,i])))) {
						tempsubdata <- subset(dataFrame[,i], dataFrame[,i] != Inf)
						if (all(is.wholenumber(tempsubdata))) { 
							dataChar[i,2] <- "integer" 
							colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), nchar(max(tempsubdata))) + 1)
						} else { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(tempsubdata,0))) + 3) + 1) }
					} else {
						if (all(is.wholenumber(dataFrame[,i]))) { 
							dataChar[i,2] <- "integer" 
							colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), nchar(max(dataFrame[,i]))) + 1)
						} else { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 3) + 1) }
					}
				}	
			}
		} else {
			if (dataChar[i,2] == "integer") { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), nchar(max(dataFrame[,i]))) + 1)
			} else { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 1) }
		}
		
	}
	
	if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")	

	for (j in (1:ncol(dataFrame))) {
		if (dataChar[j,2] == "numeric" || dataChar[j,2] == "integer") { 
			cat(formatC(dataChar[j,1], width = colWidth[j], format = "s"), sep = "")
		} else { cat(formatC(dataChar[j,1], width = colWidth[j], format = "s", flag = "-"), sep = "") }
		if (j == ncol(dataFrame)) { 
			cat("\n")
			#if (nrow(dataFrame) >= 2) cat("\n\n") else cat("\n")
		} else { cat(formatC("", width = 1, format = "s")) }
	}

	if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")	

	for (i in (1:nrow(dataFrame))) {
		for (j in (1:ncol(dataFrame))) {
			if (dataChar[j,2] == "numeric") {
				if (is.na(match(j, onlyDecimal))) { 
					if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") }
					else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = 2), sep = "") }
				} else { 
					if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") }
					else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = 4), sep = "") }
				}
			} else {
				if (dataChar[j,2] == "integer") { 
					if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") } 
					else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "d"), sep = "") }
				} else { 
					if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s", flag = "-"), sep = "") }
					else { cat(formatC(as.character(dataFrame[i,j]), width = colWidth[j], format = "s", flag = "-"), sep = "") }	
				}
			}
			if (j == ncol(dataFrame)) { cat("\n") } else { cat(formatC("", width = 1, format = "s")) }
		}
	}
	if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")	

}
