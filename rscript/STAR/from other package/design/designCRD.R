# -------------------------------------------------------------------------------------
# designCRD: Generate randomization for complete block design.
# Created by: Alaine A. Gulles 04.11.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.22.2012
# -------------------------------------------------------------------------------------

designCRD <- function(generate, r = 2, trial = 1, display = TRUE, file = NULL) UseMethod("designCRD")

designCRD.default <- function(generate, r = 2, trial = 1, display = TRUE, file = NULL) {

	if (is.null(trial) || trial < 1 || is.character(trial) || length(trial) > 1) { stop("The argument 'trial' should be a single numeric value greater than or equal to 1.") }
	if (is.null(r) || r < 2 || is.character(r) || length(r) > 1) { stop("The argument 'r' should be a single numeric value greater than or equal to 2.") }
	if (missing(generate)) { stop("The argument 'generate' is missing.") }
	if (!is.list(generate)) { stop("The argument 'generate' must be a list.") }
	
	tempComb <- GenerateFactor(generate, times = r)
	randomize <- NULL
	for (i in (1:trial)) {
		temp <- data.frame(Trial = as.character(i), tempComb, PlotNum = sample(nrow(tempComb), nrow(tempComb), replace = FALSE))
		randomize <- rbind(randomize, temp[order(temp[,"PlotNum"]),])
	}
	rownames(randomize) <- 1:nrow(randomize)

	## display in the console the output
	if (display) { 
		cat(toupper("Design Properties:"),"\n",sep = "")
		if (ncol(tempComb) == 1) { cat("\t","Single Factor","\n",sep = "") } else { cat("\t","Factorial Design","\n",sep = "") }
		cat("\t","Completely Randomized Design","\n\n",sep = "")
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		cat("\t","Number of Replicates = ", r, "\n",sep = "")
		if (ncol(tempComb) == 1) {
			cat("\t","Treatment Name = ", names(tempComb)[1], "\n",sep = "")
			cat("\t","Treatment Levels = ", sep = "")
			if (nlevels(tempComb[,1]) <= 5) { cat(paste(levels(tempComb[,1]), collapse = ", ", sep = ""), sep = "")
			} else {
				cat(paste(levels(tempComb[,1])[1:3], collapse = ", ", sep = ""), sep = "")
				cat(paste(", ...,", levels(tempComb[,1])[nlevels(tempComb[,1])]), sep = "")
			}
			cat("\n\n")
		} else {
			for (i in (1:ncol(tempComb))) {
				cat("\t","Factor ",i," = ", names(tempComb)[i], "\n",sep = "")
				cat("\t","Levels = ", sep = "")
				if (nlevels(tempComb[,i]) <= 5) { cat(paste(levels(tempComb[,i]), collapse = ", ", sep = ""), sep = "")
				} else {
					cat(paste(levels(tempComb[,i])[1:3], collapse = ", ", sep = ""), sep = "")
					cat(paste(", ...,", levels(tempComb[,i])[nlevels(tempComb[,i])]), sep = "")
				}
				cat("\n")
			}
			cat("\n")
		}
		#cat("Results of Randomization:\n")
		#printDataFrame(randomize)
	}

	if (!is.null(file)) {
		tempFile <- strsplit(file, split = "\\.")[[1]]
		tempExt <- tolower(tempFile[length(tempFile)])
		if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
		newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
		newFile <- paste(newFile, tempExt, sep = ".")
		if (tempExt == "csv") { write.csv(randomize, file = newFile, row.names = FALSE)
		} else { save(randomize, file = newFile) }
	} else {
		cat("Results of Randomization:\n")
		printDataFrame(randomize)	
	}
	return(invisible(randomize))
}
