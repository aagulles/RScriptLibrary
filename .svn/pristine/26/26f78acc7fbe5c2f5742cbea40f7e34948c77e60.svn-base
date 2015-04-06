# -------------------------------------------------------------------------------------
# designStrip: Generate randomization for strip plot family design.
# Created by: Alaine A. Gulles 09.21.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.13.2012
# -------------------------------------------------------------------------------------

designStrip <- function(vertical, horizontal, sub = NULL, ssub = NULL, r = 2, trial = 1, display = TRUE, file = NULL) UseMethod("designStrip")

designStrip.default <- function(vertical, horizontal, sub = NULL, ssub = NULL, r = 2, trial = 1, display = TRUE, file = NULL) {

	if (missing(vertical)) { stop("The argument 'vertical' is missing.") }
	if (missing(horizontal)) { stop("The argument 'horizontal' is missing.") }

	factorList <- c(FactorList(vertical),FactorList(horizontal))
	randomize <- NULL
	numItems <- prod(length(factorList[[1]]), length(factorList[[2]]))
	for (i in (1:trial)) {
		for (j in (1:r)) {
			vertical <- sample(factorList[[1]], length(factorList[[1]]))
			horizontal <- sample(factorList[[2]], length(factorList[[2]]))
			randomize <- rbind(randomize, data.frame(row.names = NULL, Trial = i, Block = j, rep(vertical, each = length(factorList[[2]])), rep(horizontal, length(factorList[[1]]))))
		}
	}

	if (!is.null(sub)) {
		factorList <- c(factorList, FactorList(sub))
		numItems <- prod(numItems, length(factorList[[3]]))
		book <- randomize
		randomize <- NULL
		for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(factorList[[3]], length(factorList[[3]])))) }
	}
	if (!is.null(ssub)) {
		factorList <- c(factorList, FactorList(ssub))
		numItems <- prod(numItems, length(factorList[[4]]))
		book <- randomize
		randomize <- NULL
		for (i in (1:nrow(book))) { randomize <- rbind(randomize, data.frame(row.names = NULL, book[i,], sample(factorList[[4]], length(factorList[[4]])))) }
	}
	randomize <- data.frame(randomize, PlotNum = 1:numItems)
	colnames(randomize) <- c("Trial", "Block", names(factorList), "PlotNum")

 	if (display) {
		stripLabel <- c("Strip Plot Design", "Strip-Split Plot Design", "Strip-Split-Split Plot Design")
		factorLabel <- c("Vertical", "Horizontal", "Subplot", "Sub-subplot")
		cat(toupper("Design Properties:"),"\n",sep = "")
		cat("\t",stripLabel[length(factorList)-1],"\n\n",sep = "") 
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
		cat("\t","Number of Blocks = ", r, "\n",sep = "")
		for (i in (1:length(factorList))) {
			cat("\t",factorLabel[i]," Factor = ", names(factorList)[i], "\n",sep = "")
			cat("\t","Levels = ", sep = "")
			if (length(factorList[[i]]) <= 5) { cat(paste(factorList[[i]], collapse = ", ", sep = ""), sep = "")
			} else {
				cat(paste(factorList[[i]][1:3], collapse = ", ", sep = ""), sep = "")
				cat(paste(", ...,", factorList[[i]][length(factorList[[i]])]), sep = "")
			}
			cat("\n")
		}
		cat("\n")
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