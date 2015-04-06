designLattice <- function(generate, r = 2, trial = 1, file = NULL) UseMethod("designLattice")

designLattice.default <- function(generate, r = 2, trial = 1, file = NULL) {
	
	tempComb <- GenerateFactor(generate, times = 1)
	randomize <- NULL
	for (i in (1:trial)) {
		temp <- randomizeLattice(sqrt(nrow(tempComb)), r)
		randomize <- rbind(randomize, data.frame(Trial = i, temp$randomization, PlotNum = 1:nrow(tempComb)))
	}
	randomize[,"Trial"] <- factor(randomize[,"Trial"])
	
	cat(toupper("Design Properties:"),"\n",sep = "")
	cat("\t","Incomplete Block Design","\n",sep = "") 
	if (r == (sqrt(nrow(tempComb)) + 1)) { cat("\t","Balanced Lattice Design","\n\n",sep = "") } else { cat("\t","Partially Balanced Lattice Design","\n\n",sep = "") } 
	cat(toupper("Design Parameters:"),"\n",sep = "")
	cat("\t","Number of Trials = ", trial, "\n",sep = "")
	cat("\t","Number of Treatments = ", nrow(tempComb), "\n",sep = "")
	cat("\t","Number of Replicates = ", r, "\n",sep = "")
	cat("\t","Plots per Block (Block Size) = ", sqrt(nrow(tempComb)), "\n",sep = "")
	cat("\t","Block per Replicate = ", sqrt(nrow(tempComb)), "\n",sep = "")
	#cat("Results of Randomization:\n")
	#printDataFrame(randomize)
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

