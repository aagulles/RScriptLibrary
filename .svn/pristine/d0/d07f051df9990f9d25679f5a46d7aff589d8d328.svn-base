# -------------------------------------------------------------------------------------
# designLSD: Generate randomization Latin Square design (LSD) for single factor or
#            factorial experiment.
# Created by: Alaine A. Gulles 04.11.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.11.2012
# -------------------------------------------------------------------------------------

designLSD <- function(generate, trial = 1, display = TRUE, file = NULL) UseMethod("designLSD")

designLSD.default <- function(generate, trial = 1, display = TRUE, file = NULL) {

	if (is.null(trial) || trial < 1 || is.character(trial) || length(trial) > 1) { stop("The argument 'trial' should be a single numeric value greater than or equal to 1.") }
	if (missing(generate)) { stop("The argument 'generate' is missing.") }
	if (!is.list(generate)) { stop("The argument '' must be a list.") } 	

	tempComb <- GenerateFactor(generate, times = 1)
	randomize <- NULL
	for (i in (1:trial)) {
		tempRandom <- sample(1:nrow(tempComb), nrow(tempComb))
		temp <- rep(tempRandom, each = nrow(tempComb))
		dim(temp) <- c(nrow(tempComb), nrow(tempComb))
		for (j in (2:nrow(temp))) { temp[j,] <- c(temp[1,j:nrow(temp)], temp[1,1:(j-1)]) }
		temp <- temp[sample(1:nrow(tempComb), nrow(tempComb)),]
		temp <- temp[, sample(1:nrow(tempComb), nrow(tempComb))]
		tempSample <- NULL
		for (j in (1:nrow(temp))) { tempSample <- rbind(tempSample, data.frame(Trial = i, Row = 1:nrow(tempComb), Column = j, temp[,j])) }
		newtrmt <- data.frame(tempComb, tempRandom)
		tempSample <- merge(tempSample, newtrmt, by.x = names(tempSample)[ncol(tempSample)], by.y = names(newtrmt)[ncol(newtrmt)])
		tempSample <- tempSample[,2:ncol(tempSample)]
		tempSample <- tempSample[order(tempSample[,c("Column")]),]		
		randomize <- rbind(randomize, tempSample[order(tempSample[,c("Row")]),])
		randomize$PlotNum <- 1:nrow(randomize)
	}
	rownames(randomize) <- 1:nrow(randomize)
	randomize$Trial <- factor(randomize$Trial)
	randomize$Row <- factor(randomize$Row)
	randomize$Column <- factor(randomize$Column)
	
	if (display) {
		cat(toupper("Design Properties:"),"\n",sep = "")
		if (ncol(tempComb) == 1) { cat("\t","Single Factor","\n",sep = "") } else { cat("\t","Factorial Design","\n",sep = "") }
		cat("\t","Latin Square Design","\n\n",sep = "")
		cat(toupper("Design Parameters:"),"\n",sep = "")
		cat("\t","Number of Trials = ", trial, "\n",sep = "")
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