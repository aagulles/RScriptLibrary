designAlphaLattice <- function(generate, blksize, r = 2, trial = 1, file = NULL) UseMethod("designAlphaLattice")

designAlphaLattice.default <- function(generate, blksize, r = 2, trial = 1, file = NULL) {
	
	tempComb <- FactorList(generate)
	numBlk <- length(tempComb[[1]])/blksize
	withAlphaDesign <- FALSE

	#if (length(tempComb[[1]]) %% blksize != 0) { stop("The size of the block is not appropriate.")	}
	if (r < 2 || r > 4) { stop("The valid value of the number of replicate should be between 2 to 4.")}
	if (r == 2) {  
		if (blksize > numBlk) { stop("For r = 2,  block size should be less than or equal to number of blocks (number of treatments/block size).") }
	}
	if (r == 3) {
		if (r%%2 == 0) { if (blksize > numBlk) { stop("For r = 3,  block size should be less than or equal to number of blocks (number of treatments/block size).") }
 		}  else { if (blksize >= numBlk) { stop("For r = 3,  block size should be less than the number of blocks (number of treatments/block size).") }}
	}
	if (r == 4) {
		if (numBlk%%2 != 0) {
			if (numBlk%%3 != 0) {
				if (blksize > numBlk) { stop("For r = 4, number of block should be an odd number but not multiple of 3 and block size should be less than or equal to number of blocks (number of treatments/block size).")	 }
			} else {
				stop("For r = 4, number of block should be an odd number but not multiple of 3 and block size should be less than or equal to number of blocks (number of treatments/block size).")	
			}
		} else { stop("For r = 4, the number of replicate should be an odd number.") }
	}

	randomize <- NULL
	for (i in (1:trial)) {
		capture.output(temp <- design.alpha(trt = tempComb[[1]], k = blksize, r))
		temp$book[,c("replication", "block", "cols", "tempComb[[1]]", "plots")]
		randomize <- rbind(randomize, data.frame(Trial = i, temp$book[,c("replication", "block", "cols", "tempComb[[1]]", "plots")]))
	}
	randomize[,"Trial"] <- factor(randomize[,"Trial"])
	randomize[,"replication"] <- factor(randomize[,"replication"])
	colnames(randomize)[5] <- names(tempComb)
	
	cat(toupper("Design Properties:"),"\n",sep = "")
	cat("\t","Incomplete Block Design","\n",sep = "") 
	cat("\t","Alpha Lattice Design","\n\n",sep = "") 
	cat(toupper("Design Parameters:"),"\n",sep = "")
	cat("\t","Number of Trials = ", trial, "\n",sep = "")
	cat("\t","Number of Treatments = ", length(tempComb[[1]]), "\n",sep = "")
	cat("\t","Number of Replicates = ", r, "\n",sep = "")
	cat("\t","Plots per Block (Block Size) = ", blksize, "\n",sep = "")
	cat("\t","Block per Replicate = ", length(tempComb[[1]])/blksize, "\n\n",sep = "")
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

