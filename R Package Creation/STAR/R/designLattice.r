# -------------------------------------------------------------------------------------
# designLattice: Generate randomization for augmented design in LSD.
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.22.2013 
# -------------------------------------------------------------------------------------

designLattice <- function(generate, r = 2, trial = 1, numFieldRow = 1, serpentine = FALSE, file = NULL) UseMethod("designLattice")

designLattice.default <- function(generate, r = 2, trial = 1, numFieldRow = 1, serpentine = FALSE, file = NULL) {
	
	tempComb <- GenerateFactor(generate, times = 1)
     #if (nrow(tempComb)%%numFieldRow != 0) {  stop("The number of field row should be divisible by the square of the number of treatments.") }
	randomize <- NULL
     plan <- list()
     plotNum <- NULL
	tempPlotNum <- NULL
	tempBlkNum <- NULL
     numRepRow <- numFieldRow/sqrt(nrow(tempComb))
     numRepCol <- r/numRepRow
     numFieldCol <- (nrow(tempComb)*r)/numFieldRow
	for (i in (1:trial)) {
          plan[[i]] <- matrix(0, nrow = numFieldRow, ncol = numFieldCol)
          if (i == 1) { 
               plotNum <- plan[[i]]
               tempPlotNum <- plan[[i]]
               tempBlkNum <- plan[[i]]
          }
		temp <- randomizeLattice(sqrt(nrow(tempComb)), r)
          for (j in (1:r)) {
               if (j%%numRepCol != 0) { colIndex <- j%%numRepCol } else { colIndex <- numRepCol }
               if (colIndex == 1) { colIndexLL <- colIndex  } else { colIndexLL <- colIndexLL + sqrt(nrow(tempComb)) }
               colIndexUL <- colIndexLL + sqrt(nrow(tempComb)) - 1
               rowIndex <- ceiling(j/numRepCol)
               rowIndexLL <- (rowIndex * sqrt(nrow(tempComb))) - sqrt(nrow(tempComb)) + 1
               rowIndexUL <- rowIndexLL + sqrt(nrow(tempComb)) - 1 
               plan[[i]][rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- temp$plan[[j]]
               if (i == 1) {
                    tPlotNum <- matrix(as.numeric(paste(j, paste(c(rep(0, max(nchar(1:nrow(tempComb))))), collapse = ""), sep = ""))+1:nrow(tempComb), 
                                       nrow = sqrt(nrow(tempComb)), ncol = sqrt(nrow(tempComb)), byrow = TRUE)
                    if (serpentine) { for (k in seq(2, nrow(tPlotNum), by = 2)) tPlotNum[k, ] <- rev(tPlotNum[k, ]) }     
                    plotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- tPlotNum
                    tempPlotNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- matrix(1:nrow(tempComb), nrow = sqrt(nrow(tempComb)), ncol = sqrt(nrow(tempComb)), byrow = TRUE)
                    tempBlkNum[rowIndexLL:rowIndexUL, colIndexLL:colIndexUL] <- matrix(j, nrow = sqrt(nrow(tempComb)), ncol = sqrt(nrow(tempComb)), byrow = TRUE)
               }
          }
          tempBook <- merge(merge(as.data.frame.table(plan[[i]]), as.data.frame.table(plotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2")),
                            as.data.frame.table(tempPlotNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2"))
          tempBook <- suppressWarnings(merge(tempBook, as.data.frame.table(tempBlkNum), by.x = c("Var1", "Var2"), by.y = c("Var1", "Var2")))
          names(tempBook) <- c("FieldRow", "FieldCol", "Treatment", "PlotNum", "tPNum", "tRep")
          #temp$randomization$tPNum <- paste(temp$randomization[,"Rep"],1:nrow(tempComb), sep = "")
          temp$randomization$tPNum <- 1:nrow(tempComb)
          tempBook <- merge(temp$randomization, tempBook, by.x = c("Rep","Treatment", "tPNum"), by.y = c("tRep","Treatment", "tPNum"))
          tempBook$FieldRow <- as.numeric(tempBook$FieldRow)
          tempBook$FieldCol <- as.numeric(tempBook$FieldCol)
          tempBook <- tempBook[order(tempBook$FieldRow, tempBook$FieldCol),]
          randomize <- rbind(randomize, cbind(Trial = i, tempBook[,c("Rep", "Block", "Treatment", "PlotNum", "FieldRow", "FieldCol")]))
          dimnames(plan[[i]]) <- list(paste("FieldRow", 1:nrow(plan[[i]]),sep = ""), paste("FieldCol", 1:ncol(plan[[i]]), sep = ""))
          
	}
	dimnames(plotNum) <- dimnames(plan[[1]])
	randomize[,"Trial"] <- factor(randomize[,"Trial"])
     rownames(randomize) <- 1:nrow(randomize)
     
	
	cat(toupper("Design Properties:"),"\n",sep = "")
	cat("\t","Incomplete Block Design","\n",sep = "") 
	#if (r == (sqrt(nrow(tempComb)) + 1)) { cat("\t","Balanced Lattice Design","\n\n",sep = "") } else { cat("\t","Partially Balanced Lattice Design","\n\n",sep = "") } 
     cat("\t", temp[[1]],"\n\n", sep = "")
	cat(toupper("Design Parameters:"),"\n",sep = "")
	cat("\t","Number of Trials = ", trial, "\n",sep = "")
	cat("\t","Number of Treatments = ", nrow(tempComb), "\n",sep = "")
	cat("\t","Number of Replicates = ", r, "\n",sep = "")
	cat("\t","Plots per Block (Block Size) = ", sqrt(nrow(tempComb)), "\n",sep = "")
	cat("\t","Block per Replicate = ", sqrt(nrow(tempComb)), "\n",sep = "")
	cat("\t","Number of Field Row = ", nrow(plotNum), "\n",sep = "")
	cat("\t","Number of Field Column = ", ncol(plotNum), "\n",sep = "")
     
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
		cat("\nResults of Randomization:\n")
		printDataFrame(randomize)
	}
	return(invisible(list(fieldbook = randomize, layout = plan, plotNum = plotNum)))
}

