# -------------------------------------------------------------------------------------
# designAugmented: Generate randomization for augmented design in RCBD and LSD.
# Created by: Alaine A. Gulles 09.21.2010 for International Rice Research Institute
# Modified by: Alaine A. Gulles 04.30.2012
# -------------------------------------------------------------------------------------

designAugmented <- function(checkTrmt, newTrmt, r = NULL, trial = 1, design = c("rcbd", "lsd"), factorName = NULL, display = TRUE, file = NULL) UseMethod("designAugmented")

designAugmented.default <- function(checkTrmt, newTrmt, r = NULL, trial = 1, design = c("rcbd", "lsd"), factorName = NULL, display = TRUE, file = NULL) {

	design <- match.arg(design)
	if (length(checkTrmt) == 1) { checkTrmt <- paste("check",1:checkTrmt, sep = "") }
	if (length(newTrmt) == 1)   { newTrmt <- paste("new",1:newTrmt, sep = "") }
	#set.seed(0)
	if (design == "rcbd") { 
          if (length(newTrmt)%%r != 0) { stop("The number of new treatment should be divisible by the r.") }
          new <- paste("*", letters[1:(length(newTrmt)/r)], sep = "") 
          tempTrmt <- c(checkTrmt, new)
          capture.output(randomize <- designRCBD(list(tempTrmt = tempTrmt),r, trial, display = FALSE))
          randomize[,"tempTrmt"] <- as.character(randomize[,"tempTrmt"])
          tempNew <- sample(newTrmt, length(newTrmt), replace = FALSE)
          if (trial > 1) {
               for (i in (2:trial)) { tempNew <- c(tempNew, sample(newTrmt, length(newTrmt), replace = FALSE))	}
          }
          randomize[randomize[,"tempTrmt"]%in%as.character(new),"tempTrmt"] <- tempNew
                       
		if (is.null(factorName) && !is.valid.name(trimStrings(factorName[1]))) { colnames(randomize)[3] <- "Treatment"
		} else { colnames(randomize)[3] <- factorName[1] }
	} ## design == "rcbd"

	if (design == "lsd")  { 
          if (length(newTrmt)%%length(checkTrmt)**2 != 0) stop("The number of new treatments should be divisble by the number of check treatment squared.")
          addItem <- length(newTrmt)/length(checkTrmt)**2
          randomize <- NULL
		for (z in (1:trial)) {
               tempNewDesign <- NULL
               capture.output(tempDesign <- designLSD(list(tempTrmt = checkTrmt), trial = 1, display = FALSE))
               tempSample <- matrix(tempDesign[,"tempTrmt"], nrow = length(checkTrmt), ncol = length(checkTrmt), byrow = TRUE)          
               for (i in (1:length(checkTrmt))) {
                    tempRow <- matrix(0, nrow = (addItem + 1), ncol = length(checkTrmt))
                    for (j in (1:length(checkTrmt))) {
                         tempRow[sample(1:(addItem+1),1),j] <- tempSample[i,j]
                    }
                    tempNewDesign <- rbind(tempNewDesign, tempRow)
               }
               tempNewDesign[tempNewDesign == "0"] <- sample(newTrmt, length(newTrmt))
               tempNewDesign <- as.data.frame.table(tempNewDesign)
               tempNewDesign <- cbind(tempNewDesign[,"Var1"],tempNewDesign) 
               #tempNewDesign[,1] <- 1:((addItem+1)*length(checkTrmt))
               tempNewDesign[,1] <- rep(1:length(checkTrmt), each = (addItem + 1))
               tempNewDesign[,2] <- rep(1:length(checkTrmt), each = ((addItem+1)*length(checkTrmt)))
               tempNewDesign[,3] <- 1:(addItem+1)
               randomize <- rbind(randomize, data.frame(Trial = z, tempNewDesign, PlotNum = 1:nrow(tempNewDesign)))
		}
		
		if (is.null(factorName) && !is.valid.name(trimStrings(factorName[1]))) { colnames(randomize) <- c("Trial", "Row", "Column", "withRxC", "Treatment", "PlotNum") 
		} else { colnames(randomize) <- c("Trial", "Row", "Column", "withRxC", factorName, "PlotNum")  }
	}

	cat(toupper("Design Properties:"),"\n",sep = "")
	if (design == "rcbd") { cat("\t","Augmented Randomized Complete Block Design (Augmented RCBD)","\n\n",sep = "") 
	} else { cat("\t","Augmented Latin Square Design (Augmented LSD)","\n\n",sep = "") }
	cat(toupper("Design Parameters:"),"\n",sep = "")
	cat("\t","Number of Trials = ", trial, "\n",sep = "")
	cat("\t","Number of Replicated Treatments = ", length(checkTrmt), "\n",sep = "")
	cat("\t","Levels of Replicated Treatments = ", sep = "")
	if (length(checkTrmt) <= 5) { cat(paste(checkTrmt, collapse = ", ", sep = ""), "\n", sep = "")
	} else {
		cat(paste(checkTrmt[1:3],collapse = ", ", sep = ""), sep = "")
		cat(", ..., ",checkTrmt[length(checkTrmt)], "\n",sep = "")
	}
	if (design == "rcbd") cat("\t","Number of Replicates = ", r, "\n",sep = "")
	cat("\t","Number of Unreplicated Treatments = ", length(newTrmt), "\n",sep = "")
	cat("\t","Levels of UnReplicated Treatments = ", sep = "")
	if (length(newTrmt) <= 5) { cat(paste(newTrmt, collapse = ", ", sep = ""), "\n\n", sep = "")
	} else {
		cat(paste(newTrmt[1:3],collapse = ", ", sep = ""), sep = "")
		cat(", ..., ",newTrmt[length(newTrmt)], "\n\n",sep = "")
	}
	#cat("Results of Randomization:\n")
	#printDataFrame(randomize)

	if (!is.null(file)) {
		tempFile <- strsplit(file, split = "\\.")[[1]]
		tempExt <- tolower(tempFile[length(tempFile)])
		if (tempExt != "csv"){ if(tempExt != "rda") tempExt <- "csv" } 
		newFile <- paste(tempFile[1:length(tempFile)-1], collapse = ".", sep = "")
		newFile <- paste(newFile, tempExt, sep = ".")
		if (tempExt == "csv") { write.csv(randomize, file = newFile, row.names = FALSE, append = FALSE)
		} else { save(randomize, file = newFile) }
	} else {
		cat("Results of Randomization:\n")
		printDataFrame(randomize)
	}
	return(invisible(randomize))
}