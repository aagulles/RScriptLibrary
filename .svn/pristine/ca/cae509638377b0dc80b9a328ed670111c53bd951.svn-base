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
	set.seed(0)
	if (design == "rcbd") { 
		new <- paste("*", letters[1:ceiling(length(newTrmt)/r)], sep = "") 
		tempTrmt <- c(checkTrmt, new)
		capture.output(randomize <- designRCBD(list(tempTrmt = tempTrmt),r, trial, display = FALSE))
		randomize[,"tempTrmt"] <- as.character(randomize[,"tempTrmt"])
		if (is.wholenumber(length(newTrmt)/r)) {
			tempNew <- sample(newTrmt, length(newTrmt), replace = FALSE)
			if (trial > 1) {
				for (i in (2:trial)) { tempNew <- c(tempNew, sample(newTrmt, length(newTrmt), replace = FALSE))	}
			}
			randomize[randomize[,"tempTrmt"]%in%as.character(new),"tempTrmt"] <- tempNew
		} else {
			addItem <- (ceiling(length(newTrmt)/r)*r) - length(newTrmt)
			tempNew <- c(sample(newTrmt, length(newTrmt), replace = FALSE), rep("", addItem))
			if (trial > 1) {
				for (i in (2:trial)) { tempNew <- c(tempNew, c(sample(newTrmt, length(newTrmt), replace = FALSE), rep("", addItem))) }
			}
			randomize[randomize[,"tempTrmt"]%in%as.character(new),"tempTrmt"] <- tempNew
			randomize <- subset(randomize, tempTrmt != "")
			randomize[,"PlotNum"] <- 1:(length(newTrmt) + length(checkTrmt)*r)
			rownames(randomize) <- 1:nrow(randomize)
		}
		if (is.null(factorName) && !is.valid.name(trim.strings(factorName[1]))) { colnames(randomize)[3] <- "Treatment"
		} else { colnames(randomize)[3] <- factorName[1] }
	} ## design == "rcbd"

	if (design == "lsd")  { 
		addItem <- floor(length(newTrmt)/length(checkTrmt))
		check <- TRUE
		moreItem <- NULL
		while (check) {
			totalItem <- (length(checkTrmt) + addItem)**2
			itemNeeded <- (length(checkTrmt)**2) + (length(checkTrmt) * addItem) + length(newTrmt)
			if (itemNeeded == totalItem) { check <- FALSE
			} else {
				if (itemNeeded < totalItem) {
					if (is.null(moreItem)) { moreItem <- "totalItem"	
					} else { if (moreItem != "totalItem") { check <- FALSE }}
					addItem <- addItem - 1
				} else {
					if (is.null(moreItem)) { moreItem <- "itemNeeded"	
					} else { if (moreItem != "itemNeeded") { check <- FALSE }}
					addItem <- addItem + 1
				}
			}
			if (addItem < 1) { addItem <- 1; check <- FALSE }
		}	

		new <- paste("*", letters[1:addItem], sep = "") 
		tempTrmt <- c(checkTrmt, new)
		capture.output(randomize <- designLSD(list(tempTrmt = tempTrmt), trial, display = FALSE))
		randomize[,"tempTrmt"] <- as.character(randomize[,"tempTrmt"])
		if (itemNeeded == totalItem) {
			tempNew <- c(sample(newTrmt, length(newTrmt), replace = FALSE))
			if (trial > 1) {
				for (i in (2:trial)) {
					tempNew <- c(tempNew, sample(newTrmt, length(newTrmt), replace = FALSE))
				}
			}
			randomize[randomize[,"tempTrmt"]%in%as.character(new),"tempTrmt"] <- tempNew

		} else {
			tempNew <- c(sample(newTrmt, length(newTrmt), replace = FALSE), rep("", (totalItem - itemNeeded)))
			if (trial > 1) {
				for (i in (2:trial)) {
					tempNew <- c(tempNew, c(sample(newTrmt, length(newTrmt), replace = FALSE), rep("", (totalItem - itemNeeded))))

				}
			}
			randomize[randomize[,"tempTrmt"]%in%as.character(new),"tempTrmt"] <- tempNew
			randomize <- subset(randomize, tempTrmt != "")
			randomize[,"PlotNum"] <- 1:itemNeeded
			rownames(randomize) <- 1:nrow(randomize)
		}
		if (is.null(factorName) && !is.valid.name(trim.strings(factorName[1]))) { colnames(randomize)[4] <- "Treatment"
		} else { colnames(randomize)[4] <- factorName[1] }
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