# ----------------------------------------------------------------------
# estimateMissingData: Function for Estimating Missing Data using
# iterative fitting of regression models 
# Created by: Dr. Guoyou Ye
#		          Alaine A. Gulles
#             Nellwyn L. Sales 24.07.2012
# ----------------------------------------------------------------------

estimateMissingData <- function(design=c("CRD", "RCB"), data=NULL , respvar, male, female, rep) UseMethod("estimateMissingData")

estimateMissingData.default <- function(design=c("CRD", "RCB"), data=NULL , respvar, male, female, rep) {

	tempFullData <- data
	tempCompleteData <- tempFullData [complete.cases(tempFullData [, respvar]),]

	# --- Generate list of crosses with at least one observation
	levelsList<-list(cross=levels(factor(tempCompleteData$Cross)),rep=levels(factor(tempCompleteData[,match(rep, names(tempCompleteData ))])))
	crossNotMissing<-fac.gen(levelsList)
	colnames(crossNotMissing)<-c("Cross", rep) 
	fullCross<-merge(crossNotMissing,tempCompleteData ,by=intersect(names(crossNotMissing),names(tempCompleteData)),all=T,order=T)

	# --- Check if there are missing reps. If there are, generate estimates for these missing values --- #
	if (nrow(tempCompleteData)<nrow(fullCross)) {

		# --- Substitute the rep means to the missing reps --- #
		repMeans <- data.frame(levels(fullCross[,rep]), tapply(fullCross[,respvar], fullCross[,rep], mean, na.rm = TRUE))
		colnames(repMeans) <- c(rep, paste(respvar, "means", sep = ""))
		remark<-NULL
		newData <- data.frame(merge(fullCross, repMeans, by = rep), remark = "obs")
		newData$remark <- as.character(newData$remark)
		newData[is.na(newData[,respvar]),"remark"] <- "estimate"
		newData[is.na(newData[,respvar]),respvar] <- newData[is.na(newData[,respvar]),paste(respvar, "means", sep = "")]
		newData<- newData[,-I(match(paste(respvar, "means", sep = ""),names(newData)))]

		# --- Start of iteration --- #	
		if (design == "RCB") { myformula <- paste(respvar, " ~ ", rep, " + Cross", sep = "")
		} else { myformula <- paste(respvar, " ~ Cross", sep = "") }

		stable <- FALSE
		iterationNumber<-1
		newData$sumEstimates<-0
		while(!stable & iterationNumber<=100) {
			result <- aov(formula(myformula), data = newData)
			newData$predval <- fitted.values(result)
			estimatedData <- subset(newData, remark == "estimate")
			if (all(abs(estimatedData[,respvar]-estimatedData[,"predval"])/estimatedData[,respvar] < 0.001)) { 
				stable <- TRUE
				newData[newData[,"remark"] == "estimate",respvar] <- newData[newData[,"remark"] == "estimate","predval"]
				newData <- newData[,-I(match("predval", names(newData)))]
				newData <- newData[,-I(match("sumEstimates", names(newData)))]
			} else {
				newData[newData[,"remark"] == "estimate",respvar] <- newData[newData[,"remark"] == "estimate","predval"]
				newData$sumEstimates<-newData$sumEstimates + newData$predval
				newData <- newData[,-I(match("predval", names(newData)))]
				iterationNumber<-iterationNumber + 1
			}	
		}
		if (stable==FALSE){
			iterationNumber<-iterationNumber-1
			newData$sumEstimates<-newData$sumEstimates/iterationNumber
			newData[newData[,"remark"] == "estimate",respvar] <- newData[newData[,"remark"] == "estimate","sumEstimates"]
			newData <- newData[,-I(match("sumEstimates", names(newData)))]
		}
		newData<-newData[,c("Cross", rep, respvar)]
		tempFullData<-tempFullData[,-I(match(respvar, names(tempFullData)))]
		tempFullData<-merge(tempFullData,newData,by=c("Cross",rep),all=T)
		tempCompleteData <- tempFullData [complete.cases(tempFullData [, respvar]),]
	}

	# --- Check if there are missing treatment combinations. If there are, generate estimates for these missing values --- #
	if (nrow(tempCompleteData)<nrow(tempFullData)) {
		
		tempMissingData <- tempFullData[!complete.cases(tempFullData[, respvar]),]
		# --- Substitute the means of crosses where each parent is involved to the missing values --- #
		for (o in (1:nrow(tempMissingData))) {
			#tempRepSubset <- subset(tempCompleteData, tempCompleteData[,match(rep, names(tempCompleteData))] == tempMissingData[o,rep])
			#tempDataSubset <- subset(tempRepSubset, tempRepSubset[,match(male, names(tempRepSubset))] == tempMissingData[o,male] | tempRepSubset[,match(female, names(tempRepSubset))] == tempMissingData[o,female])
			tempRepSubset <- tempCompleteData[tempCompleteData[,match(rep, names(tempCompleteData))] == tempMissingData[o,rep],]
      tempDataSubset <- tempRepSubset[(tempRepSubset[,match(male, names(tempRepSubset))] == tempMissingData[o,male] | tempRepSubset[,match(female, names(tempRepSubset))] == tempMissingData[o,female])==TRUE,]
      tempMissingData[o,respvar] <- mean(tempDataSubset[,match(respvar,names(tempDataSubset))])
		}
		colnames(tempMissingData)[match(respvar,names(tempMissingData))]<-"SubsetMeans"
		newData<-data.frame(merge(tempFullData, tempMissingData, by=c("Cross", male, female, rep), all=TRUE),remark="obs")
		newData$remark <- as.character(newData$remark)
		newData[is.na(newData[,respvar]),"remark"] <- "estimate"
		newData[is.na(newData[,respvar]),respvar] <- newData[is.na(newData[,respvar]),"SubsetMeans"]
		newData<- newData[,-I(match("SubsetMeans",names(newData)))]

		# --- Start of iteration --- #	
		if (design == "RCB") { myformula <- paste(respvar, " ~ ", rep, " + Cross", sep = "")
		} else { myformula <- paste(respvar, " ~ Cross", sep = "") }

		stable <- FALSE
		while(!stable) {
			result <- aov(formula(myformula), data = newData)
			newData$predval <- fitted.values(result)
			estimatedData <- subset(newData, remark == "estimate")
			if (all(abs(estimatedData[,respvar]-estimatedData[,"predval"])/estimatedData[,respvar] < 0.001)) { 
				stable <- TRUE
				newData[newData[,"remark"] == "estimate",respvar] <- newData[newData[,"remark"] == "estimate","predval"]
				newData <- newData[,-I(match("predval", names(newData)))]
			} else {
				newData[newData[,"remark"] == "estimate",respvar] <- newData[newData[,"remark"] == "estimate","predval"]
				newData <- newData[,-I(match("predval", names(newData)))]
			}	
		}
		tempFullData<-newData
	}
	return(tempFullData)
} #end of function 