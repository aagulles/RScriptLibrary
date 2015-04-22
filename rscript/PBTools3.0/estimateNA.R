# ----------------------------------------------------------------------
# This function estimates missing values for alpha and row column designs for data without missing treatment combination
# Author: Nellwyn Sales
# ----------------------------------------------------------------------

estimateNA <- function(design=c("Alpha", "RowColumn"), fullData, respvar, parent1, parent2, rep, block=NULL, row=NULL, column=NULL) UseMethod("estimateNA")

estimateNA.default <- function(design=c("Alpha", "RowColumn"), fullData, respvar, parent1, parent2, rep, block=NULL, row=NULL, column=NULL) {
  
  # --- Substitute the rep means to the missing reps --- #
  repMeans <- data.frame(levels(fullData[,rep]), tapply(fullData[,respvar], fullData[,rep], mean, na.rm = TRUE))
  colnames(repMeans) <- c(rep, paste(respvar, "means", sep = ""))
  remark<-NULL
  newData <- data.frame(merge(fullData, repMeans, by = rep), remark = "obs")
  newData$remark <- as.character(newData$remark)
  newData[is.na(newData[,respvar]),"remark"] <- "estimate"
  newData[is.na(newData[,respvar]),respvar] <- newData[is.na(newData[,respvar]),paste(respvar, "means", sep = "")]
  newData<- newData[,-I(match(paste(respvar, "means", sep = ""),names(newData)))]
  
  # --- Start of iteration --- #
  newData$Cross <- newData[,match(parent1, names(newData))]:newData[,match(parent2, names(newData))]
  if (design == "Alpha") {myformula <- paste(respvar, " ~ Cross + ", rep, " + ", rep, ":", block, sep = "") }
  if (design == "RowColumn") {myformula <- paste(respvar, " ~ Cross + ", rep, " + ", row, ":", rep, " + ", column, ":", rep, sep = "") }

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
  
  #return data set with estimates of missing values
  if (design == "Alpha") {newData<-newData[,c(rep, block, parent1, parent2, respvar)] }
  if (design == "RowColumn") {newData<-newData[,c(rep, row, column, parent1, parent2, respvar)] }
                              
  return(newData)
} #end of function 
