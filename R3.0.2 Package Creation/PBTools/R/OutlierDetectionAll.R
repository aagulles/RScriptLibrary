# -------------------------------------------------------------------------
# Functions Included: OutlierDetection, OutlierDetectionMethod1,
#                     OutlierDetectionMethod2, OutlierGraph2,
#                     OutlierGraph2Sub
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# OutlierDetection
# Created by: Alaine A. Gulles for International Rice Research Institute
# Description: Main function for detecting outliers
# Modified by: Alaine Gulles for IRRI 11.24.2014
# -------------------------------------------------------------------------

OutlierDetection <- function(data, var, grp = NULL, rep = NULL, path = NULL, method = c("method1", "method2")) UseMethod("OutlierDetection")

OutlierDetection.default <- function(data, var, grp = NULL, rep = NULL, path = NULL, method = c("method1", "method2")) {
     if (is.character(data)) {
          nameData <- data
          if (!exists(data)) { stop(paste("The object '", nameData, "' not found.",sep ="")) }
          data <- eval(parse(text = data))
     } else { nameData <- paste(deparse(substitute(data))) }
     
     method <- match.arg(method) 
     if (is.null(path)) {  path <- getwd() }
     
     # check entry
     if (any(is.na(match(var, names(data))))) { stop("At least one of item in the character vector 'var' does not match any variable name in the dataset.") }
     if (!is.null(grp)) {
          if (any(is.na(match(grp, names(data))))) { stop("The item in the character vector 'grp' does not match any variable name in the dataset.") }
     }
     if (!is.null(rep)) {
          if (any(is.na(match(rep, names(data))))) { stop("The item in the character vector 'rep' does not match any variable name in the dataset.") }
     }
     
     # method 1:
     if (method == "method1") {
          result <- list()
          grp <- grp[1]
          if (is.null(grp) || length(grp) == 1) {
               for (i in 1:length(var)) {
                    result[[i]] <- OutlierDetectionMethod1(data, var[i], grp)
                    out <- result[[i]]
                    result[[i]][1] <- var[i]
                    result[[i]][2:3] <- out
                    names(result[[i]]) <- c("respvar", names(out))
               }     
          } 
          return(result)
     } ## end -- if (method == "method1")
     
     # method 2:
     if (method == "method2") {
          if (is.null(grp)) { stop("The parameter 'grp' is required.") }
          if (is.null(rep)) { stop("The parameter 'rep' is required.") }
          grp <- grp[1]
          result <- list()
          for (i in (1:length(var))) {
               result[[i]] <- OutlierDetectionMethod2(data = data, var = var[i], grp)
               out <- result[[i]] 
               OutlierGraph2(data = out[[1]], var = var[i], grp, rep, path)
               result[[i]][1] <- var[i]
               result[[i]][2:3] <- out
               names(result[[i]]) <- c("respvar", names(out))
          }
          return(result)
     } ## end stmt -- if (method == "method2") 
}


# -------------------------------------------------------------------------
# OutlierDetectionMethod1
# Created by: Alaine A. Gulles for International Rice Research Institute
# Description: Detect outliers using the IQR of the raw data
# Modified by: Alaine Gulles for IRRI 11.24.2014
# -------------------------------------------------------------------------

OutlierDetectionMethod1 <- function(data, var, grp = NULL) {
     
     outliers <- NULL
     capture.output(statResult <- DescriptiveStatistics(data, var, grp, statistics = c("min", "max", "iqr", "q1", "q3")))
     
     # to detection possible outlier
     IQR15 = statResult[,"IQR"]*1.5
     LLIQR15 = ifelse(statResult[,"Q1"] - IQR15 < statResult[,"Min"], statResult[,"Min"], statResult[,"Q1"] - IQR15)
     ULIQR15 = ifelse(statResult[,"Q3"] + IQR15 > statResult[,"Max"], statResult[,"Max"], statResult[,"Q3"] + IQR15)
     
     # to detect outlier
     IQR3 = statResult[,"IQR"]*3
     LLIQR3 = ifelse(statResult[,"Q1"] - IQR3 < statResult[,"Min"], statResult[,"Min"], statResult[,"Q1"] - IQR3)
     ULIQR3 = ifelse(statResult[,"Q3"] + IQR3 > statResult[,"Max"], statResult[,"Max"], statResult[,"Q3"] + IQR3)
     
     # combine Dataset
     statResult <- cbind(statResult, IQR15, LLIQR15, ULIQR15, IQR3, LLIQR3, ULIQR3)
     
     MergeData <- merge(data, statResult)
     MergeData$Variable <- NULL
     
     MergeData$OutlierCode <- 0
     MergeData$OutlierCode <- ifelse(MergeData[,"ULIQR15"] < MergeData[,var], 
                                     ifelse(MergeData[,"ULIQR3"] < MergeData[,var], 2, 1), 
                                     ifelse(MergeData[,"LLIQR15"] > MergeData[,var], 
                                            ifelse(MergeData[,"LLIQR3"] > MergeData[,var], 2,1),0))
     
     if(any(MergeData[,"OutlierCode"] > 0)) { 
          outliers <- MergeData[MergeData[,"OutlierCode"] > 0, c(names(data), "OutlierCode")] 
     } 
     
     newData <- MergeData[,c(names(data), "OutlierCode")]
     
     return(list(data = newData, outlier = outliers))
}

# -------------------------------------------------------------------------
# OutlierDetectionMethod2
# Script by: Violeta I. Bartolome for International Rice Research Institute
# Function Created by: Alaine Gulles for IRRI
# Modified by: Alaine Gulles for IRRI 11.24.2014
# Description: Function was develop for BIMMS. Detecting outliers based
#              from the IQR of the range per treatment/genotype
# -------------------------------------------------------------------------

OutlierDetectionMethod2 <- function(data, var, grp) {
     
     outliers <- NULL
     capture.output(suppressWarnings(statResultGrp <- DescriptiveStatistics(data, var, grp, statistics = c("nnmiss","min", "max", "range", "mean"))))
     statResultGrpNew <- statResultGrp[statResultGrp[,"N_NonMissObs"] >= 1, c("Variable", grp, "Min", "Max", "Mean", "Range")]
     capture.output(statResultRange <- DescriptiveStatistics(data = statResultGrpNew, var = "Range", statistics = c("min", "max", "iqr", "q1", "q3")))
     names(statResultRange)[2:ncol(statResultRange)] <- paste(names(statResultRange)[2:ncol(statResultRange)], paste("Range"), sep = "_")
     
     # to detection possible outlier
     IQR15 = statResultRange[,"IQR_Range"]*1.5
     LLIQR15 = ifelse(statResultRange[,"Q1_Range"] - IQR15 < statResultRange[,"Min_Range"], statResultRange[,"Min_Range"], statResultRange[,"Q1_Range"] - IQR15)
     ULIQR15 = ifelse(statResultRange[,"Q3_Range"] + IQR15 > statResultRange[,"Max_Range"], statResultRange[,"Max_Range"], statResultRange[,"Q3_Range"] + IQR15)
     
     # to detect outlier
     IQR3 = statResultRange[,"IQR_Range"]*3
     LLIQR3 = ifelse(statResultRange[,"Q1_Range"] - IQR3 < statResultRange[,"Min_Range"], statResultRange[,"Min_Range"], statResultRange[,"Q1_Range"] - IQR3)
     ULIQR3 = ifelse(statResultRange[,"Q3_Range"] + IQR3 > statResultRange[,"Max_Range"], statResultRange[,"Max_Range"], statResultRange[,"Q3_Range"] + IQR3)
     
     #statResultGrp$OutlierCode <- ifelse(statResultGrp[,"Range"] > ULIQR15,
     #                                    ifelse(statResultGrp[,"Range"] > ULIQR3,2,1),
     #                                   ifelse(statResultGrp[,"Range"] < LLIQR15,
     #                                          ifelse(statResultGrp[,"Range"] < LLIQR3,2,1),0))
     
     statResultGrp$OutlierCode <- ifelse(statResultGrp[,"N_NonMissObs"] > 0, 
                                         ifelse(statResultGrp[,"Range"] > ULIQR15,
                                                ifelse(statResultGrp[,"Range"] > ULIQR3,2,1),
                                                ifelse(statResultGrp[,"Range"] < LLIQR15,
                                                       ifelse(statResultGrp[,"Range"] < LLIQR3,2,1),0)), NA)
     
     statResultGrp$Variable <- NULL
     MergeData <- merge(data, statResultGrp)
          
     if(any(statResultGrp[,"OutlierCode"] > 0)) {
          tmp1 <- MergeData[!is.na(MergeData[,"OutlierCode"]),]
          outliers <- tmp1[tmp1[,"OutlierCode"] > 0, c(names(data), "OutlierCode")]
          # outliers <- MergeData[MergeData[,"OutlierCode"] > 0, c(names(data), "OutlierCode")] 
     } 
     return(list(data = MergeData, outlier = outliers))
}

# -------------------------------------------------------------------------
# OutlierGraph2
# Script by: Violeta I. Bartolome for International Rice Research Institute
# Function Created by: Alaine Gulles for IRRI
# Modified by: Alaine Gulles for IRRI 11.24.2014
# -------------------------------------------------------------------------

OutlierGraph2 <- function(data, var, grp, rep, path) {
     
     capture.output(statData <- DescriptiveStatistics(data, var, statistics = c("min", "max", "mean")))
     capture.output(statDataGrp <- DescriptiveStatistics(data, "Mean", grp, statistics = c("mean")))
     statDataGrp <- statDataGrp[!is.nan(statDataGrp[,"Mean"]),]
     statDataGrp <- statDataGrp[order(-statDataGrp[,"Mean"]),]
     tempName <- make.unique(c(names(data), "MeanRank"))[length(make.unique(c(names(data), "MeanRank")))]
     statDataGrp[,tempName] <- 1:nrow(statDataGrp)
     statDataGrp <- statDataGrp[,-I(match(c("Variable", "Mean"), names(statDataGrp)))]

     tempData <- merge(data, statDataGrp, by = grp)
     tempData <- tempData[!is.na(tempData[,var]),]
     tempGrpName <- make.unique(c(names(tempData), paste(grp,"New", sep = "")))[length(make.unique(c(names(tempData), paste(grp,"New", sep = ""))))]
     tempData[,tempGrpName] <- trimStrings(as.character(tempData[,grp]))
     tempData[,tempGrpName] <- paste(tempData[,tempGrpName], tempData[,"OutlierCode"], sep = " ")
     tempData <- tempData[order(tempData[,tempName]),]
     
     # determine how many png files will be created
     numLines <- 20
     numjpegfiles <- floor(nrow(statDataGrp)/numLines)
     pagesWithExtraLines <- nrow(statDataGrp) %% numLines
     #extraLines <- nrow(statDataGrp) - (numjpegfiles*numLines)
     LL <- 1
     graphNum <- 1
     while (LL <= nrow(statDataGrp)) {
          numEntries <- ifelse(graphNum <= pagesWithExtraLines ,numLines + 1, numLines)
          UL <- LL + numEntries - 1
          # UL <- ifelse((numEntries + LL - 1) > nrow(statDataGrp), nrow(statDataGrp), numEntries + LL - 1)
          tempData1 <- tempData[tempData[,tempName] >= LL & tempData[,tempName] <= UL,]
          png(filename = paste(path, "/",var,"_Graph", graphNum,".png", sep = ""), res = 150, height = 800, width = 1000)
          OutlierGraph2Sub(data = tempData1, var, grp = tempGrpName, ylabel = tempName, rep, 
                           min = statData[,"Min"], max = statData[,"Max"], mean = statData[,"Mean"])  
          dev.off()
          LL <- UL + 1
          graphNum <- graphNum + 1
     }
}

# -------------------------------------------------------------------------
# OutlierGraph2Sub
# Script by: Violeta I. Bartolome for International Rice Research Institute
# Function Created by: Alaine Gulles for IRRI
# Modified by: Alaine Gulles for IRRI 11.24.2014
# -------------------------------------------------------------------------

OutlierGraph2Sub <- function(data, var, grp, ylabel, rep, min, max, mean) {
     
     marginNumLinesTop <- ifelse(nrow(data) < 100, 2, 0.5)
     printTop <- ifelse(nrow(data) < 100, 1, 0)
     
     graphLabel <- data[,c(grp, ylabel)]
     graphLabel <- graphLabel[!duplicated(data[,ylabel]),]
     
     maxCharLabel <- max(nchar(graphLabel[,grp]))
     marginNumLinesLeft <- as.numeric(ifelse(ceiling(maxCharLabel/10) <= 5,5,
                                             ifelse(ceiling(maxCharLabel/10) >= 6, 6, 6 - ((50 - maxCharLabel)/40))))
     
     marginSizeLeft <- ifelse(maxCharLabel <= 10, 1,
                              ifelse(maxCharLabel >= 50, 4, 4 - (3 * (50 - maxCharLabel))/40))
          
     par(mar = c(1,marginNumLinesLeft,marginNumLinesTop,1), mai = c(1,marginSizeLeft,1,1), cex.axis = 0.7)
     plot(data[,var], data[,ylabel], xlim = c(min, max), pch = " ", axes = F,
          ylim = rev(range(data[,ylabel])), xaxt = "n", yaxt = "n", ylab = " ", xlab = " ")
     LL <- data[1,ylabel]
     UL <- data[nrow(data), ylabel]
          
     for (j in seq((LL + 0.5), (UL - 0.5), 2)) rect(min, j, max, j + 1, col = 'gray', border = NA)
     
     points(data[,var], data[,ylabel], yaxt = "n", pch = 16, col = "black", cex = 2)
     text(data[,var], data[,ylabel], label = data[,rep], col='white',cex=.6)
     axis(3, at = c(min, mean, max), pos = LL-printTop)
     axis(2, at = seq(LL,UL,1), labels=graphLabel[,grp], las = 1, tick = FALSE, hadj = 1, padj = 1)
     segments(x0 = mean, y0 = (LL - 1), x1 = mean, y1 = UL+1, lty = 2, col = "red")
     
}
