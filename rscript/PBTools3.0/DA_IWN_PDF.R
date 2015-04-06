# ----------------------------------------------------------------------------------------------------
# RCropStat Beta Version: Function
# ----------------------------------------------------------------------------------------------------
# DataAttribute
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles 05.02.2012
# ----------------------------------------------------------------------------------------------------
#
#' Function for displaying data attributes
#' @param data dataframe containing the data
#'
#' @return table tempTable

DataAttribute <- function(data) {
  if(is.character(data)) { 
    nameData <- data
    if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
    tempData <- eval(parse(text = data))
  } else { 
    nameData <- paste(deparse(substitute(data)))	
    #if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
    tempData <- data
  }
  if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
  tempTable <- data.frame(names(tempData))
  withFactor <- FALSE
  for (i in (1:nrow(tempTable))) {
    if (is.factor(tempData[,i])) {
      withFactor <- TRUE
      if (is.ordered(tempData[,i])) tempTable[i,2] <- "ordered factor"
      else tempTable[i,2] <- "factor"
      tempTable[i,3] <- as.character(nlevels(tempData[,i]))
      if (nlevels(tempData[,i]) > 5) tempTable[i,4] <- paste(c(paste(levels(tempData[,i])[1:2], collapse = ", ", sep = ""),levels(tempData[,i])[nlevels(tempData[,i])]), collapse = ", ..., ", sep = "")
      else tempTable[i,4] <- paste(levels(tempData[,i]), collapse = ", ", sep = "")
    } else {
      if (typeof(tempData[,i]) == "double") { tempTable[i,2] <- "numeric" }
      else tempTable[i,2] <- typeof(tempData[,i])
      tempTable[i,3] <- ""
      tempTable[i,4] <- ""
    }
  }
  names(tempTable) <- c("VAR NAME", "TYPE", "NLEVELS", "LEVELS")
  if (!withFactor) { tempTable <- tempTable[,1:2] }
  return(tempTable)
}


# -------------------------------------
# RCropStat Utilities: Function
# -------------------------------------
# isWholenumber/is.wholenumber:
# Reference: FROM FULLREFMAN.PDF
# -------------------------------------
#
#' Function for determining whether input is a whole number
#' @param x input
#' @param tol cutoff for absolute difference between x and round(x)

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) { 
  abs(x - round(x)) < tol 
}


# RCropStat Utilities

# -----------------------------------------------------------------------------------
# printDataFrame: Print the dataframe
# Created by: Alaine A. Gulles 04.30.2012 for International Rice Research Institute
# Modified by: Alaine A. Gulles 07.10.2012
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# ARGUMENTS:
# dataFrame = the data frame to be printed
# -----------------------------------------------------------------------------------
#
#' Function for formatted printing of a data frame
#' @param dataFrame data to be printed
#' @param border logical; whether a border is to be printed or not
#' @param digits number of digits used in formatting contents of the table contents

printDataFrame <- function(dataFrame, border = TRUE, digits = NULL) {
  
  if (!is.data.frame(dataFrame)) { stop("The argument 'dataFrame' should be of class data.frame.") }
  dataChar <- DataAttribute(dataFrame)[,1:2]
  dataChar[,1] <- as.character(dataChar[,1])
  
  colWidth <- NULL
  onlyDecimal <- NULL
  for (i in (1:nrow(dataChar))) {
    if (dataChar[i,2] == "numeric") {
      if (all(is.na(dataFrame[,i]))) { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 1)
      } else {
        if (all(dataFrame[,i] < 1) && all(dataFrame[,i] > -1)) { 
          colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 5) + 1)
          onlyDecimal <- c(onlyDecimal, i)
        } else {
          if (any(is.na(is.wholenumber(dataFrame[,i])))) {
            tempsubdata <- subset(dataFrame[,i], dataFrame[,i] != Inf)
            if (nrow(tempsubdata) == 0 || is.null(nrow(tempsubdata))) {
              colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), 3) + 1) 
            } else {
              if (all(is.wholenumber(tempsubdata))) { 
                dataChar[i,2] <- "integer" 
                colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), nchar(max(tempsubdata))) + 1)
              } else {
                if (is.null(digits)) { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(tempsubdata,0))) + 3) + 1)
                } else { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(tempsubdata,0))) + digits) + 1) }
              }
            }
            
          } else {
            if (all(is.wholenumber(dataFrame[,i]))) { 
              dataChar[i,2] <- "integer" 
              colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), nchar(max(dataFrame[,i]))) + 1)
            } else { 
              if (is.null(digits)) { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + 3) + 1)
              } else { colWidth <- c(colWidth, max(3, nchar(as.character(dataChar[i,1])), max(nchar(round(dataFrame[,i],0))) + digits) + 1) }
            }
          }
        }	
      }
    } else {
      if (dataChar[i,2] == "integer") { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), nchar(max(dataFrame[,i]))) + 1)
      } else { colWidth <- c(colWidth, max(nchar(as.character(dataChar[i,1])), max(nchar(as.character(dataFrame[,i])))) + 1) }
    }
    
  }
  
  if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")	
  
  for (j in (1:ncol(dataFrame))) {
    if (dataChar[j,2] == "numeric" || dataChar[j,2] == "integer") { 
      cat(formatC(dataChar[j,1], width = colWidth[j], format = "s"), sep = "")
    } else { cat(formatC(dataChar[j,1], width = colWidth[j], format = "s", flag = "-"), sep = "") }
    if (j == ncol(dataFrame)) { 
      cat("\n")
      #if (nrow(dataFrame) >= 2) cat("\n\n") else cat("\n")
    } else { cat(formatC("", width = 1, format = "s")) }
  }
  
  if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")	
  
  for (i in (1:nrow(dataFrame))) {
    for (j in (1:ncol(dataFrame))) {
      if (dataChar[j,2] == "numeric") {
        if (is.na(match(j, onlyDecimal))) { 
          if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") 
          } else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = 2), sep = "") }
        } else { 
          if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") 
          } else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "f", digits = 4), sep = "") }
        }
      } else {
        if (dataChar[j,2] == "integer") { 
          if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s"), sep = "") } 
          else { cat(formatC(dataFrame[i,j], width = colWidth[j], format = "d"), sep = "") }
        } else { 
          if (is.na(dataFrame[i,j])) { cat(formatC("", width = colWidth[j], format = "s", flag = "-"), sep = "") }
          else { cat(formatC(as.character(dataFrame[i,j]), width = colWidth[j], format = "s", flag = "-"), sep = "") }	
        }
      }
      if (j == ncol(dataFrame)) { cat("\n") } else { cat(formatC("", width = 1, format = "s")) }
    }
  }
  if (border) cat(formatC(paste(rep("-",sum(colWidth) + nrow(dataChar) - 1), collapse = ""), width = sum(colWidth) + nrow(dataChar) - 1, format = "s"), "\n")	
  
}
