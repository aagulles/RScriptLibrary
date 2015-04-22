# ------------------------------------------------------------------
# displayErrMsg
# Description: display the error message
# Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------

displayErrMsg <- function(object) UseMethod("displayErrMsg")

displayErrMsg.default <- function(object) {
     if (class(object) != "try-error") { stop("Object should be of class 'try-error'.") }
     errMsg <- trimStrings(strsplit(object, ":")[[1]])
     errMsg <- trimStrings(paste(strsplit(errMsg, "\n")[[length(errMsg)]], collapse = " "))
     errMsg <- gsub("\"", "", errMsg)
     return(errMsg)
} ## end of function -- displayErrMsg