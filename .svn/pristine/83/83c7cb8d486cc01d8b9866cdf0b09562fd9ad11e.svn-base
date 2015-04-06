#-----------------------------------------
# Created by: Alaine Gulles
#-----------------------------------------

class.information2 <-function(x, data) UseMethod("class.information2")

class.information2.default <-function(x, data) {
	data.info <- NULL
	for (i in (1:length(x))) {
		num.levels <- nlevels(data[,match(x[i], names(data))])
		if (num.levels > 0) {
			factor.levels <- levels(data[,match(x[i], names(data))])
			if (num.levels > 7) list.factor <- paste(paste(factor.levels[1:3], collapse = "  "), " ... ", factor.levels[num.levels], sep = "")
			else list.factor <- paste(factor.levels, collapse = "  ", sep = "")
			newdata<-c(x[i], num.levels, list.factor)
			data.info <- as.table(rbind(data.info, newdata))
		}
	}
	colnames(data.info) <- c("Factors  ", "No of Levels  ", "Levels")
	rnames<-rep("",nrow(data.info))
	rownames(data.info) <- rnames
	return(data.info)
}

