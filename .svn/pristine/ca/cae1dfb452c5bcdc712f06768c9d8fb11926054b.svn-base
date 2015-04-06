###################################################################
#																  #
# 			R-CROPSTAT UTILITIES: STATISTICS					  #
#																  #
###################################################################
# DOCUMENTATION:												  #
# - LAST UPDATED: OCT 7, 2010									  #
###################################################################

skew <- function(x, na.rm = FALSE) {
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	m3 <- sum((x - mean(x))^3)/length(x)
	s3 <- sqrt(var(x))^3
	m3/s3
}

kurtosis.var <- function(x, na.rm = FALSE) {
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	m4 <- sum((x - mean(x))^4)/length(x)
	s4 <- var(x)^2
	(m4/s4) - 3
}

se.mean <- function(x, na.rm = FALSE) { 
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	sqrt(var(x)/length(x))
}

se.skewness <- function(x, na.rm = FALSE) { 
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	skew(x)/sqrt(6/length(x))
}

se.kurtosis <- function(x, na.rm = FALSE) { 
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	kurtosis(x)/sqrt(24/length(x))
}

mode.var <- function(x, na.rm = FALSE) {
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	f <- table(x)
	f <- sort(f)
	if (length(f) == 1) mode <- x[1]
	else {
		if (f[[1]] == f[[length(f)]]) mode <- ""
		else {
			i <- length(f) - 1
			mode <- as.character(rownames(f)[length(f)])
			while(f[[i]] == f[[length(f)]]) {
				mode <- c(mode,rownames(f)[i])
				i <- i-1
			}
			mode <- sort(mode)	
		}
	}
	return(mode)
}

coef.var <- function(x, na.rm = FALSE) { 
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	(sd(x)/mean(x))*100 
}

freq.var <- function(x) {
	f <- table(x)
	out <- data.frame(f)
	names(out)[1] <- "CATEGORY"
	return(out)
}

ucss <- function(x, na.rm = FALSE) { 
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	sum(x**2)	
}

css <- function(x, na.rm = FALSE) { 
	if (is.factor(x)) stop ("need numeric data")
	if (na.rm) x <- x[!is.na(x)] else if(any(is.na(x))) return(x[FALSE][NA])
	sum(x**2) - ((sum(x)**2)/length(x))
}


summaryStat <- function(data, grp = NULL, stat){

	avail.stat <-c("n", "nnmiss", "nmiss", "sum", "css", "ucss", "se.skew", "se.kurtosis", "range", "iqr", "var", "sd", "se.mean", 
		  "cv", "mean", "median", "mode", "min", "max", "q1", "q3", "skew", "kurtosis")
	stat.to.compute <- match(stat, avail.stat)
	statcode <- c(rep("0", each = 23)) 
	for (i in (1:23)) {
		for (j in (1:length(stat.to.compute))) {
			if (i == stat.to.compute[j]) statcode[i] <-"1"
		}
	}
	if (is.null(grp)) grp <- c(rep(1, each = nrow(data)))
	summaryTable <- NULL
	out.label <- NULL

	if(statcode[1] == "1") { 
		a <- NULL
		for (i in (1:ncol(data))) {	a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, length)))	}
		summaryTable <- a
		colnames(summaryTable)[ncol(summaryTable)] <- "Nobs"
		out.label <- c(out.label, "No. of Obs.")
	}

	if(statcode[2] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(apply(data[i], 2, is.na), grp, length)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "NNMissObs"
		out.label <- c(out.label,"No. of Non-Missing Obs.")
	}
	if(statcode[3] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, length) - tapply(apply(data[i], 2, is.na), grp, length)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "NMObs"
		out.label <- c(out.label,"No. of Missing Obs.")
	}
	if(statcode[18] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, min, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "min"
		out.label <- c(out.label,"Minimum")
	}
	if(statcode[19] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, max, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "max"
		out.label <- c(out.label,"Maximum")
	}
	if(statcode[4] == "1") {
		a <- NULL		
			for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, sum, na.rm = TRUE)))
			ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
			colnames(summaryTable)[ncol(summaryTable)] <- "sum"
			out.label <- c(out.label,"Sum")
	}
	if(statcode[15] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, mean, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "mean"
		out.label <- c(out.label,"Mean")
	}
	if(statcode[16] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, median, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "median"
		out.label <- c(out.label,"Median")
	}
	if(statcode[17] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, data.frame(Var = i, Freq = paste(tapply(data[[i]], grp, mode.var, na.rm = TRUE)[[1]], collapse = ", ", sep = "")))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "mode"
		out.label <- c(out.label,"Mode")
	}
	if(statcode[20] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, quantile, probs = 0.25, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Q1"
		out.label <- c(out.label,"1st Quartile")
	}
	if(statcode[21] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, quantile, probs = 0.75, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "Q3"
		out.label <- c(out.label,"3rd Quartile")
	}
	if(statcode[9] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, max, na.rm = TRUE) - tapply(data[[i]], grp, min, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "range"
		out.label <- c(out.label,"Range")
	}
	if(statcode[10] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, quantile,probs = 0.75, na.rm = TRUE) - tapply(data[[i]], grp, quantile, probs = 0.25, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "iqr"
		out.label <- c(out.label,"Inter Quartile Range")
	}
	if(statcode[11] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, var, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "var"
		out.label <- c(out.label,"Variance")
	}
	if(statcode[12] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, sd, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "stdev"
		out.label <- c(out.label,"Standard Deviation")
	}
	if(statcode[13] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, se.mean, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "se.Mean"
		out.label <- c(out.label,"Std. Error of the Mean")
	}
	if(statcode[14] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, coef.var, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "cv"
		out.label <- c(out.label,"Coefficient of Variation")
	}
	if(statcode[5] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, css, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "CSS"
		out.label <- c(out.label,"Corrected Sum of Squares")
	}
	if(statcode[6] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, ucss, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "UCSS"
		out.label <- c(out.label,"Uncorrected Sum of Squares")
	}
	if(statcode[22] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, skew, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "skew"
		out.label <- c(out.label,"Skewness")
	}
	if(statcode[7] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, se.skewness, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE.Skew"
		out.label <- c(out.label,"Std. Error of Skewness")
	}
	if(statcode[23] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, kurtosis.var, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "kurtosis"
		out.label <- c(out.label,"Kurtosis")
	}
	if(statcode[8] == "1") {
		a <- NULL
		for (i in (1:ncol(data))) a <- rbind(a, as.data.frame.table(tapply(data[[i]], grp, se.kurtosis, na.rm = TRUE)))
		ifelse(is.null(summaryTable), summaryTable <- a, summaryTable <- cbind(summaryTable, a[ncol(a)]))
		colnames(summaryTable)[ncol(summaryTable)] <- "SE.Kurtosis"
		out.label <- c(out.label,"Std. Error of Kurtosis")
	}
	
	if (is.null(ncol(grp))) {
		summaryTable[,1] <- names(data)
	}
	else {
		variable <- c(rep(names(data), each = nrow(as.data.frame.table(table(grp)))))
		summaryTable <- data.frame(variable,summaryTable)
	}
	colnames(summaryTable)[1] <- "Variable"
	return(summaryTable)
}

