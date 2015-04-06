# --------------------------------------------------
# BALANCED INCOMPLETE BLOCK DESIGN
# File Created by: Alaine A. Gulles 04.13.2011
# File Modified by: Alaine A. Gulles 04.18.2011
# Note: included in RCROPSTAT_UTILITIES_STAT_NEW
# --------------------------------------------------

BIBDTest <- function(data, respvar, trmt, block, method = NULL, descriptive = FALSE, alpha = 0.05) UseMethod("BIBDTest")

BIBDTest.default <- function(data, respvar, trmt, block, method = NULL, descriptive = FALSE, alpha = 0.05) {

	if (is.character(data)) { data <- eval(parse(text = data)) } 
	data[,block] <- factor(data[,block]) 
	data[,trmt] <- factor(data[,trmt]) 

	a <- nlevels(data[,trmt]) # -- trmt levels
	b <- nlevels(data[,block]) # -- blk levels
	r <- unique(table(data[,trmt]))
	k <- unique(table(data[,block]))
	lambda <- (r * (k - 1))/(a - 1)

	procedure <- c("lsd", "tukey", "snk", "duncan")
	method <- procedure[na.omit(match(method, procedure))]
	if (length(method) == 0) method <- NULL

	LSDTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff) {
		Tprob <- qt(1 - alpha/2, dfError)
		thestring <- paste(paste(capture.output(print(output1 <- order.group(levels(data[,trmt]),
		      	      TrmtMean.adj, rep(1, a), MSError, Tprob, StdErr.adjtrtmean,
					k /(lambda * a), StdErr.diff))), collapse = "\n"), collapse = "\n\n")
		output1[,4] <- r
		output1[,5] <- StdErr.diff
		colnames(output1)[1] <- trmt
		output1 <- cbind(output1[,1:2], std.err = output1[,5], N = output1[,4], grouping = output1[,3])
		return(list(method = "Least Significant Difference (LSD) Test", cval = Tprob * StdErr.diff, stderr = StdErr.diff, pw = output1))				
	}

	TukeyTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff) {
		Tprob <- qtukey(1 - alpha, a, dfError)
		thestring <- paste(paste(capture.output(print(output1 <- order.group(levels(data[,trmt]),
				      TrmtMean.adj, rep(1, a), MSError, Tprob, StdErr.adjtrtmean,
					(k /(lambda * a))/2, StdErr.diff))), collapse = "\n"), collapse = "\n\n")
		output1[,4] <- r
		output1[,5] <- StdErr.diff/sqrt(2)
		colnames(output1)[1] <- trmt
		output1 <- cbind(output1[,1:2], std.err = output1[,5], N = output1[,4], grouping = output1[,3])
		return(list(method = "Tukey's Honestly Significant Difference (HSD) Test", 
				cval = Tprob * StdErr.diff/sqrt(2), stderr = StdErr.diff/sqrt(2), pw = output1))
	}

	DuncanTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff){
		Tprob <- qtukey((1 - alpha)^(1:(a-1)), 2:a, dfError)
		duncan <- Tprob * StdErr.diff/sqrt(2)
		critical.range <- rbind(tabular.value = Tprob, critical.value = duncan)
		colnames(critical.range) <- 2:a
		thestring <- paste(paste(capture.output(print(output1 <- order.group(levels(data[,trmt]),
			            TrmtMean.adj, rep(1, a), MSError, Tprob, StdErr.adjtrtmean,
					k /(lambda * a), 2, dfError, alpha, StdErr.diff/sqrt(2)))), collapse = "\n"), collapse = "\n\n")
		output1[,4] <- r
		output1[,5] <- StdErr.diff/sqrt(2)
		colnames(output1)[1] <- trmt
		output1 <- cbind(output1[,1:2], std.err = output1[,5], N = output1[,4], grouping = output1[,3])
		return(list(method = "Duncan's Multiple Range test", crange = critical.range, stderr = StdErr.diff/sqrt(2), pw = output1))
	}

	SNKTest <- function(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff){
		Tprob <- qtukey((1 - alpha), 2:a, dfError)
		SNK <- Tprob * StdErr.diff/sqrt(2)
		critical.range <- rbind(tabular.value = Tprob, critical.value = SNK)
		colnames(critical.range) <- 2:a
		thestring <- paste(paste(capture.output(print(output1 <- order.group(levels(data[,trmt]),
			           TrmtMean.adj, rep(1, a), MSError, Tprob, StdErr.adjtrtmean,
					k /(lambda * a), 2, dfError, alpha, StdErr.diff/sqrt(2)))), collapse = "\n"), collapse = "\n\n")
		output1[,4] <- r
		output1[,5] <- StdErr.diff/sqrt(2)
		colnames(output1)[1] <- trmt
		output1 <- cbind(output1[,1:2], std.err = output1[,5], N = output1[,4], grouping = output1[,3])
		return(list(method = "Student Newman Keuls (SNK) Test", crange = critical.range, stderr = StdErr.diff/sqrt(2), pw = output1))
	}

	anova.table <- list()
	trmtStat <- NULL
	blkStat <- NULL
	tempStat <- NULL
	pw <- list()

	for (i in (1:length(respvar))) {
		result1 <- lm(formula(paste(respvar[i], "~", block, "+", trmt)), data)
		temp.anova1 <- anova(result1)
		rownames(temp.anova1)[1] <- paste(rownames(temp.anova1)[1], "(unadj)", sep = "")
		rownames(temp.anova1)[2] <- paste(rownames(temp.anova1)[2], "(adj)", sep = "")
		rownames(temp.anova1)[nrow(temp.anova1)] <- "Error"
		temp.anova1[nrow(temp.anova1)+1, 1:2] <- c(sum(temp.anova1[,1]), sum(temp.anova1[,2]))
		rownames(temp.anova1)[nrow(temp.anova1)] <- "Total"
		temp.anova1[1,3:5] <- NA

		result2 <- lm(formula(paste(respvar[i], "~", trmt, "+", block)), data)
		temp.anova2 <- anova(result2)	
		rownames(temp.anova2)[1] <- paste(rownames(temp.anova2)[1], "(unadj)", sep = "")
		rownames(temp.anova2)[2] <- paste(rownames(temp.anova2)[2], "(adj)", sep = "")
		temp.anova2[1,3:5] <- NA
		anova.table[[i]] <- rbind(temp.anova1[2,], temp.anova2[1,], temp.anova1[1,], temp.anova2[2,], temp.anova1[3,], temp.anova1[4,])
		dfError <- temp.anova1[nrow(temp.anova1)-1,1]
		MSError <- temp.anova1[nrow(temp.anova1)-1,3]
		options(width = 5000)
		cat("Balanced Incomplete Block Design\n")
		ClassInformation(data[,c(block, trmt, respvar[i])], respvar = respvar[i])
          
          if (descriptive) {
               cat("\n\n")
               DescriptiveStatistics(data, var = respvar[i], statistics = c("n", "mean", "sd", "min", "max"))
          }
          
     	cat("\n\nAnalysis of Variance Table\n")
		cat("Response Variable:", respvar[i],"\n")		
		printAOVTable(anova.table[[i]])

		# --- COMPUTE SOME STATISTICS --- #
		blk.sum <- tapply(data[,respvar[i]], data[,block], sum)
		trt.sum <- tapply(data[,respvar[i]], data[,trmt], sum)
		trt.mean <- tapply(data[,respvar[i]], data[,trmt], mean, na.rm = TRUE)
		cmb.mean <- tapply(data[,respvar[i]], list(data[,trmt],data[,block]), mean, na.rm = TRUE)
		num <- (!is.na(cmb.mean)) %*% matrix(blk.sum)
		trmtTotal.adj <- matrix(trt.sum) - (as.numeric(num)/k)
		SSTrmt.adj <- (colSums(trmtTotal.adj * trmtTotal.adj) * k)/(lambda * a)
		MSTrmt.adj <- SSTrmt.adj/(a - 1)
		StdErr.diff <- sqrt((2 * k * MSError)/(lambda * a))
		StdErr.adjtrteff <- sqrt((k * MSError)/(lambda * a))
		StdErr.adjtrtmean <- sqrt((MSError * (1 + ((k * r * (a - 1))/(lambda * a))))/(a * r))
		TrmtEffect.adj <- trmtTotal.adj * (k /(lambda * a))    # --- TrmtEffect.adj = intrablock estimate
		interblk.est <- (num - (k * r * mean(data[,respvar[i]])))/(r - lambda)
		TrmtMean.adj <- mean(data[,respvar[i]]) + TrmtEffect.adj

		cat("\nTable of Means\n", sep = "")
		tempTable <- data.frame(RespVar = respvar[i], TrmtMean = trt.mean, AdjTrmtMean = TrmtMean.adj, AdjTrmtEffect = TrmtEffect.adj)
		tempTable <- data.frame(rownames(tempTable), tempTable)
		colnames(tempTable)[1] <- trmt
		printDataFrame(tempTable)
		trmtStat <- rbind(trmtStat, tempTable)
		cat("\nStandard Error of the ...\n", sep = "")
		tempTable <- data.frame(Diff = StdErr.diff, AdjTrmtEffect = StdErr.adjtrteff, AdjTrmtMean = StdErr.adjtrtmean)
		printDataFrame(tempTable)
		tempStat <- rbind(tempStat, data.frame(RespVar = respvar[i], tempTable))
		cat("\nInterblock Estimate:\n", sep = "")
		tempTable <- data.frame(RespVar = respvar[i], t(interblk.est))
		printDataFrame(tempTable)
		blkStat <- rbind(blkStat, tempTable)
          
		if (anova.table[[i]][match(paste(trmt, "(adj)", sep = ""), rownames(anova.table[[i]])),5] < alpha) {
		     pw.result <- list()
		     cat("\nPairwise Mean Comparison for Variable", respvar[i],"\n\n")
		     if (!is.null(method)) {
		          for (j in (1:length(method))) {
		               if (method[j] == "lsd") {
		                    pw.result[[j]] <- LSDTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff) 
		                    StatLabel <- "LSD"
		               }
		               if (method[j] == "tukey") { 
		                    pw.result[[j]] <- TukeyTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff) 
		                    StatLabel <- "HSD"
		               }
		               if (method[j] == "snk") { 
		                    pw.result[[j]] <- SNKTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff)
		                    StatLabel <- "SNK"
		               }
		               if (method[j] == "duncan") { 
		                    pw.result[[j]] <- DuncanTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff) 
		                    StatLabel <- "duncan"
		               }
		               pw.result[[j]]$pw <- pw.result[[j]]$pw[order(pw.result[[j]]$pw[,trmt]),]
		               rownames(pw.result[[j]]$pw) <- 1:nrow(pw.result[[j]]$pw)
		               
		               cat(pw.result[[j]]$method, "\n\n")
		               cat(formatC("Alpha", format = "s", width = 24, flag = "-"), formatC(alpha, format = "f", digits = 2, width = 10, flag = "#"), "\n", sep = "")
		               cat(formatC("Error Degrees of Freedom", format = "s", width = 24, flag = "-"), formatC(dfError, format = "d", width = 10, flag = "#"), "\n", sep = "")
		               cat(formatC("Error Mean Square", format = "s", width = 24, flag = "-"), formatC(MSError, format = "f", digits = 4, width = 10, flag = "#"), "\n", sep = "")
		               if (method[j] == "lsd" || method[j] == "tukey") {
		                    cat(formatC("Critical Value", format = "s", width = 24, flag = "-"), formatC(pw.result[[j]]$cval/pw.result[[j]]$stderr, format = "f", digits = 4, width = 10, flag = "#"), "\n", sep = "")
		                    cat(formatC(paste(StatLabel, "Test Statistics"), format = "s", width = 24, flag = "-"), formatC(pw.result[[j]]$cval, format = "f", digits = 4, width = 10, flag = "#"), "\n\n", sep = "")
		               } else {
		                    cat("\n")
		                    cat(formatC("Number of Means", format = "s", width = 24, flag = "-"), formatC(colnames(pw.result[[j]]$crange), format = "d", width = 10, flag = "#"), "\n") 
		                    cat(formatC("Critical Value", format = "s", width = 24, flag = "-"), formatC(pw.result[[j]]$crange[1,], format = "f", digits = 4, width = 10, flag = "#"), "\n")
		                    cat(formatC(paste(StatLabel,"Test Statistics"), format = "s", width = 24, flag = "-"), formatC(pw.result[[j]]$crange[2,], format = "f", digits = 4, width = 10, flag = "#"), "\n\n")
		               }
		               printDataFrame(pw.result[[j]]$pw, border = TRUE)
		               cat("* Means with the same letter are not significantly different. \n\n")
		          }
		          pw[[i]] <- pw.result
		     } else {
		          if (nlevels(data[,trmt]) <= 5) {
		               pw.result[[1]] <- LSDTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff)
		               StatLabel <- "LSD"
		          } else {
		               pw.result[[1]] <- TukeyTest(dfError, MSError, TrmtMean.adj, StdErr.adjtrtmean, StdErr.diff) 
		               StatLabel <- "HSD"
		          }
		          pw.result[[1]]$pw <- pw.result[[1]]$pw[order(pw.result[[1]]$pw[,trmt]),]
		          rownames(pw.result[[1]]$pw) <- 1:nrow(pw.result[[1]]$pw)
		          
		          cat(pw.result[[1]]$method, "\n\n")
		          cat(formatC("Alpha", format = "s", width = 24, flag = "-"), formatC(alpha, format = "f", digits = 2, width = 10, flag = "#"), "\n", sep = "")
		          cat(formatC("Error Degrees of Freedom", format = "s", width = 24, flag = "-"), formatC(dfError, format = "d", width = 10, flag = "#"), "\n", sep = "")
		          cat(formatC("Error Mean Square", format = "s", width = 24, flag = "-"), formatC(MSError, format = "f", digits = 4, width = 10, flag = "#"), "\n", sep = "")
		          cat(formatC("Critical Value", format = "s", width = 24, flag = "-"), formatC(pw.result[[1]]$cval/pw.result[[1]]$stderr, format = "f", digits = 4, width = 10, flag = "#"), "\n", sep = "")
		          cat(formatC(paste(StatLabel, "Test Statistics"), format = "s", width = 24, flag = "-"), formatC(pw.result[[1]]$cval, format = "f", digits = 4, width = 10, flag = "#"), "\n\n", sep = "")
		          printDataFrame(pw.result[[1]]$pw, border = TRUE)
		          cat("* Means with the same letter are not significantly different. \n\n")
		          pw[[i]] <- pw.result[[1]]
		     }
		} else { pw[[i]] <- NULL }
	}
	rownames(tempStat) <- 1:nrow(tempStat)
	rownames(trmtStat) <- 1:nrow(trmtStat)
	rownames(blkStat) <- 1:nrow(blkStat)

	# --- OUTPUT --- #
	if (length(pw) >= 1) {
		return(invisible(list(anovaTable = anova.table,
				trmtStatistic = trmtStat,
				stdErr = tempStat,
				InterBlkEst = blkStat,
				alpha = alpha,
				pairwise = pw)))
	} else {
		return(invisible(list(anovaTable = anova.table,
				trmtStatistic = trmtStat,
				stdErr = tempStat,
				InterBlkEst = blkStat)))
	}

}
