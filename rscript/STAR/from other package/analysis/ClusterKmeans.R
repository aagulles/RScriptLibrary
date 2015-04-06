#--------------------------------------------------------------------------------------------------------#
#     CLUSTER ANALYSIS FUNCTION (KMEANS)09.25.2012      						               #
#     ClusterKmeans <- function(data, var, clusterMem, descriptiveStat, kgraph, clusterNum, outputPath)  #
#--------------------------------------------------------------------------------------------------------#

ClusterKmeans <- function(data, var, clusterMem = TRUE, descriptiveStat= TRUE,kgraph = TRUE, clusterNum=2, outputPath = NULL)
{
  	if (is.character(data)) { 
		nameData <- data
		if (!exists(nameData)) { stop(paste("The object '", nameData,"' does not exists.", sep = "")) }
		tempData <- eval(parse(text = data)) 
  	}else{
  	if (is.data.frame(data)) { 
		nameData <- paste(deparse(substitute(data)))	
		tempData <- data	
		} else { stop ("The argument should either be a data frame or a character string indicating the name of the data frame.") }
  	}  
	if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
	if (!is.character(var)) 	{ stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
  
	options(width = 5000, digits = 6)

	km <- kmeans(tempData[,var], centers = as.numeric(clusterNum))
	cat("K-MEANS CLUSTER ANALYSIS\n")
	Membership <- km$cluster
	memData <- data.frame(Membership)
	memSummary <- cbind(tempData, memData)
	memberList <- list()


	if (clusterMem){
	cat("\nCLUSTER MEMBERSHIP\n")
		for (i in (1:as.numeric(clusterNum))) {
			cat("\nMember of Cluster =",i,"\n")
			temp <- rownames(subset(memData, Membership == i))
			memberList[[i]] <- rownames(subset(memData, Membership == i))
	    		names(memberList)[i] <- paste("Cluster Number:", i)
			index <- 1
		for (j in (1:ceiling(length(temp)/15))) {
			if(index+14 > length(temp)) { cat(temp[index:length(temp)], "\n")
			} else { cat(temp[index:(index+14)], "\n") }
			index <- index + 15
			}
		}
		cat("\n")
		cat("MEANS FOR EACH CLUSTER\n")
		print(km$centers)
		cat("\n")
		cat("WITHIN SUM OF SQUARES FOR EACH CLUSTER\n")
		print(km$withinss)
		cat("\n")
		cat("SIZE FOR EACH CLUSTER\n")
		print(km$size)   	
		cat("\n")
	}
	
	if (descriptiveStat){
	all <- cbind(tempData[,var], Membership)
	DescriptiveStatistics(data = all, var = var, grp = "Membership", statistics = c("sd", "min", "max") )
	}

	if (!is.null(outputPath)) {
		if (kgraph){
			png(paste(outputPath, "Kmeans.png", sep = ""))
	  		plot(tempData[,var], col=Membership, pch=Membership)
			dev.off()
		}
	}
	return(list(ClusterMethod = km, Membership = memberList, MemberSummary = memSummary))

}#### end statement ClusterKmeans####