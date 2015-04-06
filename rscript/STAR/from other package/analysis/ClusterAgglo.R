#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#     CLUSTER ANALYSIS FUNCTION (AGGLOMERATIVE) 09.25.2012                      								     				     	           #
#     ClusterAgglo <- function(data, var, stand= TRUE, distance , clusmethod, distMatrix, clusterMem, descriptiveStat, dendrogram, clusterBox, clusterNum, outputPath) #
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------#

ClusterAgglo <- function(data, var, stand= TRUE, distance = c("Euclidean", "Maximum", "Manhattan", "Minkowski", "Canberra"), 
			clusmethod = c("Single", "Complete", "Average", "Ward", "Centroid"),distMatrix = FALSE, clusterMem = TRUE, 
			descriptiveStat = TRUE, dendrogram = TRUE, clusterBox = TRUE, clusterNum=2, outputPath = NULL)
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
  
	distance <- match.arg(distance)
	clusmethod <- match.arg(clusmethod)
	options(width = 5000, digits = 6)

	if (distance == "Euclidean"){
    			if (stand){
      			a <- data.Normalization(tempData[,var], type = "n1")
      			adist <- dist(a, method = "euclidean")
    			}else
      			adist <- dist(tempData[,var], method = "euclidean")	
  		}
	if (distance == "Maximum"){
    			if (stand){
      			a <- data.Normalization(tempData[,var], type = "n1")
      			adist <- dist(a, method = "maximum")
    			}else
      			adist <- dist(tempData[,var], method = "maximum")	
  		}
	if (distance == "Manhattan"){
    			if (stand){
      			a <- data.Normalization(tempData[,var], type = "n1")
      			adist <- dist(a, method = "manhattan")
    			}else
      			adist <- dist(tempData[,var], method = "manhattan")	
  		}
	if (distance == "Minkowski"){
    			if (stand){
      			a <- data.Normalization(tempData[,var], type = "n1")
      			adist <- dist(a, method = "minkowski")
    			}else
      			adist <- dist(tempData[,var], method = "minkowski")	
  		}
	if (distance == "Canberra"){
    			if (stand){
      			a <- data.Normalization(tempData[,var], type = "n1")
      			adist <- dist(a, method = "canberra")
    			}else
      			adist <- dist(tempData[,var], method = "canberra")	
  		}

	cat("\nAGGLOMERATIVE CLUSTER ANALYSIS\n")

	if (distMatrix){
	cat("\nDISTANCE MATRIX\n")
	print(adist)
  	}

	if (clusmethod == "Complete"){	
    		amethod <- hclust(adist, method = "complete")}
  	if (clusmethod == "Single"){	
    		amethod <- hclust(adist, method = "single")}
  	if (clusmethod == "Average"){	
    		amethod <- hclust(adist, method = "average")}
  	if (clusmethod == "Ward"){	
  		amethod <- hclust(adist, method = "ward")}
  	if (clusmethod == "Centroid"){	
   	 	amethod <- hclust(adist, method = "centroid")}
  
	Membership <- cutree(amethod, k = as.numeric(clusterNum))
	memData <- data.frame(Membership)
	memSummary <- cbind(tempData, memData)
	memberList <-list()

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
	    	cat("\nNumber of members in each cluster\n\n")
    		cat("Cluster number\t")
	    	for (i in (1:attributes(table(Membership))$dim)) {
    		  	cat("\t",attributes(table(Membership))$dimnames$Membership[i])
	    	}
    		cat("\n")
	    	cat("Cluster Size\t")
    		for (i in (1:attributes(table(Membership))$dim)) {
      		cat("\t",table(Membership)[i])

	    	}
	}
	if (!is.null(outputPath)) {
		if (dendrogram){
			png(paste(outputPath, "AggloGraph.png", sep = ""))
			dendro <- as.dendrogram(amethod)
			plot(amethod)
			plot(dendro, horiz = FALSE, center = TRUE,
				nodePar = list(lab.cex = 0.6, lab.col = "forest green", pch = NA),
				main = "Dendrogram using Agglomerative Clustering Method")
  	
			if(clusterBox){
				rect.hclust(tree=amethod, k= as.numeric(clusterNum), border=c("red", "blue", "green", "purple"))
			}
			dev.off()
		}
	}

	if (descriptiveStat){
		cat("\n")
		cat("\n")
		all <- cbind(tempData[,var], Membership)
		DescriptiveStatistics(data = all, var = var, grp = "Membership", statistics = c("mean", "sd", "min", "max") )
	}

  	cat("\nCOPHENETIC CORRELATION\n")
  	Copcorr <- cophenetic(amethod)
	cop <- cor(adist, Copcorr)  
  	print(cop)

	return(list(ClusterMethod = amethod, Membership = memberList, MemberSummary = memSummary, DistMatrix = adist, CopheneticCorrelation = cop))

}#### end statement ClusterAgglo####
