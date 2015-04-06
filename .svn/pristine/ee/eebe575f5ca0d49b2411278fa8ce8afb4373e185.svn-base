#-----------------------------------------------------------------------------------------------------------------------#
#     PRINCIPAL COMPONENT ANALYSIS(PCA) FUNCTION 09.25.2012 MODIFIED 03.18.13                                           #
#     PCA(data, var, idVar, descriptiveStat, matx, transform, saveScore, outputPath, scree, biplot, pcaplot, useIdVar)	#
#-----------------------------------------------------------------------------------------------------------------------#

PCA <-function(data, var, idVar = NULL, descriptiveStat = TRUE, matx = c("corr","cov"), transform = c("zerocenter", "unitvar", "none"), 
               saveScore = TRUE, outputPath = NULL,scree = TRUE, biplot = TRUE, pcaplot = TRUE, useIdVar = FALSE, pChars, pSizes, pCol, 
               showLeg = FALSE, legTitle, legPos = "bottomright", legNcol = 1, axesNum = 3) 

UseMethod("PCA")

PCA.default <-function(data, var, idVar = NULL, descriptiveStat = TRUE, matx = c("corr","cov"), transform = c("zerocenter", "unitvar", "none"), 
                 saveScore = TRUE, outputPath = NULL, scree = TRUE, biplot = TRUE, pcaplot = TRUE, useIdVar = FALSE, pChars, pSizes, pCol, 
                       showLeg = FALSE, legTitle, legPos = "bottomright", legNcol = 1, axesNum = 3) 
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
  	numObsNM = nrow((na.omit(tempData[,c(var,idVar)])))
  	numObs = nrow((tempData[,c(var,idVar)]))
  	
  	if (numObs == numObsNM) { cat("Number of Observations:", numObs,"\n")
  	} else {
  	  cat("Number of Observations:", numObs,"\n")
  	  cat("Number of Observations Used:", numObsNM,"\n\n")
  	}
    
  	tempData <- na.omit(tempData[,c(var,idVar)])
    
	if (!is.data.frame(tempData)) { stop("The object should be of type data frame.") }
	if (!is.character(var)) 	{ stop(paste("The object 'var' should be a character vector.", sep = "")) }
	if (any(is.na(match(var, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }


  	if (!is.null(idVar)) 		{ 
  	  if (any(is.na(match(idVar, names(tempData))))) { stop("At least one item in the character vector 'var' does not match any variable name in the dataset.") }
  	  tempGrp <- tempData[idVar]
  	  for (i in (1:ncol(tempGrp))) { tempGrp[,i] <- factor(tempGrp[,i]) }
  	} else { tempGrp <- rep(1, each = nrow(tempData)) }
  	tempVar <- tempData[var]
    
  	transform <- match.arg(transform)
  	options(width = 6000, digits = 4)

	#### ---- Output of the analysis ---- ####  

	if (matx == "corr"){
	pc <- prcomp(tempData[,var], scale = TRUE)
	}
	if (matx == "cov"){
		
		if(transform == "zerocenter"){pc <- prcomp(tempData[,var], center = TRUE, scale = FALSE)}
		if(transform == "unitvar"){pc <- prcomp(tempData[,var], center = FALSE, scale = TRUE)}
		if(transform == "none"){pc <- prcomp(tempData[,var], center = FALSE, scale = FALSE)}		
	}

	#### ---- Function for displaying Descriptive Statistics ---- ####
  if(descriptiveStat){
  	DescriptiveStatistics(data = tempData, var = var, grp = NULL, statistics = c("nnmiss", "mean", "sd", "min", "max") )
  	cat("\nCORRELATION MATRIX\n")
  	cat("\n")
    cor(tempVar, method = "pearson")
  	outCor <- round(cor(tempVar, method = "pearson"),4)
    print(outCor)
  	cat("\nCOVARIANCE MATRIX\n")
  	cat("\n")
    cov(tempVar)
    outCov<- round(cov(tempVar),4)
  	print(outCov)
  }
    
	  pctable<- rbind(summary(pc)$importance, pc$sdev^2) 
  	rownames(pctable)[nrow(pctable)] <- "EigenValues"

  	cat("\nPRINCIPAL COMPONENT ANALYSIS\n")
		pctable<-data.frame(pctable)
		pcname<-data.frame(dimnames(pctable)[[1]])
		colnames(pcname) <- "Statistics"
		pcstat<-cbind(pcname, pctable)
		printDataFrame(pcstat, digits=4)
  	cat("\n")

  	cat("\nEIGENVECTORS\n")
  	vect<-data.frame(pc$rotation)	
		pcvname<-data.frame(dimnames(vect)[[1]])
		colnames(pcvname) <- "Variables"
		pcvect<-cbind(pcvname, vect)
		printDataFrame(pcvect, digits=4)
  	cat("\n")

# 	  cat("\nPRINCIPAL COMPONENTS SCORES\n")
# 	  score<- data.frame(pc$x)
  	score<-data.frame(pc$x[1:nrow(tempData), 1:axesNum])
    pcscore <-score
	  if(!is.null(idVar)){
      pcscore <- cbind(data.frame(tempGrp), score)
#       printDataFrame(pcscore, digits= 4)
	  } ##else 
# 	  printDataFrame(pcscore, digits= 4)
 
    
	  if (!is.null(outputPath)) {  
	    if(scree){
	      png(filename = paste(outputPath,"ScreePlot.png",sep=""))  
	      screeplot(pc, type="l", ylab="EigenValue", main="Scree Plot")
	      dev.off()
	    } 
      if(biplot){
        for(i in 1:(axesNum-1)){
          a=i+1
          for(j in a:axesNum){
            png(filename = paste(outputPath,"Biplot", i, "and", j, ".png",sep=""))
            biplot(pc, choices = c(i,j), cex=0.8, expand = 1, main = "Biplot")
            dev.off()
          }
        }
	    }
      if (pcaplot){   
        if(!is.null(idVar)){
	            idLevels=levels(factor(tempData[, idVar]))
	            pCharNames <- make.unique(c(names(tempData),"pCharCode"), sep = "")
	            pCharName <- pCharNames[length(pCharNames)]
	            pColName <- make.unique(c(pCharNames,"pColCode"), sep = "")
	            pColName <- pColName[length(pColName)]
	           
              tempData[,pCharName]<-NA
	            tempData[,pColName]<-NA
	            if(!useIdVar) {
                grpCode<- data.frame(idLevels, pChars,pCol)
	              tempData[,pCharName] <- grpCode$pChars[match(tempData[,idVar], grpCode[,"idLevels"])]
                tempData[,pColName]<- grpCode$pCol[match(tempData[,idVar], grpCode[,"idLevels"])]
	            } 
              
	            x<- pc$x
	            y<- x[1:nrow(x), 1:axesNum]
	            if (useIdVar) { 
	              plotType = "n" 
	              pChars = idLevels
	              pCharData = rownames(x)
	            } else plotType = "p"
	                   pCharData = tempData[,pCharName]
	            for(i in 1:(ncol(y)-1)){
	              a=i+1
	              for(j in a:(ncol(y))){
	                png(filename = paste(outputPath,"PC", i, "and", j, ".png",sep=""))
	                plot(y[,i], y[,j], xlab = paste("PC", i), ylab = paste("PC", j),  main = "PCA Plot",
	                type = plotType, pch = pCharData, cex = pSizes, col = as.vector(tempData[,pColName]))
	                if(useIdVar)text(y[,i], y[,j], labels =  tempData[, idVar], cex = pSizes, col = as.vector(tempData[,pColName]))
	                if (showLeg) legend(legPos, title = legTitle, legend = idLevels, pch = pChars, cex = pSizes, col = pCol, ncol = legNcol)	                
	                dev.off()
	              }
	             }
	          }
          else{
            x<- pc$x
            y<- x[1:nrow(x), 1:axesNum]
            for(i in 1:(ncol(y)-1)){
              a=i+1
              for(j in a:(ncol(y))){
                png(filename = paste(outputPath,"PC", i, "and", j, ".png",sep=""))
                plot(y[,i], y[,j], xlab = paste("PC", i), ylab = paste("PC", j),  main = "PCA Plot")
                dev.off()
              }
            }
          }	      
	    }
      
		  if (saveScore){
		  	write.csv(data.frame(Rownames= rownames(pcscore), pcscore), row.names = FALSE, file=paste(outputPath, "PCScores.csv", sep = ""))
		  }
	  }       
  	return(list(PC = pc, SummaryTable = pctable))
  
}### end stmt (PCA Function)