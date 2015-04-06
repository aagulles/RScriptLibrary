#----------------------------------------------------------------
# Generation Mean Analysis                
# uses Weighted Least Square Method Regression 
# starts with the full regression model                            
#                  
# Script Created by: Nellwyn L. Sales  05.07.2012                                     
#---------------------------------------------------------------- 

#----------------------------------------------------------------
# ARGUMENTS:
# data - a dataframe that contains summary statistics of the generations
# meanVar - string; name of the variable for the generation means
# weightsVar - string; name of the variable for the generation weights
# stdDevVar - string; name of the variable for the generation standard deviations
# numObsVar - string; name of the variable for the number of observations per generation
# generationVar - string; name of the variable for the generation 
# usersNotationVector - a vector of the user's notation for different generations
# generalNotationVector - a vector of the general notation for different generations
# alpha - a string, user's desired level of significance
#---------------------------------------------------------------

generationMeanSummaryStats.test<-function(data, meanVar, weightsVar, stdDevVar, numObsVar, generationVar, usersNotationVector, generalNotationVector, alpha) {
	
	#trim the strings
	meanVar <- trim.strings(meanVar)
	generationVar <- trim.strings(generationVar)
	usersNotationVector <- trim.strings(usersNotationVector)
	generalNotationVector <- trim.strings(generalNotationVector)
	if (!(is.null(weightsVar))) weightsVar <- trim.strings(weightsVar)
	if (!(is.null(stdDevVar))) stdDevVar <- trim.strings(stdDevVar)
	if (!(is.null(numObsVar))) numObsVar <- trim.strings(numObsVar)
	alpha <- trim.strings(alpha) 

	alphaValue<-as.double(alpha)
	
	#create new dataframe 'summaryStats' containing the mean and weights
	if (is.null(weightsVar)) {
		if (is.null(numObsVar)) {
			weights <-1/(data[,match(stdDevVar, names(data))]^2)
		} else {
			weights <-data[,match(numObsVar, names(data))]/(data[,match(stdDevVar, names(data))]^2) 
		}
		summaryStats<-data.frame(generation=data[,match(generationVar, names(data))], mean=data[,match(meanVar, names(data))],weights)

	} else {
		summaryStats<-data.frame(generation=data[,match(generationVar, names(data))], mean=data[,match(meanVar, names(data))],weights=data[,match(weightsVar, names(data))])
	}
	result <- list()
	result$summaryStats<-summaryStats

	#get rows in summaryStats that is included in the usersNotationVector
	#and merge the each row with the generalNotation, bSub, fSub, pSub and m column
	selectedGen<-NULL
	for (i in (1:length(usersNotationVector))) {
		#get subscripts of B,F and P
		b<-substr(generalNotationVector[i],2,2)
		f<-substr(generalNotationVector[i],4,4)
		p<-substr(generalNotationVector[i],6,6)
		
		#for each generation, assign corresponding values of a and d
		#check if P1
		if (b=="0" && f=="0" && p=="1") {
			a<-1
			d<-0
		}else {
			#check if P2
			if (b=="0" && f=="0" && p=="2") {
				a<- -1
				d<-0
			}else {
				#check if Fn
				if (b=="0" && p=="0") {
					a<-0
					d<-(0.5^(as.numeric(f)-1))
				}else {
					#check if backcrossed to P1
					if (f=="1" && p=="1") {
						a<-(1-(0.5^as.numeric(b)))
						d<-0.5^as.numeric(b)
					} else {
						#check if backcrossed to P2
						if (f=="1" && p=="2") {
							a<-(1-(0.5^as.numeric(b)))*(-1)
							d<-0.5^as.numeric(b)
						}
					}
				}
			}
		}
		aa<-a*a
		ad<-a*d
		dd<-d*d
		newRow <- cbind(subset(summaryStats, summaryStats$generation == usersNotationVector[i]), generalNotation=generalNotationVector[i], bSubscript=b, fSubscript=f, pSubscript=p, m = 1, a=a, d=d, aa=aa, ad=ad, dd=dd)
		selectedGen<-rbind(selectedGen,newRow)
	}
	selectedGen<-data.frame(selectedGen)
 	
	cat("\nGENERATION MEAN ANALYSIS")

	# Fitting Set A Full Model
	cat("\n\n-------------------------------------------------")
	cat("\nSET A FULL MODEL: mean ~ 0 + m + a + d + aa + ad")
	cat("\n-------------------------------------------------")
	
	# check if X matrix is full rank
	XdataFrame1<-data.frame(m=selectedGen$m, a=selectedGen$a, d=selectedGen$d, aa=selectedGen$aa, ad=selectedGen$ad)
	Xmatrix1<-as.matrix(XdataFrame1)
	library(Matrix)
	rankMatrix1<-rankMatrix(Xmatrix1)[1]
	result$Xmatrix1<-Xmatrix1
	result$rankXmatrix1<-rankMatrix1

	if (rankMatrix1 != min(nrow(Xmatrix1), ncol(Xmatrix1))) {
		cat("\n\n X matrix is not full rank. Full model cannot be fitted.")
		cat("\n Finding a reduced model with full rank X matrix...")
	}
	
	while (rankMatrix1 != min(nrow(Xmatrix1), ncol(Xmatrix1))) {
		corrMatrix<-cor(Xmatrix1)
		diag(corrMatrix)<-NA
		if (sum(is.na(corrMatrix)) < (nrow(corrMatrix)*ncol(corrMatrix))) {
			maxCorr<-max(abs(corrMatrix),na.rm=T)
			dependentColumnsList<-which(abs(corrMatrix)==maxCorr, arr.ind=TRUE)
			maxColumnIndex<-max(dependentColumnsList)
			Xmatrix1<-Xmatrix1[,-maxColumnIndex]
			rankMatrix1<-rankMatrix(Xmatrix1)[1]
		} else {
			break
		}
	}
	
	if (rankMatrix1 == min(nrow(Xmatrix1), ncol(Xmatrix1))) {
		myFormula1<-"mean ~ 0"
		for (i in 1:length(colnames(Xmatrix1))) {
			myFormula1<-paste(myFormula1, "+", colnames(Xmatrix1)[i]) 
		}
		runBackwardRegression(selectedGen, myFormula1, alphaValue, "A")
	} else {
		cat("\n\n ERROR: Cannot find a reduced model with full rank X matrix.\n")
	}

	# Fitting Set B Full Model	
	cat("\n\n-------------------------------------------------")
	cat("\nSET B FULL MODEL: mean ~ 0 + m + a + d + aa + dd")
	cat("\n-------------------------------------------------")

	# check if X matrix is full rank
	XdataFrame2<-data.frame(m=selectedGen$m, a=selectedGen$a, d=selectedGen$d, aa=selectedGen$aa, dd=selectedGen$dd)
	Xmatrix2<-as.matrix(XdataFrame2)
	library(Matrix)
	rankMatrix2<-rankMatrix(Xmatrix2)[1]
	result$Xmatrix2<-Xmatrix2
	result$rankXmatrix2<-rankMatrix2

	if (rankMatrix2 != min(nrow(Xmatrix2), ncol(Xmatrix2))) {
		cat("\n\n X matrix is not full rank. Full model cannot be fitted.")
		cat("\n Finding a reduced model with full rank X matrix...")
	}

	
	while (rankMatrix2 != min(nrow(Xmatrix2), ncol(Xmatrix2))) {
		corrMatrix<-cor(Xmatrix2)
		diag(corrMatrix)<-NA
		if (sum(is.na(corrMatrix)) < (nrow(corrMatrix)*ncol(corrMatrix))) {
			maxCorr<-max(abs(corrMatrix),na.rm=T)
			dependentColumnsList<-which(abs(corrMatrix)==maxCorr, arr.ind=TRUE)
			maxColumnIndex<-max(dependentColumnsList)
			Xmatrix2<-Xmatrix2[,-maxColumnIndex]
			rankMatrix2<-rankMatrix(Xmatrix2)[1]
		} else {
			break
		}
	}
	if (rankMatrix2 == min(nrow(Xmatrix2), ncol(Xmatrix2))) {
		myFormula2<-"mean ~ 0"
		for (i in 1:length(colnames(Xmatrix2))) {
			myFormula2<-paste(myFormula2, "+", colnames(Xmatrix2)[i]) 
		}	
		runBackwardRegression(selectedGen, myFormula2, alphaValue, "B")	
	} else {
		cat("\n\n ERROR: Cannot find a reduced model with full rank X matrix.\n")
	}

	return(list(output = result))
}
