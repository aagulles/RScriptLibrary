# ------------------------------------------------------------------
# ssaTestPrep
# Description: Performs Single Site Analysis for P-Rep Design
# Scripts Created by: Violeta I Bartolome for IRRI
# Scripts Modified by: Alaine A Gulles
# Functions Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------
# Arguments: (ssaTestPrep)
# data
# respvar
# geno
# row
# column
# env = NULL
# is.random = FALSE
# excludeCheck = FALSE
# checkList = NULL 
# moransTest = FALSE
# spatialStruc = c("none", "CompSymm", "Gaus", "Exp", "Spher")
# descriptive = FALSE
# varCorr = FALSE
# heatmap = FALSE
# diagplot = FALSE
# outputPath = NULL
# ------------------------------------------------------------------

ssaTestPrep <- function(data, respvar, geno, row, column, env = NULL, is.random = FALSE, excludeCheck = FALSE, checkList = NULL, 
                        moransTest = FALSE, spatialStruc = c("none", "CompSymm", "Gaus", "Exp", "Spher"), descriptive = FALSE, 
                        varCorr = FALSE, heatmap = FALSE, diagplot = FALSE, histogram = FALSE, boxplot = FALSE,
                        outputPath = NULL) UseMethod("ssaTestPrep")     

ssaTestPrep.default <- function(data, respvar, geno, row, column, env = NULL, is.random = FALSE, excludeCheck = FALSE, checkList = NULL, 
                        moransTest = FALSE, spatialStruc = c("none", "CompSymm", "Gaus", "Exp", "Spher"), descriptive = FALSE, 
                        varCorr = FALSE, heatmap = FALSE, diagplot = FALSE, histogram = FALSE, boxplot = FALSE, 
                        outputPath = NULL) {
     # load the library needed for the analysis
     library(nlme)
     #library(ape)
     #library(gdata)
     
     # define variables and initializes
     result <- list()

     
     # if env column is not specified, create EnvLevel column
     if (is.null(env)) {
          env <- make.unique(c(names(data), "Env"))[length(make.unique(c(names(data), "Env")))]
          data <- cbind(data, 1)
          colnames(data)[ncol(data)] <- env
          addEnv <- TRUE
     } else { addEnv <- FALSE }
     
     # set environment to factor
     data[,env] <- factor(data[,env])     

     genoFactorType <- ifelse(is.random, "Random", "Fixed")
     cat(paste(rep("=", 25), collapse = ""), "\n")
     cat("GENOTYPE AS:", genoFactorType,"\n")
     cat(paste(rep("=", 25), collapse = ""), "\n\n")
     
     if (is.null(outputPath)) { outputPath <- paste(getwd(), "/", sep = "") }
     
     for (i in (1:length(respvar))) {
          width1 <- 25 + nchar(respvar[i])
          cat(paste(rep("-", width1), collapse = ""), "\n")
          cat("RESPONSE VARIABLE:", respvar[i],"\n")
          cat(paste(rep("-", width1), collapse = ""), "\n\n")
          if (addEnv) { capture.output(descStat <- DescriptiveStatistics(data, var = respvar[i], grp = NULL, statistics = c("nnmiss", "mean","sd"))) 
          } else { capture.output(descStat <- DescriptiveStatistics(data, var = respvar[i], grp = env, statistics = c("nnmiss", "mean","sd"))) }
          cat("DESCRIPTIVE STATISTICS:\n\n", sep = "")
          printDataFrame(descStat, border = FALSE)
          cat("\n\n")
          
          result[[i]] <- list()
          result[[i]]$design <- "p-rep design"
          result[[i]]$respvar <- respvar[i]
          
          for (j in (1:nlevels(data[,match(env, names(data))]))) {
               
               result[[i]]$site[[j]] <- list()
               result[[i]]$site[[j]]$env <- env
               result[[i]]$site[[j]]$envLevel <- levels(data[,match(env, names(data))])[j]
               result[[i]]$site[[j]]$envLabel <- addEnv
               
               if (!addEnv) {
                    width2 <- 23 + nchar(levels(data[,match(env, names(data))])[j]) + nchar(env)
                    cat(paste(rep("-", width2), collapse = ""), "\n")
                    cat("ANALYSIS FOR:", env, "=", levels(data[,match(env, names(data))])[j],"\n")
                    cat(paste(rep("-", width2), collapse = ""), "\n\n")
               }
               
               # --- create a temp data with one environment level only
               tempData <- subset(data, data[,match(env, names(data))] == levels(data[,match(env, names(data))])[j])
               if (all(is.na(tempData[,match(respvar[i], names(tempData))]))) { next }
               SSAResult <- try(ssaTestPrepMain(data = tempData, respvar = respvar[i], geno, row, column, is.random, excludeCheck, 
                                                checkList, moransTest, spatialStruc, varCorr), silent = TRUE) 
               
               if (class(SSAResult) == "try-error") {
                    errMsg <- displayErrMsg(SSAResult)
                    cat("ERROR: ", errMsg, "\n\n", sep = "")
                    next
               }
               
               if (diagplot) {
                    myfilename <- paste(outputPath, "diagPlotSEA_", respvar[i], sep = "")
                    if (!addEnv) myfilename <- paste(myfilename, "_", levels(data[,match(env, names(data))])[j], sep = "")
                    if (is.random) { myfilename <- paste(myfilename, "_random.png", sep = "")
                    } else { myfilename <- paste(myfilename, "_fixed.png", sep = "") }
                    
                    png(filename = myfilename)
                    if (addEnv) { graphSEAdiagplot(residuals(SSAResult$output$model), fitted(SSAResult$output$model), respvar[i], is.random) 
                    } else {
                         graphSEAdiagplot(residuals(SSAResult$output$model), 
                                          fitted(SSAResult$output$model), respvar[i], is.random, env, 
                                          envLevel = levels(data[,match(env, names(data))])[j]) 
                    }
                    dev.off()
               } # end stmt -- if (diagplot)
               
               if (heatmap) {
                    prevDir <- getwd()
                    setwd(outputPath)
                    heatmapCmd <- paste("Heatmap(SSAResult$residual", sep = "")
                    if (is.random) { heatmapCmd <- paste(heatmapCmd, ", genAs = 'random'" , sep = "")
                    } else { heatmapCmd <- paste(heatmapCmd, ", genAs = 'fixed'" , sep = "") }
                    heatmapCmd <- paste(heatmapCmd, ", row = '", row, "', column = '", column, "'", sep = "")
                    heatmapCmd <- paste(heatmapCmd, ", respvar = '", respvar[i], "', model = 'p-rep'", sep = "")
                    if (addEnv) { heatmapCmd <- paste(heatmapCmd, ", env = NULL)", sep = "") 
                    } else { heatmapCmd <- paste(heatmapCmd, ", env = '", env, "')", sep = "")  }
                    # heatmapCmd <- paste(heatmapCmd, ", env = '", env, "')", sep = "") 
                    
                    heatmap1 <- eval(parse(text = heatmapCmd))  
                    if (heatmap1[[1]]$site[[1]] == "unique") {
                    } else{
                         if (heatmap1[[1]]$site[[1]] == "empty") {
                         } else {
                              cat("Error: ", heatmap1[[1]]$site[[1]], "\n\n", sep = "")   
                         }
                    } 
                    setwd(prevDir)
               } # end stmt -- if (heatmap)
               
               if (histogram) {
                    histFilename <- paste(outputPath, "histSEA_", respvar[i], sep = "")
                    mainLabel <- paste("Histogram of ", respvar[i], sep = "")
                    if (!addEnv) { 
                         histFilename <- paste(histFilename, "_", levels(data[,match(env, names(data))])[j], sep = "")
                         mainLabel <- paste(mainLabel, " for ", env, " = ", levels(data[,match(env, names(data))])[j], sep = "")
                    }
                    histFilename <- paste(histFilename, ".png", sep = "")
                    
                    png(histFilename)
                    hist(tempData[,match(respvar[i], names(tempData))], xlab = respvar[i], main = mainLabel)
                    dev.off()
               } # end stmt -- if (histogram)
               
               # remove the variable env in the residual file if it was created 
               if (addEnv) {
                    SSAResult$residual <- SSAResult$residual[,-I(match(env, names(SSAResult$residual)))]
               }
               
               # remove the variables Test, Check & TestCheck in the residual file if check is excluded
               if (is.random) {
                    if (excludeCheck) {
                         SSAResult$residual <- SSAResult$residual[,-I(match(c("Test", "Check", "TestCheck"), names(SSAResult$residual)))]
                    }
               }
               result[[i]]$site[[j]]$ssaResult <- SSAResult
               
          } ## end stmt --- for (j in (1:nlevels(data[,match(env, names(data))])))
          
          # boxplot
          if (boxplot) {
               boxplotFilename <- paste(outputPath, "boxplotSEA_", respvar[i], ".png", sep = "")
               png(boxplotFilename)
               boxplot(data[,match(respvar[i], names(data))] ~ data[,match(env, names(data))], 
                       data = data, xlab = respvar[i], main = paste("Boxplot of ", respvar[i], sep = ""))
               dev.off()
          } # end stmt -- if (boxplot)
     } # end stmt --- for (i in (1:length(respvar)))
     
     detach("package:nlme")
     
     return(result)     
}

# ------------------------------------------------------------------
# ssaTestPrepMain
# Description: Performs Single Site Analysis for P-Rep Design
# Scripts Created by: Violeta I Bartolome for IRRI
# Scripts Modified by: Alaine A Gulles
# Functions Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------

ssaTestPrepMain <- function(data, respvar, geno, row, column, is.random = FALSE, excludeCheck = FALSE, 
                            checkList = NULL, moransTest = FALSE, 
                            spatialStruc = c("none", "CompSymm", "Gaus", "Exp", "Spher"), varCorr = FALSE) {
     
     # --- count the number of observations read and used --- #
     nobsRead <- nrow(data) # determine the number of observation read
     numPlots <- length(unique(data[,row]))*length(unique(data[,column]))
     tempData <- subset(data, subset = (is.na(data[,respvar]) == F))
     nobsUsed <- nrow(tempData) # determine the number of observation to be used in the analysis
     
     # --- compute response rate --- #
     responseRate <- nobsUsed/numPlots
     if (responseRate < 0.8) { stop("Too many missing observations. Cannot proceed with the analysis.") }
     
     # --- convert to factor:
     tempData[,match(geno, names(tempData))] <- factor(trimStrings(tempData[,match(geno, names(tempData))]))
     tempData[,match(row, names(tempData))] <- factor(trimStrings(tempData[,match(row, names(tempData))]))
     tempData[,match(column, names(tempData))] <- factor(trimStrings(tempData[,match(column, names(tempData))]))
     nLevelsGeno <- nlevels(tempData[,geno])
     levelsGeno <- levels(tempData[,geno])
     
     # --- determine if there is at least two entry with at least 2 levels
     capture.output(numObsPerEntry <- DescriptiveStatistics(data, var = respvar, grp = geno, statistics = "nnmiss"))
     EntryWithRep <- numObsPerEntry[numObsPerEntry[,"N_NonMissObs"] >= 2,]
     if (nrow(EntryWithRep)<= 1) { stop("At least two genotype should have at least two observation for the analysis to proceed.") }

     # --- print data summary
     cat("DATA SUMMARY:\n\n")
     cat("Number of observations read: ", nobsRead,"\n")
     cat("Number of observations used: ", nobsUsed,"\n\n")
     classInfo <- class.information(geno, tempData)
     print(classInfo)
     cat("\n\n")
     
     # compute harmonic means
     numreps <- data.frame(n = tapply(eval(parse(text = paste("tempData$", respvar,sep = ""))), eval(parse(text = paste("tempData$", geno, sep = ""))), FUN = length))
     numreps <- as.numeric(1/sapply(1/numreps, mean))
     
     # --- create a dummy variable 
     tempData <- data.frame(tempData, dummy = rep(1, nobsUsed)) 

     # if genotype is random and excluding checks in the variance estimation,
     # --- create a check and test variable if genotype is random
     if (is.random) {
          checkTestWarning <- FALSE
          if (excludeCheck) {
               # check if is checkList is not NULL and with correct entry
               if (is.null(checkList)) { excludeCheck <- FALSE 
               } else {
                    # check
                    checkList <- trimStrings(checkList)
                    correctCheckList <- all(checkList %in% levelsGeno)
                    if (!correctCheckList) {
                         deletedCheck <- checkList[!(checkList %in% levelsGeno)]
                         if (length(deletedCheck) == length(checkList)) {
                              checkList <- NULL
                              excludeCheck <- FALSE
                              checkTestWarning <- "All specified control/check levels are not present in this trial. Nothing is excluded in the estimation of genotypic variance."
                         } else {
                              checkList <- checkList[(checkList %in% levelsGeno)]
                              checkTestWarning <- paste("The following control/check level(s) is(are) not present in this trial and deleted from the list of genotype levels to exclude: ", paste(deletedCheck, collapse = ", "), sep = "")
                         }
                    }
               }
          } ## end stmt -- if (excludeCheck)
          
          if (excludeCheck) {
               # --- create Test and Check columns
               tempData$Check <- NULL
               tempData$Test <- NULL
               tempData[,"Check"] <- factor(ifelse(!is.na(match(tempData[,geno], checkList)), checkList[match(tempData[,geno], checkList)], 0))
               tempData[,"Test"] <- factor(ifelse(tempData[,"Check"] == 0, levelsGeno[match(tempData[,geno], levelsGeno)], 0))
               tempData$TestCheck <- factor(tempData[,"Test"]:tempData[,"Check"])
          } ## end stmt -- if (excludeCheck)
     }

     # --- construct the model
     if (!is.random) { myformula1 <- paste("fixed = ", respvar, " ~ ", geno, ", random = ~1|dummy", sep = "") # if genotype is fixed
     } else {
          # if genotype is random
          if (excludeCheck) { 
               myformula1 <- paste("fixed = ", respvar, " ~ Check, random = ~1|TestCheck", sep = "")
               myformula2 <- paste("fixed = ", respvar, " ~ Check, random = ~1|dummy", sep = "")
          } else {
               myformula1 <- paste("fixed = ", respvar, " ~ 1, random = ~1|", geno, sep = "")
               myformula2 <- paste("fixed = ", respvar, " ~ 1, random = ~1|dummy", sep = "")
          }
     }
     
     # --- call the lme function 
     result1 <- try(ssaModeling(data = tempData, row, column, formula = myformula1, spatialStruc), silent = TRUE)
     if (class(result1) == "try-error") {
          errMsg <- displayErrMsg(result1)
          stop(paste("ERROR: ", errMsg, sep = ""))
     }

     if (is.random) { 
          result2 <- ssaModeling(data = tempData, row, column, formula = myformula2, spatialStruc) 
          if (class(result2) == "try-error") {
               errMsg <- displayErrMsg(result2)
               stop(paste("ERROR: ", errMsg, sep = ""))
          }
     }

     corStruc <- c("none","CompSymm", "Gaus", "Exp", "Spher")
     corStrucLabel <- c("Zero Spatial Correlation","Compound Symmetry Correlation Structure", 
                       "Gaussian Spatial Correlation Structure", "Exponential Spatial Correlation Structure", 
                       "Spherical Spatial Correlation Structure")

     # --- get the variance component -- added by AAG
     varcompTable <- data.frame(Group = rownames(VarCorr(result1$model)),
                                Variance = as.numeric(VarCorr(result1$model)[,"Variance"]), 
                                Std.Dev. = as.numeric(VarCorr(result1$model)[,"StdDev"]))
     
     if (!is.random) { varcompTable <- varcompTable[-1,]
     } else {
          varcompTable[,"Group"] <- as.character(varcompTable[,"Group"])
          if (excludeCheck) { varcompTable[1,"Group"] <- "Test:Check" } else { varcompTable[1,"Group"] <- geno }
     }

     # print the variance component
     if (varCorr) {
          cat("VARIANCE COMPONENTS TABLE:\n\n")
          printDataFrame(varcompTable, border = FALSE, digits = 4)
          cat("\n\n")
     }
     
     # --- perform test of significance for genotypic effect --- #
     if (!is.random) {
          anovaTable <- anova(result1$model)
          anovaTable <- data.frame(rownames(anovaTable)[match(geno, rownames(anovaTable))],
                                   anovaTable[match(geno, rownames(anovaTable)),])
          rownames(anovaTable) <- NULL
          colnames(anovaTable)[1] <- "" 
          cat("TESTING FOR THE SIGNIFICANCE OF GENOTYPIC EFFECT:\n", sep = "")
          cat(corStrucLabel[match(result1$strucUsed, corStruc)], "\n\n", sep = "")
          cat("Analysis of Variance Table\n", sep = "")
          printDataFrame(anovaTable, border = FALSE, digits = 5)
          cat("\n\n")
     } else {
          cat("TESTING FOR THE SIGNIFICANCE OF GENOTYPIC EFFECT USING -2 LOGLIKELIHOOD RATIO TEST:\n", sep = "")
          cat(corStrucLabel[match(result1$strucUsed, corStruc)], "\n\n", sep = "")
          cat("Formula for Model1: ", myformula1,"\n",sep = "")
          cat("Formula for Model2: ", myformula2,"\n\n",sep = "")
          modelTable <- modelComparisonTable(result1$model, result2$model)
          modelTableTemp <- data.frame(rownames(modelTable),modelTable)
          colnames(modelTableTemp)[c(1, 6, 7)] <- c("", "Df", "Pr(>Chisq)")
          printDataFrame(modelTableTemp, border = FALSE, digits = 5)
          cat("\n\n")
          
          if (excludeCheck) {
               anovaTable <- anova(result1$model)
               anovaTable <- data.frame(rownames(anovaTable)[match("Check", rownames(anovaTable))],
                                        anovaTable[match("Check", rownames(anovaTable)),])
               rownames(anovaTable) <- NULL
               colnames(anovaTable)[1] <- "" 
               cat("TESTING FOR THE SIGNIFICANCE OF CHECK EFFECT:\n", sep = "")
               cat(corStrucLabel[match(result1$strucUsed, corStruc)], "\n\n", sep = "")
               cat("Analysis of Variance Table\n", sep = "")
               printDataFrame(anovaTable, border = FALSE, digits = 5)
               cat("\n\n")
          }
     }
     
     # --- estimate and print the lsmeans/predicted means
     if (!is.random) { 
          sumStat.table <- ssaEstimateLSMeans(result1$model)
          cat("GENOTYPE LSMEANS:\n\n")
          printDataFrame(sumStat.table, border = FALSE, digits = 4)
          cat("\n\n")
          colnames(sumStat.table)[match("LSMean", names(sumStat.table))] <- paste(respvar,"_Means", sep = "")
     } else { 
          if (excludeCheck) {
               sumStat.table <- coef(result1$model)
               sumStat.table <- data.frame(Test = t(as.data.frame(strsplit(rownames(sumStat.table), ":")))[,1],
                                           sumStat.table[,"(Intercept)"])
               sumStat.table <- subset(sumStat.table, Test != "0")
               rownames(sumStat.table) <- NULL
               colnames(sumStat.table) <- c(geno, "Means")
               attr(sumStat.table, "heading") <- paste("Predicted Means of ", geno,": \n", sep = "")
               cat("PREDICTED MEANS: \n\n")
               printDataFrame(sumStat.table, border = FALSE, digits = 4)
               cat("\n\n")
               colnames(sumStat.table)[match("Means", names(sumStat.table))] <- paste(respvar,"_PredictedMeans", sep = "")
               
               sumStat.lsmeans <- ssaEstimateLSMeans(result1$model)
               sumStat.lsmeans <- sumStat.lsmeans[sumStat.lsmeans[,"Check"] != "0",]
               rownames(sumStat.lsmeans) <- 1:nrow(sumStat.lsmeans)
               colnames(sumStat.lsmeans)[1] <- geno
               attr(sumStat.table, "heading") <- paste("LS Means of ", geno,": \n", sep = "")
               cat("CHECK/CONTROL LSMEANS: \n\n")
               printDataFrame(sumStat.lsmeans, border = FALSE, digits = 4)
               cat("\n\n")
               colnames(sumStat.lsmeans)[match("LSMean", names(sumStat.lsmeans))] <- paste(respvar,"_Means", sep = "")
          } else {
               sumStat.table <- coef(result1$model)
               sumStat.table <- data.frame(rownames(sumStat.table), sumStat.table)
               rownames(sumStat.table) <- NULL
               colnames(sumStat.table) <- c(geno, "Means")
               attr(sumStat.table, "heading") <- paste("Predicted Means of ", geno,": \n", sep = "")
               cat("PREDICTED MEANS: \n\n")
               printDataFrame(sumStat.table, border = FALSE, digits = 4)
               cat("\n\n")
               colnames(sumStat.table)[match("Means", names(sumStat.table))] <- paste(respvar,"_PredictedMeans", sep = "")
          }
     }
     
     # --- estimate heritability
     if (is.random) {
          if (excludeCheck) {
               herit <- varcompTable[varcompTable[,"Group"] == "Test:Check", "Variance"]/(varcompTable[varcompTable[,"Group"] == "Test:Check", "Variance"] + (varcompTable[varcompTable[,"Group"] == "Residual", "Variance"]/numreps)) 
          } else { 
               herit <- varcompTable[varcompTable[,"Group"] == geno, "Variance"]/(varcompTable[varcompTable[,"Group"] == geno, "Variance"] + (varcompTable[varcompTable[,"Group"] == "Residual", "Variance"]/numreps)) 
          }
          cat("HERITABILITY:\n\n")
          cat(formatC(herit, digits = 2, format = "f"),"\n\n\n")
          
     }
     
     # --- compute SED
     if (!is.random) SEDTable <- ssaSpatialSED(result1$model)
     
     # --- residual 
     tempDataResid <- data.frame(tempData, Residuals = residuals(result1$model))
     
          
     # --- COMPUTE MORAN'S I TO TEST IF SPATIAL CORRELATION EXISTS and print
     if (moransTest) {
          #tempDataResid <- data.frame(tempData, Residuals = residuals(result1$model))
          dataDists <- as.matrix(dist(cbind(tempDataResid[, row], tempDataResid[, column])))
          dataDistsInv <- 1/dataDists
          diag(dataDistsInv) <- 0
          moranResult <- data.frame(Moran.I(tempDataResid$Residuals, dataDistsInv))
          attr(moranResult, "heading") <- "Moran's I Autocorrelation Index"
          cat(toupper(attr(moranResult, "heading")), "\n\n", sep = "")
          printDataFrame(moranResult, border = FALSE, digits = 4)
          cat("\n\n")
     }
     
     if (!is.random) {
          colnames(tempDataResid)[ncol(tempDataResid)] <- paste(respvar, "_resid_fixed", sep = "")
          tempDataResid <- tempDataResid[,-I(match("dummy", names(tempDataResid)))]
          return(list(output = result1, formula = myformula1, varcomp.table = varcompTable, 
                      summary.statistics = sumStat.table, SEDTable = SEDTable,
                      residual = tempDataResid, numreps = numreps))
     } else {
          colnames(tempDataResid)[ncol(tempDataResid)] <- paste(respvar, "_resid_random", sep = "")
          tempDataResid <- tempDataResid[,-I(match("dummy", names(tempDataResid)))]
          if (excludeCheck) {
               return(list(output = result1, formula = myformula1, varcomp.table = varcompTable, 
                           summary.statistics = sumStat.table, lsmeans = sumStat.lsmeans, heritability = herit,
                           residual = tempDataResid, numreps = numreps, checkTestWarning = checkTestWarning))
          } else {
               return(list(output = result1, formula = myformula1, varcomp.table = varcompTable, 
                           summary.statistics = sumStat.table, heritability = herit,
                           residual = tempDataResid, numreps = numreps, checkTestWarning = checkTestWarning))
          }
          
     } 
} ## end of function --- ssaTestPrepMain

# ------------------------------------------------------------------
# ssaModeling
# Description: Perform the SSA for p-rep design using nlme package
# Scripts Created by: Violeta I Bartolome for IRRI
# Scripts Modified by: Alaine A Gulles
# Functions Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------

ssaModeling <- function(data, row, column, formula, spatialStruc) {
     # Single-Site analysis modelling for P-Rep design
     strucChoices <- c("CompSymm", "Gaus", "Exp", "Spher")
     spatialFormula <- c(paste("corCompSymm(0, form = ~ ", row, " + ", column, ")", sep = ""),
                         paste("corGaus(1, form = ~ ", row, " + ", column, ")", sep = ""),
                         paste("corExp(1, form = ~ ", row, " + ", column, ")", sep = ""),
                         paste("corSpher(form = ~ ", row, " + ", column, ")", sep = ""))
     
     mymodel <- list()
     resultAIC <- list()
     myformula <- list()
     errMsg <- list()
     
     for (i in (1:length(spatialStruc))) {
          if (spatialStruc[i] == "none") { myformula[[i]] <- formula
          } else {
               index <- match(spatialStruc[i], strucChoices)
               myformula[[i]] <- paste(formula, ", correlation = ", spatialFormula[index], sep = "")
          }
          mymodel[[i]] <- try(eval(parse(text = paste("lme(",myformula[[i]],", data, method = 'REML')", sep = ""))), silent = TRUE)
          if (class(mymodel[[i]]) == "try-error") { resultAIC[[i]] <- NA
          } else { resultAIC[[i]] <- AIC(mymodel[[i]]) }
     }
     
     if(length(spatialStruc) == 1) { 
          if (is.na(resultAIC[[1]])) { return(mymodel[[i]])
          } else { return(list(model = mymodel[[1]], formula = myformula[[1]], strucUsed = spatialStruc)) }
          
     } else {
          if (all(is.na(unlist(resultAIC)))) { return(mymodel[[i]])
          } else {
               index <- match(max(round(unlist(resultAIC), digits = 4), na.rm = TRUE), round(unlist(resultAIC), digits = 4))
               return(list(model = mymodel[[index]], formula = myformula[[i]], strucUsed = spatialStruc[index]))          
          }
     }
} ## end of function -- ssaModeling

# ------------------------------------------------------------------
# ssaSpatialSED
# Description: Compute the SED for SSA p-rep design
# Scripts Created by: Violeta I Bartolome for IRRI
# Scripts Modified by: Alaine A Gulles
# Functions Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------

ssaSpatialSED <- function(model) {
     # STANDARD ERROR OF DIFFERENCE
     covs <- as.matrix(vcov(model))
     vars <- diag(covs)
     vdiff <- outer(vars, vars, "+") - 2 * covs
     sed <- sqrt(vdiff[upper.tri(vdiff)])
     SEDTable <- data.frame(c("Minimum", "Average", "Maximum"), t(data.frame(min(sed), mean(sed), max(sed))))
     rownames(SEDTable) <- 1:nrow(SEDTable)
     colnames(SEDTable) <- c("","Estimate")
     cat("STANDARD ERROR OF THE DIFFERENCE (SED): \n\n", sep = "")
     printDataFrame(SEDTable, border = FALSE, digits = 4)
     cat("\n\n")
     return(SEDTable)
} ## end of function --- ssaSpatialSED


# ------------------------------------------------------------------
# ssaEstimateLSMeans
# Description: Estimate the LSMeans for SSA p-rep design
# Scripts Created by: Violeta I Bartolome for IRRI
# Scripts Modified by: Alaine A Gulles
# Functions Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------

ssaEstimateLSMeans <- function(model) {
     # --- estimate and print the lsmeans
     fixeffect <- data.frame(fixef(model))
     sumStat.table <- data.frame(levels(model$data[,names(model$contrast)[1]]), c(fixeffect[1,], fixeffect[-1,]+fixeffect[1,]))
     colnames(sumStat.table) <- c(names(model$contrast)[1], "LSMean")
     #cat("GENOTYPE LSMEANS\n\n")
     #printDataFrame(sumStat.table, border = FALSE, digits = 4)
     #cat("\n\n")
     return(sumStat.table)
} ## end of function -- ssaEstimateLSMeans


# ------------------------------------------------------------------
# ssaTestPrepResid
# Description: Consolidate the residuals for SSA p-rep design
# Scripts and Function Created by: Alaine A Gulles for IRRI
# Scripts and Function Modified by: Alaine A Gulles
# ------------------------------------------------------------------

ssaTestPrepResid <- function(object) {
     tempResidAll <- NULL
     for (i in (1:length(object))) {
          tempResidEnv <- NULL
          for (j in (1:length(object[[i]]$site))) {
               # create the dataframe of the residuals for all levels of env
               if (is.null(tempResidEnv)) { 
                    if (nrow(object[[i]]$site[[j]]$ssaResult$residual) > 0) { tempResidEnv <- object[[i]]$site[[j]]$ssaResult$residual }
               } else {
                    if (nrow(object[[i]]$site[[j]]$ssaResult$residual) > 0) { tempResidEnv <- rbind(tempResidEnv, object[[i]]$site[[j]]$ssaResult$residual) }
               }
               
          }
          
          # create data frame of the residuals for all response variables
          if (is.null(tempResidAll)) {
               if (nrow(tempResidEnv) > 0) { tempResidAll <- tempResidEnv }
          } else {
               if (nrow(tempResidEnv) > 0) { tempResidAll <- merge(tempResidAll, tempResidEnv, all = TRUE, sort = TRUE) }
          }
     }
     return(tempResidAll)
} # end of function ssaTestPrepResid

# ------------------------------------------------------------------
# ssaTestPrepSummary
# Description: Consolidate the summary statistics for SSA p-rep design
# Scripts and Function Created by: Alaine A Gulles for IRRI
# Scripts and Function Modified by: Alaine A Gulles
# ------------------------------------------------------------------


ssaTestPrepSummary <- function(object) {
     tempStatAll <- NULL
     for (i in (1:length(object))) {
          tempStatEnv <- NULL
          for (j in (1:length(object[[i]]$site))) {
               # create the dataframe of the summary statistics for all levels of env
               if (is.null(tempStatEnv)) { 
                    if (nrow(object[[i]]$site[[j]]$ssaResult$summary.statistics) > 0) { 
                         if (object[[1]]$site[[1]]$envLabel) { tempStatEnv <- object[[i]]$site[[j]]$ssaResult$summary.statistics 
                         } else {
                              tempStatEnv <- cbind(rep(object[[i]]$site[[j]]$envLevel, nrow(object[[i]]$site[[j]]$ssaResult$summary.statistics)),
                                                   object[[i]]$site[[j]]$ssaResult$summary.statistics)
                              colnames(tempStatEnv)[1] <- object[[i]]$site[[j]]$env
                         }                                      
                    }
               } else {
                    if (nrow(object[[i]]$site[[j]]$ssaResult$residual) > 0) { 
                         if (object[[1]]$site[[1]]$envLabel) { tempStatEnv <- rbind(tempStatEnv, object[[i]]$site[[j]]$ssaResult$residual) }
                         
                    } else {
                         tempStatEnv1 <- cbind(rep(object[[i]]$site[[j]]$envLevel, nrow(object[[i]]$site[[j]]$ssaResult$summary.statistics)),
                                               object[[i]]$site[[j]]$ssaResult$summary.statistics)
                         colnames(tempStatEnv1)[1] <- object[[i]]$site[[j]]$env
                         tempStatEnv <- rbind(tempStatEnv, tempStatEnv1)
                    }
               }
          }
          
          # create data frame of the residuals for all response variables
          if (is.null(tempStatAll)) {
               if (nrow(tempStatEnv) > 0) { tempStatAll <- tempStatEnv }
          } else {
               if (nrow(tempStatEnv) > 0) { tempStatAll <- merge(tempStatAll, tempStatEnv, all = TRUE, sort = TRUE) }
          }
     }
     return(tempStatAll)
} # end of function ssaTestPrepSummary


