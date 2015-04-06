#--------------------------------------------------------------------------------
# This function performs partial least squares regression and saves a biplot in the working directory
#
# ARGUMENTS:
# data - data frame of GxE means
# outFileName - path and filename of the text file to be created
# respvar - vector of the names of response variable
# genotype - string; name of genotype factor
# environment - string; name of environment factor
# covariateData - data frame of covariate data
# covariateEnvironment - string; name of environment factor in the covariate data; NULL if genotype characteristics are used as covariate
# covariateGenotype - string; name of genotype factor in the covariate data; NULL if environment characteristics are used as covariate
#
# Created by: Nellwyn L. Sales
#--------------------------------------------------------------------------------

plsRegression<-function(data, outFileName, respvar, genotype, environment, covariateData, covariateEnvironment, covariateGenotype) {
  
  #trim the strings
  outFileName<-trimStrings(outFileName)
  respvar<-trimStrings(respvar)
  genotype<-trimStrings(genotype)
  environment<-trimStrings(environment)
  if (!is.null(covariateEnvironment)) covariateEnvironment<-trimStrings(covariateEnvironment)
  if (!is.null(covariateGenotype)) covariateGenotype<-trimStrings(covariateGenotype)
  
  #add title
  capture.output(cat("\nPARTIAL LEAST SQUARES REGRESSION\n\n\n"),file=outFileName,append = FALSE)
  
  #call needed libraries
  library(pls)
  library (lme4)
  
  result<-list()
  
  #set genotype and environment to factors
  data[,environment]<-factor(data[,environment])
  data[,genotype]<-factor(data[,genotype])
  levelsGeno<-levels(data[,genotype])
  levelsEnv<-levels(data[,environment])
  
  commonLevels<-intersect(levelsGeno,levelsEnv)
  if (length(commonLevels)>0) {
    withCommonLevels<-TRUE
  } else {
    withCommonLevels<-FALSE
  }
  
  # --- if max length of the characters of levelsGeno or levelsEnv greater than 4 or Geno and Env have common levels, recode the levels
  if (max(nchar(levelsGeno))>4 || max(nchar(levelsEnv))>4 || withCommonLevels) {
    
    # --- recode genotype and environment levels --- #
    newCodingGeno<-data.frame(Genotype=levelsGeno, Code=paste("G",seq(1:length(levelsGeno)), sep=""))
    newCodingEnv<-data.frame(Environment=levelsEnv, Code=paste("E",seq(1:length(levelsEnv)), sep=""))
    
    suppressWarnings(result$newCodingGeno <- newCodingGeno)
    suppressWarnings(result$newCodingEnv <- newCodingEnv)
    recodedLevels <- TRUE
    temp.data <- data
    
    # --- attach the new labels to temp.data --- #
    temp.data$CodedGeno <- newCodingGeno$Code[match(temp.data[,genotype], newCodingGeno$Genotype)]
    temp.data$CodedEnv <- newCodingEnv$Code[match(temp.data[,environment], newCodingEnv$Environment)]
    
    # --- attach the new labels to covariateData
    if (!is.null(covariateEnvironment)) covariateData$CodedEnv <- newCodingEnv$Code[match(covariateData[,covariateEnvironment], newCodingEnv$Environment)]
    if (!is.null(covariateGenotype)) covariateData$CodedGeno <- newCodingGeno$Code[match(covariateData[,covariateGenotype], newCodingGeno$Genotype)]
    
    recodedLevels <- TRUE
    
  } else {
    temp.data <- data
    temp.data$CodedGeno <- temp.data[,match(genotype, names(temp.data))]
    temp.data$CodedEnv <- temp.data[,match(environment, names(temp.data))]
    if (!is.null(covariateEnvironment)) covariateData$CodedEnv <- covariateData[,match(covariateEnvironment, names(covariateData))]
    if (!is.null(covariateGenotype)) covariateData$CodedGeno <- covariateData[,match(covariateGenotype, names(covariateData))]
    
    recodedLevels <- FALSE
  }
  
  #set CodedEnv and CodedGeno to factors
  temp.data$CodedGeno <- factor(temp.data$CodedGeno)
  temp.data$CodedEnv <- factor(temp.data$CodedEnv)
  levelsGeno<-levels(temp.data[,"CodedGeno"])
  levelsEnv<-levels(temp.data[,"CodedEnv"])
  e <- length(levelsEnv)
  g <- length(levelsGeno)
  temp.dataAll<-temp.data
  
  for (i in 1:length(respvar)) {
    result$respvar[[i]]<-list()
    
    capture.output(cat("------------------------------\n"),file=outFileName,append = TRUE)
    capture.output(cat("RESPONSE VARIABLE: ", respvar[i], "\n"),file=outFileName,append = TRUE)
    capture.output(cat("------------------------------\n\n\n"),file=outFileName,append = TRUE)
    result$respvar[[i]]<-respvar[i]
    
    temp.data<-temp.dataAll
    
    if (nlevels(temp.data[,"CodedEnv"])>2) {
      
      #check if data contains one value for each env-geno combination
      rep<-tapply(temp.data[, respvar[i]] , temp.data[,c("CodedEnv", "CodedGeno")], function(x) length(which(!is.na(x))))
      
      #if data contains more than 1 observation per env-geno combination, compute for the means
      if (any(rep>1, na.rm=TRUE)) {
        meansWide<-as.data.frame(tapply(temp.data[, respvar[i]] , temp.data[,c("CodedEnv", "CodedGeno")], function(x) mean(x,na.rm=TRUE)))
        meansWide<-cbind(envTemp=rownames(meansWide), meansWide)
        colnames(meansWide)[1]<-"CodedEnv"
        temp.data<-data.frame(reshape(meansWide, direction="long", varying=colnames(meansWide)[2:ncol(meansWide)], v.names=respvar[i], timevar="CodedGeno", idvar="CodedEnv", times=colnames(meansWide)[2:ncol(meansWide)]), row.names = NULL)
        
        rep<-tapply(temp.data[, respvar[i]] , temp.data[,c("CodedEnv", "CodedGeno")], function(x) length(which(!is.na(x))))
        
        #set CodedEnv and CodedGeno to factors
        temp.data$CodedGeno <- factor(temp.data$CodedGeno)
        temp.data$CodedEnv <- factor(temp.data$CodedEnv)
      }
      
      #compute response rate
      responseRate<-1-((sum(rep==0, na.rm=TRUE)+sum(is.na(rep)))/(nlevels(temp.data[,"CodedEnv"])*nlevels(temp.data[,"CodedGeno"])))
      
      if (responseRate < 0.80) {
        
        capture.output(cat("*** \nERROR: Too many missing observations. Cannot proceed with the analysis.\n***\n\n\n"),file=outFileName,append = TRUE)
        suppressWarnings(result$respvar[[i]]$Error<-"ERROR: Too many missing observations. Cannot proceed with the analysis.")
        
      } else {
        #get only needed columns
        temp.data<-temp.data[, c("CodedEnv", "CodedGeno", respvar[i])]
        
        #run lm
        myformula1<-paste(respvar, " ~ CodedEnv + CodedGeno")
        model2 <- lm(formula(myformula1), data=temp.data)
        RESID <- resid(model2)
        
        #column bind temp.data with RESID
        RESID <- cbind(temp.data, RESID)
        
        #reshape RESID
        if (!is.null(covariateEnvironment)) {
          RESIDUAL <- reshape(data = RESID, v.names = "RESID", idvar = "CodedEnv", timevar = "CodedGeno", direction = "wide", drop=respvar[i])
          colnames(RESIDUAL) <-  gsub("RESID.", "G", colnames(RESIDUAL))
          
          Covariate <- covariateData[-match(covariateEnvironment, names(covariateData))]
          Covariate <- Covariate[order(Covariate[,"CodedEnv"]),]
          Covariate2 <-cbind(Covariate[match("CodedEnv", names(Covariate))],scale(Covariate[-match("CodedEnv", names(Covariate))]))
          
          residual2 <- cbind(RESIDUAL[1],scale(RESIDUAL[-1]))
          
          Varname <- names(Covariate[-match("CodedEnv", names(Covariate))])
          
          plsdata <- merge(residual2, Covariate2, by.x="CodedEnv", by.y="CodedEnv")
          
          Ydata <-  as.matrix(plsdata[,2:(g+1)])
          Xdata <- as.matrix(plsdata[,(g+2):ncol(plsdata)])
        }
        
        if (!is.null(covariateGenotype)) {
          RESIDUAL <- reshape(data = RESID, v.names = "RESID", idvar = "CodedGeno", timevar = "CodedEnv", direction = "wide", drop=respvar[i])
          colnames(RESIDUAL) <-  gsub("RESID.", "E", colnames(RESIDUAL))
          
          Covariate <- covariateData[-match(covariateGenotype, names(covariateData))]
          Covariate <- Covariate[order(Covariate[,"CodedGeno"]),]
          Covariate2 <-cbind(Covariate[match("CodedGeno", names(Covariate))],scale(Covariate[-match("CodedGeno", names(Covariate))]))
          
          residual2 <- cbind(RESIDUAL[1],scale(RESIDUAL[-1]))
          
          Varname <- names(Covariate[-match("CodedGeno", names(Covariate))])
          
          plsdata <- merge(residual2, Covariate2, by.x="CodedGeno", by.y="CodedGeno")
          
          Ydata <-  as.matrix(plsdata[,2:(e+1)])
          Xdata <- as.matrix(plsdata[,(e+2):ncol(plsdata)])
        }
        
        model1 <- plsr(Ydata~Xdata, scale=FALSE)
        capture.output(summary(model1, what="training"),file=outFileName,append = TRUE)
        capture.output(cat("\n"),file=outFileName,append = TRUE)
        #suppressWarnings(result$respvar[[i]]$percentExplained<-summary(model1, what="training")
        
        scores <- scores(model1)[,1:2]
        xloadin2 <- loadings(model1)[,1:3]
        yloadin2 <- Yloadings(model1)[,1:3]
        
        mscores <- as.matrix(scores)
        factor1 <- max(abs(mscores))
        mscores2 <- (1/factor1)*mscores
        scores3 <- as.data.frame(mscores2)
        colnames(scores3) <- c("dim1","dim2")  
        
        xload <- as.matrix(xloadin2)
        xload <- xload[,1:2]
        factor2 <- max(abs(xload))
        xload2 <- (1/factor2)*xload
        
        xload3 <- as.data.frame(xload2)
        colnames(xload3) <- c("dim1","dim2")
        
        yload <- as.matrix(yloadin2)
        yload <- yload[,1:2]
        factor3 <- max(abs(yload))
        yload2 <- (1/factor3)*yload
        
        yload3 <- as.data.frame(yload2)
        colnames(yload3) <- c("dim1","dim2")
        
        if (!is.null(covariateEnvironment)) {
          scores4 <- cbind(levelsEnv, scores3)
          names(scores4) <- c("name", "dim1", "dim2")
          scores4$type <- "ENV"
          rownames(scores4) <- NULL
          
          xload4 <- cbind(Varname, xload3)
          names(xload4) <- c("name", "dim1", "dim2")
          xload4$type <- "VAR"
          rownames(xload4) <- NULL
          
          yload4 <- cbind(levelsGeno, yload3)
          names(yload4) <- c("name", "dim1", "dim2")
          yload4$type <- "GENO"
          rownames(yload4) <- NULL
        }
        if (!is.null(covariateGenotype)) {
          scores4 <- cbind(levelsGeno, scores3)
          names(scores4) <- c("name", "dim1", "dim2")
          scores4$type <- "GENO"
          rownames(scores4) <- NULL
          
          xload4 <- cbind(Varname, xload3)
          names(xload4) <- c("name", "dim1", "dim2")
          xload4$type <- "TRAIT"
          rownames(xload4) <- NULL
          
          yload4 <- cbind(levelsEnv, yload3)
          names(yload4) <- c("name", "dim1", "dim2")
          yload4$type <- "ENV"
          rownames(yload4) <- NULL
        }
          
        # BIPLOT
        BIPLOT <- rbind(scores4, xload4, yload4)
        suppressWarnings(result$respvar[[i]]$scores<-BIPLOT)
        
        # save scores to csv file
        pcScoreFilename<-paste(getwd(),"/Scores_",respvar[i],".csv",sep = "")
        bplot_print<-data.frame(BIPLOT)
        write.csv(bplot_print, file=pcScoreFilename, row.names = FALSE)
        
        # set minimum and maximum values of x- and y-axis
        ymin <- min(BIPLOT$dim2) - abs(.07*(max(BIPLOT$dim2)-min(BIPLOT$dim2)))
        ymax <- max(BIPLOT$dim2) + abs(.07*(max(BIPLOT$dim2)-min(BIPLOT$dim2)))
        
        xmin <- min(BIPLOT$dim1) - abs(.07*(max(BIPLOT$dim1)-min(BIPLOT$dim1)))
        xmax <- max(BIPLOT$dim1) + abs(.07*(max(BIPLOT$dim1)-min(BIPLOT$dim1)))
        
        #start plotting
        png(filename = paste(getwd(),"/pls_biplot_",respvar,".png",sep = ""))
        par(cex=0.9)
        plot(xload4$dim1,xload4$dim2,cex=0.8, pch=" ", ylim=c(ymin, ymax), xlim=c(xmin,xmax), main="Biplot PLS", frame=TRUE,xlab="Factor 1", ylab="Factor 2")
        
        points(yload4$dim1,yload4$dim2,cex=0.8, pch=" ")
        points(scores4$dim1,scores4$dim2,cex=0.8, pch=" ")
        
        text(xload4$dim1, xload4$dim2, labels=xload4$name,col="black", cex=0.8)
        text(yload4$dim1, yload4$dim2, labels=yload4$name,col="red",cex=0.8)
        text(scores4$dim1, scores4$dim2, labels=scores4$name, col="blue", cex=0.8)
        
        abline(h=0,v=0,lty=2,col="black")
        if (!is.null(covariateEnvironment)) s <- seq(e)
        if (!is.null(covariateGenotype)) s <- seq(g)
        arrows(0, 0, 0.98*(scores4$dim1[s]), 0.99*scores4$dim2[s], col="blue",lwd=1,length=0.05)
        
        dev.off()
        
      } #end of else for if (responseRate < 0.80)
    
    } else {
      capture.output(cat("*** \nERROR: The environment factor should have at least three levels.\n***\n\n\n"),file=outFileName,append = TRUE)
      suppressWarnings(result$respvar[[i]]$Error<-"ERROR: The environment factor should have at least three levels.")
    } #end of else for if (nlevel(temp.data[,environment])>2)
    
    capture.output(cat("\n"),file=outFileName,append = TRUE)
    
  } #end of for (i in 1:length(respvar))
  
  # if levels were recoded and graphs were created, display codes used
  if (recodedLevels) {
    capture.output(cat("------------------------------\n"),file=outFileName,append = TRUE)
    capture.output(cat("\nCODES USED IN GRAPH:\n\n"),file=outFileName,append = TRUE)
    capture.output(newCodingGeno,file=outFileName,append = TRUE)
    capture.output(cat("\n"),file=outFileName,append = TRUE)
    capture.output(newCodingEnv,file=outFileName,append = TRUE)
  }
  
  capture.output(cat("\n==============================\n"),file=outFileName,append = TRUE)
  
  return(result)
}








