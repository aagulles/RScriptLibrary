#--------------------------------------------------------------------------------
# This function performs factorial regression and saves results in a text file
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

factorialRegression<-function(data, outFileName, respvar, genotype, environment, covariateData, covariateEnvironment, covariateGenotype) {
  
  options(show.signif.stars=FALSE)
  
  #trim the strings
  outFileName<-trimStrings(outFileName)
  respvar<-trimStrings(respvar)
  genotype<-trimStrings(genotype)
  environment<-trimStrings(environment)
  if (!is.null(covariateEnvironment)) covariateEnvironment<-trimStrings(covariateEnvironment)
  if (!is.null(covariateGenotype)) covariateGenotype<-trimStrings(covariateGenotype)
  
  #add title
  capture.output(cat("\nFACTORIAL REGRESSION\n\n\n"),file=outFileName,append = FALSE)
  
  #call needed libraries
  library (lme4)
  
  result<-list()
  
  #set genotype and environment to factors
  data[,environment]<-factor(data[,environment])
  data[,genotype]<-factor(data[,genotype])
  levelsGeno<-levels(data[,genotype])
  levelsEnv<-levels(data[,environment])
  e <- length(levelsEnv)
  g <- length(levelsGeno)
  
  for (i in 1:length(respvar)) {
    result$respvar[[i]]<-list()
    
    capture.output(cat("------------------------------\n"),file=outFileName,append = TRUE)
    capture.output(cat("RESPONSE VARIABLE: ", respvar[i], "\n"),file=outFileName,append = TRUE)
    capture.output(cat("------------------------------\n\n\n"),file=outFileName,append = TRUE)
    result$respvar[[i]]<-respvar[i]
    
    temp.data<-data
    
    if (nlevels(temp.data[,environment])>2) {
      
      #check if data contains one value for each env-geno combination
      rep<-tapply(temp.data[, respvar[i]] , temp.data[,c(environment, genotype)], function(x) length(which(!is.na(x))))
      
      #if data contains more than 1 observation per env-geno combination, compute for the means
      if (any(rep>1, na.rm=TRUE)) {
        meansWide<-as.data.frame(tapply(temp.data[, respvar[i]] , temp.data[,c(environment, genotype)], function(x) mean(x,na.rm=TRUE)))
        meansWide<-cbind(envTemp=rownames(meansWide), meansWide)
        colnames(meansWide)[1]<-environment
        temp.data<-data.frame(reshape(meansWide, direction="long", varying=colnames(meansWide)[2:ncol(meansWide)], v.names=respvar[i], timevar=genotype, idvar=environment, times=colnames(meansWide)[2:ncol(meansWide)]), row.names = NULL)
        
        rep<-tapply(temp.data[, respvar[i]] , temp.data[,c(environment, genotype)], function(x) length(which(!is.na(x))))
        
        #set environment and genotype to factors
        temp.data[, genotype] <- factor(temp.data[, genotype])
        temp.data[, environment] <- factor(temp.data[, environment])
      }
      
      #compute response rate
      responseRate<-1-((sum(rep==0, na.rm=TRUE)+sum(is.na(rep)))/(nlevels(temp.data[,environment])*nlevels(temp.data[,genotype])))
      
      if (responseRate < 0.80) {
        
        capture.output(cat("*** \nERROR: Too many missing observations. Cannot proceed with the analysis.\n***\n\n\n"),file=outFileName,append = TRUE)
        suppressWarnings(result$respvar[[i]]$Error<-"ERROR: Too many missing observations. Cannot proceed with the analysis.")
        
      } else {
        #get only needed columns
        temp.data<-temp.data[, c(environment, genotype, respvar[i])]
        if (!is.null(covariateEnvironment)) {
          Varname <- names(covariateData[-match(covariateEnvironment, names(covariateData))])
          Both_Data <- merge(temp.data,covariateData, by.x=environment, by.y=covariateEnvironment)
        }
        if (!is.null(covariateGenotype)) {
          Varname <- names(covariateData[-match(covariateGenotype, names(covariateData))])
          Both_Data <- merge(temp.data,covariateData, by.x=genotype, by.y=covariateGenotype)
        }
        
        ## Two way Analysis of Variance
        Both_Data[,environment] <- factor(Both_Data[,environment])
        Both_Data[,genotype] <- factor(Both_Data[,genotype])
        
        myformula1<-paste(respvar[i], " ~ ", environment, " + ", genotype, sep="")
        model0 <- lm(formula(myformula1), data=Both_Data)
        
        ## Factorial Regression with Env and Gen retained in the model
        if (!is.null(covariateEnvironment)) interactions<-paste(genotype, ":", Varname, sep="")
        if (!is.null(covariateGenotype)) interactions<-paste(environment, ":", Varname, sep="")
        
        addedInteractions<-paste(interactions,collapse=" + ")
        myformulaUpper<-paste(" ~ ", environment, " + ", genotype, " + ", addedInteractions, sep="")
        myformulaLower<-paste(" ~ ", environment, " + ", genotype, sep="")
        
        model.step <- step(model0, scope=list(upper = formula(myformulaUpper), lower= formula(myformulaLower)), direction="both", trace=FALSE, data = Both_Data)
                
        anovaTable<-anova(model.step)
        capture.output(anovaTable,file=outFileName,append = TRUE)
        suppressWarnings(result$respvar[[i]]$anovaTable<-anovaTable)
        
      } #end of else for if (responseRate < 0.80)
    
    } else {
      capture.output(cat("*** \nERROR: The environment factor should have at least three levels.\n***\n\n\n"),file=outFileName,append = TRUE)
      suppressWarnings(result$respvar[[i]]$Error<-"ERROR: The environment factor should have at least three levels.")
    } #end of else for if (nlevel(temp.data[,environment])>2)
    
    capture.output(cat("\n"),file=outFileName,append = TRUE)
    
  } #end of for (i in 1:length(respvar))
  
  capture.output(cat("\n==============================\n"),file=outFileName,append = TRUE)
  
  return(result)
}








