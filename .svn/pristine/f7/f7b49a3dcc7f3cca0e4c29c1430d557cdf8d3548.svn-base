# -----------------------------------------------------------------
# GENERAL LINEAR MODEL USING ASREML
# File Created by: Alaine A. Gulles 01.3.2014
#                  for International Rice Research Institute
# Note: uses the asreml package
# -----------------------------------------------------------------


GenLinearModelASR <- function(data, respvar, fixedmodel, randommodel, errorterm) UseMethod("GenLinearModelASR")

GenLinearModelASR.default <- function(data, respvar, fixedmodel, randommodel, errorterm) {

     # determine if asreml package is installed and loads it
     isPkgInstalled <- require(asreml)

     if (!isPkgInstalled) { stop("Error: ASReml package is not installed in your computer. Cannot perform analysis.") }      
     
     # if the package is installed
     nameData <- paste(deparse(substitute(data)))
     myformula <- paste("asreml(", respvar, " ~ ", fixedmodel, sep = "")
     if (!is.null(randommodel)) {
          myformula <- paste(myformula, ", random = ~ ", randommodel, sep = "")
     }
     myformula <- paste(myformula, ", rcov = ~ ", errorterm, sep = "")
     myformula <- paste(myformula, ", data = ", nameData, ", na.method.X = 'include')", sep = "")
     #myformula <- paste(myformula, ", data = eval(parse(text =", nameData,")), na.method.X = 'include')", sep = "")
     #formula(myformula)
     #asrmodel <- asreml(eval(parse(text = myformula)), data, na.method.X = "include")
     #capture.output(rcb.asr <- asreml(yield ~ Variety, random = ~ Rep, rcov = ~ units, 
     #                  data = newdata,
     #                  na.method.X = "include"))
     result <- capture.output(asrmodel <- eval(parse(text = myformula)))

     if (result[3] == "ABORT") stop("Error")

     # display the descriptive statistics of the response variable
     DescriptiveStatistics(data, var = respvar, statistics = c("nnmiss", "min", "max", "mean", "sd"))

     # testing the fixed effect in the model
     waldresult <- wald(asrmodel)
     if (rownames(waldresult)[1] == "(Intercept)") {
          # remove the ()
          rownames(waldresult)[1] <- "Intercept"
     }

     # print the table for testing the fixed effect in the model
     ConstructAOVTable(waldresult)

     # display the variance component of the model
     printDataFrame(cbind(item = rownames(summary(asrmodel)$varcomp), summary(asrmodel)$varcomp))

     # save/display fixed coefficients
     fixedcoef <- as.data.frame(asrmodel$coefficients$fixed)
     if (!is.na(match("(Intercept)", rownames(fixedcoef)))) {
          rownames(fixedcoef)[match("(Intercept)", rownames(fixedcoef))] <- "Intercept"
     }

     # save/display random coefficients, if any
     #if (!is.null(asrmodel$coefficients$random)) { cat("print the coefficient of the random component") }

     
     # display the 
     #addends <- strsplit(fixedmodel, split = " + ", fixed = TRUE)[[1]]
     #for (i in (1:length(addends))) {
     #     predVal <- predict(asrmodel, classify = addends[i])     
     #     predVal <- predict(asrmodel, classify = list("Rep", "Variety"))     
     #}
     

}