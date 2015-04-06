#----------------------------------------------------------------
# Description: This function constructs model comparison table when model1 and model2 are fitted using REML=TRUE
#              and the term to be tested is random.
# Output: model comparison table
# Author: Nellwyn L. Sales  08.07.2013
# Script Modified by: Alaine A. Gulles 07.24.2014
#
# Arguments:
#  model1 - resulting object when lmer function is called using REML=TRUE (full model)
#  model2 - resulting object when lmer function is called using REML=TRUE (reduced model)
#----------------------------------------------------------------

modelComparisonTable<-function (model1, model2) UseMethod("modelComparisonTable")

modelComparisonTable.default<-function (model1, model2) {
  
     models.table=NULL
     
     # the following codes was hidden by AAGulles 07.24.2014 for R Version 3.0.2
     #models.table<-rbind(summary(model2)@AICtab[,1:3],summary(model1)@AICtab[1:3])
     #rownames(models.table)<-c("Model2", "Model1")
     #Chisq<-c("", round(2*abs(models.table[1,3]-models.table[2,3]), digits=4))
     #Df<-c("", 1)
     #probValue<-pchisq(2*abs(models.table[1,3]-models.table[2,3]), 1, lower.tail=FALSE)
     #pvalue<-c("",formatC(as.numeric(format(probValue, scientific=FALSE)), format="f"))
     #models.table<-cbind(models.table, Chisq, Df, pvalue)
     #colnames(models.table) <- c("AIC", "BIC", "logLik", "Chisq", "Df", "Pr(>Chisq)")
     
     # the following codes was modified by AAGulles 07.24.2014 for R Version 3.0.2
     models.table <- cbind(rbind(AIC(model2), AIC(model1)),
                           rbind(BIC(model2), BIC(model1)),
                           rbind(logLik(model2), logLik(model1)),
                           rbind(NA, round(2* abs(logLik(model2) - logLik(model1)), digits = 4)),
                           rbind(NA, 1),
                           rbind(NA, as.numeric(format(pchisq(2*abs(logLik(model2) - logLik(model1)), 1, lower.tail=FALSE), scientific = FALSE))))
     colnames(models.table) <- c("AIC", "BIC", "logLik", "Chisq", "Chi Df", "Pr(>Chisq)")
     rownames(models.table) <- c("Model2", "Model1")
     
     return(models.table)
}