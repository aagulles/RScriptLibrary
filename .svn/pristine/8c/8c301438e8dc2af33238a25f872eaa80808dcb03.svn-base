#----------------------------------------------------------------
# This function performs Test on Hardy-Weinberg Equilibrium. This uses the HWE.test.genind function of adegenet package.                      
# PBTools uses this for popgen analysis 
#
# ARGUMENTS:
# genindObject - genind object
# display - value is "full" if the list of complete tests will be displayed; value is "matrix" if matrix of p-values will be displayed
#                                                     
# Script Created by: Nellwyn L. Sales
#----------------------------------------------------------------

HWETest<-function(genindObject, display=c("full", "matrix")) UseMethod("HWETest")

HWETest.default<-function(genindObject, display=c("full", "matrix")) {
  library(adegenet)
  
  if (display=="matrix") {
    result1<-HWE.test.genind(genindObject, res.type="matrix")
    summaryTable<-format(round(as.data.frame(result1),4), digits=4, nsmall=4, scientific=FALSE)
    
    #cat("HARDY-WEINBERG EQUILIBRIUM TEST\n\n")
    #cat("MATRIX OF P-VALUES:\n\n")
    #print(summaryTable)
  }
  
  if (display=="full") {
    result2<-HWE.test.genind(genindObject, res.type="full")
    
    summaryTable<-NULL
    for (i in 1:length(result2)) {
      locusList<-eval(parse(text=paste("result2$", names(result2)[i], sep="")))
      for (j in 1:length(locusList)) {
        popnList<-eval(parse(text=paste("locusList$", names(locusList)[j], sep="")))
        newRow<-data.frame(Locus=names(result2)[i], Population=names(locusList)[j], Chi_value=popnList$statistic, df=popnList$parameter, p_value=popnList$p.value)
        summaryTable<-rbind(summaryTable, newRow)
      }
    }
    
    #format summaryTable
    rownames(summaryTable)<-NULL
    colnames(summaryTable)<-c("Locus", "Population", "Chisq Value", "df", "Pr(>Chisq)")
    summaryTable[, "Chisq Value"]<- formatC(as.numeric(format(summaryTable[, "Chisq Value"], scientific=FALSE)), format="f")
    summaryTable[, "Pr(>Chisq)"]<- formatC(as.numeric(format(summaryTable[, "Pr(>Chisq)"], scientific=FALSE)), format="f")
    
    #cat("HARDY-WEINBERG EQUILIBRIUM TEST\n\n")
    #cat("SUMMARY OF TESTS:\n\n")
    #print(summaryTable)
  }
    
  return(summaryTable)
  
}