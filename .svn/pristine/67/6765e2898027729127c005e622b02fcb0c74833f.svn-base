#----------------------------------------------------------------
# This function imports csv data and converts it to genind object.                      
# PBTools uses this for popgen analysis
#
# ARGUMENTS:
# filename - the full path and filename of the file to be imported
# extension	- extension of file, e.g. csv
# population - column name which corresponds to the population
# individual - column name which corresponds to the individual
# ploidyDegree - ploidy
#                                                     
# Script Created by: Nellwyn L. Sales
#----------------------------------------------------------------

importGenFile<-function(filename, extension=c("csv", "gtx", "str", "stru", "dat", "gen"), population=NULL, individual=NULL, ploidyDegree=2, sep=NULL) UseMethod("importGenFile")

importGenFile.default<-function(filename, extension=c("csv", "gtx", "str", "stru", "dat", "gen"), population=NULL, individual=NULL, ploidyDegree=2, sep=NULL) {
  library(adegenet)
  
  if (extension=="csv") {
    data1<-read.csv(filename, header=TRUE)
    
    if (is.null(sep)) {
      data1<-df2genind(data1[-c(match(population, names(data1)),match(individual, names(data1)))], ploidy=ploidyDegree, pop=data1[,population], ind.names=data1[,individual], sep="")
    } else {
      if (sep=="/") {
        data1<-df2genind(data1[-c(match(population, names(data1)),match(individual, names(data1)))], ploidy=ploidyDegree, pop=data1[,population], ind.names=data1[,individual], sep="/")
      }
    }
  }
  
  if (extension=="gtx") {
    data1<-read.genetix(filename)
  }
  
  if (extension=="str" || extension=="stru") {
    data1<-read.structure(filename)
  }
  
  if (extension=="dat") {
    data1<-read.fstat(filename)
  }
  
  if (extension=="gen") {
    data1<-read.genepop(filename)
  }
  
  return(data1)
}