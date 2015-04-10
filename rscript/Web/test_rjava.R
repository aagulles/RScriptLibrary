##RJava

###ALL
start GBLUP
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
library(synbreed)
library(PBTools)
library(Matrix)
library(pbtgs)
sink(paste("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/GBLUPOut.txt", sep = ""))
gsGBLUP <- doGBLUP("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData2V.csv", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv", 3, FALSE, NULL, "t1", NULL, c("Trait.1"), NULL, TRUE, "random", NULL, 5, 2)
sink()
reached end.
TestGBLUP end
###
#GBLUP
#one trait
library(synbreed)
library(PBTools)
library(Matrix)
library(pbtgs)
sink(paste("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/GBLUPOut.txt", sep = ""))
gsGBLUP <- doGBLUP("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData2V.csv", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv", 3, FALSE, NULL, "t1", NULL, c("Trait.1"), NULL, TRUE, "random", NULL, 5, 2)
sink()

#two traits
start GBLUP
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
library(synbreed)
library(PBTools)
library(Matrix)
library(pbtgs)
sink(paste("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/GBLUPOut.txt", sep = ""))
gsGBLUP <- doGBLUP("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData2V.csv", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv", 3, FALSE, NULL, "t1", NULL, c("Trait.1", "Trait2"), NULL, TRUE, "random", NULL, 5, 2)
sink()
reached end.
TestGBLUP end

#Test2
start GBLUP
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
library(synbreed)
source("E:/StarPbtools/GS/script/BLUP_synbreed_gv.R")
library(Matrix)
source("E:/StarPbtools/GS/script/BLUP_synbreed_cv.R")
source('E:/StarPbtools/GS/script/doGBLUP.R')
source("E:/StarPbtools/GS/script/marker_relationship.R")
source('E:/StarPbtools/GS/script/createGSPlots.R')
sink(paste("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/GBLUPOut.txt", sep = ""))
gsGBLUP <- doGBLUP("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData2V.csv", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv", 3, FALSE, NULL, "t1", NULL, c("Trait.1"), NULL, FALSE, "BL", "random", NULL, 5, 2)
sink()
reached end.
TestGBLUP end


outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
pheno_file = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData2V.csv"
geno_file = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv"
markerFormat = 3
importRel = FALSE
rel_file = NULL
rMatType = "t1"
map_file = NULL
traitNames = c("Trait.1", "Trait.2")
covariates = NULL
doCV = FALSE #TRUE #
varCompEst = "BL"
samplingStrat = "random"
popStruc_file = NULL
nfolds = 5
nrep = 2

#graph, after running doGBLUP
gpData1 = gsObject
relData = gen_RM$relMat
# traitNames
# outputPath, gsObject, gen_RM$relMat, traitNames


#Test1
start GBLUP
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  library(synbreed)
source("E:/StarPbtools/GS/script/BLUP_synbreed_gv.R")
library(Matrix)
source("E:/StarPbtools/GS/script/BLUP_synbreed_cv.R")
source('E:/StarPbtools/GS/script/doGBLUP.R')
source("E:/StarPbtools/GS/script/marker_relationship.R")
source('E:/StarPbtools/GS/script/createGSPlots.R')
sink(paste("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/GBLUPOut.txt", sep = ""))
gsGBLUP <- doGBLUP("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData2V.csv", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv", 3, FALSE, NULL, "t1", NULL, c("Trait.1","Trait2"), NULL, FALSE, "BL", "random", NULL, 5, 2)
sink()
reached end.
TestGBLUP end



outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
sink(paste(outputPath, "GBLUPOut.txt", sep = ""))
gsGBLUP <- doGBLUP("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData.csv", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv", 3, FALSE, NULL, "t1", NULL, c("Trait.1"), NULL, TRUE, "BL", "random", 5, 2)
sink()


doGBLUP <- function(outputPath, pheno_file, geno_file = NULL, markerFormat = c(1, 2, 3), 
                    importRel = FALSE, rel_file = NULL, rMatType = c("t1", "t2", "t3", "t4"), 
                    map_file = NULL, # ped_file = NULL, #peFormat = NULL, #data quality check options, ...,
                    traitNames, covariates = NULL, doCV = FALSE, varCompEst = c("BL", "BRR"), samplingStrat = c("random","within popStruc"), nfolds = 2, nrep = 1) {
  gsGBLUP <- doGBLUP(, , , , , , , , , , , , , 2, 1)
  
doGBLUP <- function(
outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
pheno_file = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/phenoData.csv"
geno_file = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/synbreed.csv"
markerFormat = 3
importRel = FALSE
rel_file = NULL
rMatType = "t1"
map_file = NULL
traitNames = c("Trait.1")
covariates = NULL
doCV = FALSE #TRUE #
varCompEst = "BL"
samplingStrat = "random"
popStruc_file = NULL
nfolds = 5
nrep = 2

rm(nrep)
rm(list=c("nrep", "rMatType"))
rm(list=c("map_file", "rMatType"))
rm(list=c(nrep, rMatType))
# for cv:
traitName = "Trait.1"

#GSDataImputation
start GSDataImputation
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
source("E:/StarPbtools/GS/script/GSDataImputation.R")
gsDataImputation <- GSDataImputation("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/genoData_qc.csv", "random", NULL, NULL, "synbreed")
reached end.
TestGSDataImputation end

maf =
gsDataImputation <- GSDataImputation(, , , , NULL, "synbreed")
GSDataImputation <- function(outputPath, geno_file, impType = c("random", "family"), pheno_file = NULL, familyTrait = NULL, packageFormat = c("synbreed", "rrBLUP", "BGLR")) {
  GSDataImputation <- function(
    outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
    geno_file = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/genoData_qc.csv"
    impType = "random"
    pheno_file = NULL; familyTrait = NULL
    packageFormat = "synbreed"
    

#GSDataCheck
start GSDataCheck
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  source("E:/StarPbtools/GS/script/GSDataCheck.R")
# gsDataCheckOut <- GSDataCheck("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets//genoData.csv", 3, 0.1, 0.05, 0.9)
gsDataCheckOut <- GSDataCheck("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/ori_file_type3_ed.csv", 3, 0.1, 0.05, 0.9)
gsDataCheckOut2 <- GSDataCheck("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/maize_geno2.csv", 5, 0.1, 0.05, 0.9)
reached end.
TestGSDataCheck end


# gsDataCheckOut <- GSDataCheck(, , )
# GSDataCheck <- function(outputPath, geno_file, type, nmiss = 0.1, maf = 0.05, cor_threshold = 0.90) 
# GSDataCheck <- function(
  outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
  geno_file = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets//genoData.csv"
  type = 3
  nmiss = 0.1
  maf = 0.05
  cor_threshold = 0.9
  


#GSDataPrep
start GSDataPrep
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  source("E:/StarPbtools/GS/script/GSDataPrep.R")
source('E:/StarPbtools/GS/script/trimStrings.R')
gsDataPrepOut <- GSDataPrep("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/StarPbtools/GS/data/maize_phenoCov2.csv", "E:/StarPbtools/GS/data/maize_geno2.csv", "E:/StarPbtools/GS/data/maize_map2.csv", "E:/StarPbtools/GS/data/relMatFile2.m", "csv", "csv", "csv", "ttxt")
reached end.
TestGSDataPrep here

# start GSDataPrep
# resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
#   source("E:/StarPbtools/GS/script/GSDataPrep.R")
# source('E:/StarPbtools/GS/script/trimStrings.R')
# gsDataPrepOut <- GSDataPrep("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/StarPbtools/GS/data/maize_phenoCov2.csv", "E:/StarPbtools/GS/data/maize_geno2.csv", "E:/StarPbtools/GS/data/maize_map2.csv", "E:/StarPbtools/GS/data/relMatFile2.m", "csv", "csv", "csv", "ttxt")
# reached end.
# TestGSDataPrep here


gsDataPrepOut <- GSDataPrep(
  , , , , , , , , )
str(gsDataPrepOut)
GSDataPrep(
  outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
    pheno_file = "E:/StarPbtools/GS/data/maize_phenoCov2.csv"
    geno_file = "E:/StarPbtools/GS/data/maize_geno2.csv"
    map_file = "E:/StarPbtools/GS/data/maize_map2.csv"
    rel_file = "E:/StarPbtools/GS/data/relMatFile2.m"
    #                        ped_file = NULL, 
    pFormat = "csv"
    gFormat = "csv"
    mFormat = "csv"
    rFormat = "ttxt"
    #                        , peFormat = NULL
)

# GenSim

start GenSim
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  outFileName: markerRelation_realized
library(synbreed)
library(pedigreemm)
library(PBTools)
library(pbtgs)
genSimOut <- doGenSim("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/ori_file_type3_ed.csv", "csv", FALSE, "realized", "markerRelation_realized", 3)
reached end.
TestGenSim end


#marker-based
start GenSim
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  outFileName: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/markerRelation_realized
library(synbreed)
library(pedigreemm)
source('E:/StarPbtools/GS/script/trimStrings.R')
source("E:/StarPbtools/GS/script/GSDataCheck.R")
source("E:/StarPbtools/GS/script/marker_relationship.R")
source("E:/StarPbtools/GS/script/pedigree_relationship.R")
source("E:/StarPbtools/GS/script/doGenSim.R")
genSimOut <- doGenSim("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/ori_file_type3_ed.csv", "csv", FALSE, "realized", "markerRelation_realized", 3)
reached end.
TestGenSim end



any(geno_file == 1)
doGenSim(
  outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
  fileName = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/ori_file_type3_ed.csv"
  fileFormat = "csv"
  doPedigree = FALSE
  relType = "realized" #"sm-smin"
  outFileName = "markerRelation_realized" #"E:/App Files/workspace_Juno/RJavaManager/sample_datasets/markerRelation_realized" #sm-smin" 
  markerFormat = 3
  

#pedigree-based
# start GenSim
# resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
#   outFileName: pedigreeRelation_dom
library(synbreed)
library(pedigreemm)
# library(GeneticsPed)
source('E:/StarPbtools/GS/script/trimStrings.R')
source("E:/StarPbtools/GS/script/GSDataCheck.R")
source("E:/StarPbtools/GS/script/marker_relationship.R")
source("E:/StarPbtools/GS/script/pedigree_relationship.R")
source("E:/StarPbtools/GS/script/doGenSim.R")
genSimOut <- doGenSim("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/pedFile.csv", "csv", TRUE, "dom", "pedigreeRelation_dom", 1)
# reached end.
# TestGenSim end

start GenSim
resultFolderPath: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  outFileName: E:/App Files/workspace_Juno/RJavaManager/sample_datasets/pedigreeRelation_dom
library(synbreed)
library(pedigreemm)
# library(GeneticsPed)
source('E:/StarPbtools/GS/script/trimStrings.R')
source("E:/StarPbtools/GS/script/GSDataCheck.R")
source("E:/StarPbtools/GS/script/marker_relationship.R")
source("E:/StarPbtools/GS/script/pedigree_relationship.R")
source("E:/StarPbtools/GS/script/doGenSim.R")
genSimOut <- doGenSim("E:/App Files/workspace_Juno/RJavaManager/sample_datasets/", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/pedFile.csv", "csv", TRUE, "dom", "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/pedigreeRelation_dom", 1)
genSimOut <- doGenSim(
  , , , , , , 1)

doGenSim <- function(
  outputPath = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/"
  fileName = "E:/App Files/workspace_Juno/RJavaManager/sample_datasets/pedFile.csv"
  fileFormat = "csv"
  doPedigree = TRUE
  relType = "dom"
  outFileName = "pedigreeRelation_dom"  #E:/App Files/workspace_Juno/RJavaManager/sample_datasets/
  markerFormat = 1
  
  gen_RM <- pedigree_relationship(outputPath, ped_file = fileName, fFormat = fileFormat, type = relType, outputName = outFileName)
  gen_RM <- pedigree_relationship(
    outputPath, 
    ped_file = fileName, 
    fFormat = fileFormat, 
    type = relType, 
    outputName = outFileName
  
reached end.
TestGenSim end




# example with maize data and imputing by family
data(maize)
# first only recode alleles
maize.coded <- codeGeno(maize,label.heter=NULL)
maize.coded2 <- codeGeno(maize,label.heter=NULL)
# set 200 random chosen values to NA
set.seed(123)
ind1 <- sample(1:nrow(maize.coded $geno),200)
ind2 <- sample(1:ncol(maize.coded $geno),200)
original <- maize.coded2$geno[cbind(ind1,ind2)]
View(maize.coded2$geno[cbind(ind1,ind2)])
View(maize.coded2$geno)
maize.coded2$geno[cbind(ind1,ind2)] <- NA
# imputing of missing values by family structure
maize.imputed <- codeGeno(maize.coded2,impute=TRUE,impute.type="family",label.heter=NULL)
maize.im <- codeGeno(maize.coded2,impute=FALSE,label.heter=NULL)
summary(maize.coded2)
summary(maize.im)

str(maize)
View(maize.coded$covar)
str(maize.coded)
str(maize.imputed)


maize$c
str(maize)
maize.coded <- codeGeno(maize,label.heter=NULL)
str(maize.coded)
maize.imputed <- codeGeno( maize.coded,impute=TRUE,impute.type="family",label.heter=NULL)

str(mice)



data(maize)
# first only recode alleles
maize.coded <- codeGeno(maize,label.heter=NULL)
maize.coded2 <- codeGeno(maize,label.heter=NULL)
maize.imputed <- codeGeno(maize.coded2,impute=TRUE,impute.type="family",label.heter=NULL)
