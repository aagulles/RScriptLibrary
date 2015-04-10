###########
#  GSDataPrep
#' Function for preparing/matching data for genomic selection
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param pheno_file name of phenotypic file
#' @param geno_file name of genotypic file
#' @param map_file name of map file
#' @param rel_file name of relationship matrix file
#  ped_file = NULL, 
#' @param pFormat format of phenotypic file
#' @param gFormat format of genotypic file
#' @param mFormat  format of map file
#' @param rFormat  format of relationship matrix file
#  , peFormat = NULL
#  Output
#   isNewPhenoCreated - logical; whether new phenotype file is created or not
#   isNewGenoCreated logical; whether new genotype file is created or not
#   isNewRelMatCreated - logical; whether new relationship matrix file is created or not
#   isNewMapCreated - logical; whether new map file is created or not
# #     isNewPedCreated = isNewPedCreated,
#   phenoData - R data containing phenotype information
#   genoData - R data containing genotype information
#   mapData - R data containing map information
#   relData - R data containing relationship matrix
###########

# # pheno_file = "E:/StarPbtools/GS/data/phenotype_ed.csv"
# 
# # geno_file = "E:/StarPbtools/GS/data/ori_file_type3_ed.csv"
# # rel_file = "E:/StarPbtools/GS/data/ori_file_type3-synbreed-markerRelation_sm-smin.m"
# # pFormat = "csv"
# # gFormat = "csv" 
# 
# outputPath <- "E:/StarPbtools/GS/sample_output/"
# pheno_file = "E:/StarPbtools/GS/data/maize_phenoCov2.csv"
# # geno_file = "E:/StarPbtools/GS/data/maize_geno2.csv"
# geno_file = "E:/StarPbtools/GS/data/maize_geno2b.csv"
# # cov_file = "E:/StarPbtools/GS/data/maize_cov.csv"
# map_file = "E:/StarPbtools/GS/data/maize_map2.csv"
# #create rel_file
# rel_file = "E:/StarPbtools/GS/data/relMatFile2.m"
# # ped_file = "E:/StarPbtools/GS/data/maize_ped.csv" # NULL
# pFormat = "csv"
# gFormat = "csv" 
# # cFormat = "csv" 
# mFormat = "csv" 
# rFormat = "csv" 
# peFormat = "csv"

###GS data preparation/matching

GSDataPrep <- function(outputPath, pheno_file, geno_file = NULL,  map_file = NULL, rel_file = NULL, 
#                        ped_file = NULL, 
                       pFormat, gFormat = NULL, mFormat = NULL, rFormat = NULL
#                        , peFormat = NULL
                       ) {

  isNewPhenoCreated <- NULL
  isNewGenoCreated <- NULL
  isNewRelMatCreated <- NULL
  isNewMapCreated <- NULL
#   isNewPedCreated <- NULL
  
  genoData <- NULL
  relData <- NULL
  mapData <- NULL
#   pedData <- NULL

  PG_diff_gid <- 0
  GP_diff_gid <- 0
  PR_diff_gid <- 0
  RP_diff_gid <- 0
  
  formatsep <- function(fformat) {
    switch(fformat, 
         csv = ",",
         ctxt = ",",
         stxt = " ",
         ttxt = "\t",
         sctxt = ";")
  }

  pSep = formatsep(pFormat)
  phenoData <- read.table(pheno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = pSep)    
  rownames(phenoData) <- trimStrings(rownames(phenoData))
  P_gid <- rownames(phenoData)

  if (!is.null(geno_file)) {
    gSep = formatsep(gFormat)
    genoData <- read.table(geno_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = gSep)    
    rownames(genoData) <- trimStrings(rownames(genoData))
    colnames(genoData) <- trimStrings(colnames(genoData))
    G_gid <- rownames(genoData)
    G_mid <- colnames(genoData)
  }

  if (!is.null(rel_file)) {
    #check that rel. matrix is symmetric matrix and row and column names match
#     relData <- as.matrix(read.table(rel_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = "\t"))
#     relData <- as.matrix(read.table(rel_file, header = TRUE, check.names = FALSE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = "\t"))
    rSep = formatsep(rFormat)
    relData <- as.matrix(read.table(rel_file, header = TRUE, check.names = FALSE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = rSep))
    # View(relData)
    # class(relData)
    # isSymmetric(relData)
    rownames(relData) <- trimStrings(rownames(relData))
    colnames(relData) <- trimStrings(colnames(relData))
    Rrow_gid <- rownames(relData)
    Rcol_gid <- colnames(relData)
    if (sum(Rrow_gid != Rcol_gid) != 0 || (!isSymmetric(relData))) stop("Invalid relationship matrix.")
    #TODO: check if contents are valid: diag: 1 or 2, off-diag: < respective diag
  }  

  if (!is.null(map_file)) {
    mSep = formatsep(mFormat)
    mapData <- read.table(map_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = mSep)    
    rownames(mapData) <- trimStrings(rownames(mapData))
    M_mid <- rownames(mapData)
  }

#   if (!is.null(ped_file)) {
#     peSep = formatsep(peFormat)
#     pedData <- read.table(ped_file, header = TRUE, na.strings = c("NA","."), row.names = 1, blank.lines.skip=TRUE, stringsAsFactors=FALSE, sep = peSep)    
#     rownames(pedData) <- trimStrings(rownames(pedData))
#     Pe_gid <- as.character(pedData$ID) #rownames(pedData)
#   }

  # trim strings of IDs
#   P_gid <- trimStrings(P_gid)
#   G_gid <- trimStrings(G_gid)
#   G_mid <- trimStrings(G_mid)
#   Rrow_gid <- trimStrings(Rrow_gid)
#   Rcol_gid <- trimStrings(Rcol_gid)
#   M_mid <- trimStrings(M_mid)
#   Pe_gid <- trimStrings(Pe_gid)

  #match pheno and geno data
  if (!is.null(geno_file)) {
    ##genotypes in phenoData w/c are not in genoData; display (if any)
    PG_diff_gid <- setdiff(P_gid, G_gid)
    
    ##genotypes in genoData w/c are not in phenoData; display (if any)
    #   GP_diff_gid <- as.character(as.matrix(setdiff(paste(G_gid), P_gid)))
    #   GP_diff_gidNoSpace<-GP_diff_gid[which(GP_diff_gid!="")]
    GP_diff_gid <- setdiff(G_gid, P_gid)
    
    ###reduce (if needed) phenoData, sort genotypes as in genoData
    #   phenoData_red <- merge(phenoData, data.frame(G_gid), by.x = "row.names", by.y = "G_gid", sort = FALSE)  
    #   rownames(phenoData_red) <- phenoData_red[,1]
    #   phenoData_red <- phenoData_red[,-1]
    #   View(phenoData_red)
    isNewPhenoCreated <- FALSE
    if (length(PG_diff_gid) != 0) {
      phenoData_red <- merge(phenoData, data.frame(G_gid), by.x = "row.names", by.y = "G_gid", sort = FALSE)  
      rownames(phenoData_red) <- phenoData_red[,1]
      phenoData_red <- phenoData_red[,-1]
      #save new pheno file
      #     write.table(phenoData_red, file = paste(outputPath, "/newPhenoData.csv", sep = ""), sep = ",", quote = FALSE, row.names = TRUE, col.names = TRUE)
#       write.csv(phenoData_red, file = paste(outputPath, "/newPhenoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
      isNewPhenoCreated <- TRUE
      phenoData <- phenoData_red
      P_gid <- rownames(phenoData)
      #     PG_diff_gid
    } 
    write.csv(phenoData, file = paste(outputPath, "/phenoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
    
    ###reduce (if needed) genoData
    #   genoData_red <- merge(genoData, data.frame(P_gid), by.x = "row.names", by.y = "P_gid", sort = FALSE)
    #   rownames(genoData_red) <- genoData_red[,1]
    #   genoData_red <- genoData_red[,-1]
    isNewGenoCreated <- FALSE
    if (length(GP_diff_gid) != 0) {
      genoData_red <- merge(genoData, data.frame(P_gid), by.x = "row.names", by.y = "P_gid", sort = FALSE)
      rownames(genoData_red) <- genoData_red[,1]
      genoData_red <- genoData_red[,-1]
#       write.csv(genoData_red, file = paste(outputPath, "/newGenoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
      isNewGenoCreated <- TRUE
      genoData <- genoData_red
      G_gid <- rownames(genoData)
    }
    write.csv(genoData, file = paste(outputPath, "/genoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
  }
  

  #match pheno data and rel matrix
  if (!is.null(rel_file)) {
    ##genotypes in phenoData w/c are not in rel matrix; display (if any)
    PR_diff_gid <- setdiff(P_gid, Rrow_gid)
    
    ##genotypes in rel matrix w/c are not in phenoData; display (if any)
    RP_diff_gid <- setdiff(Rrow_gid, P_gid)
    
    ###reduce (if needed) phenoData, sort genotypes as in rel matrix
    #   phenoData_red <- merge(phenoData, data.frame(Rrow_gid), by.x = "row.names", by.y = "Rrow_gid", sort = FALSE)  
    #   rownames(phenoData_red) <- phenoData_red[,1]
    #   phenoData_red <- phenoData_red[,-1]
    #   View(phenoData_red)
    if (is.null(isNewPhenoCreated)) isNewPhenoCreated <- FALSE
    if (length(PR_diff_gid) != 0) {
      phenoData_red <- merge(phenoData, data.frame(Rrow_gid), by.x = "row.names", by.y = "Rrow_gid", sort = FALSE)  
      rownames(phenoData_red) <- phenoData_red[,1]
      phenoData_red <- phenoData_red[,-1]
      #save new pheno file
      #     write.table(phenoData_red, file = paste(outputPath, "/newPhenoData.csv", sep = ""), sep = ",", quote = FALSE, row.names = TRUE, col.names = TRUE)
#       write.csv(phenoData_red, file = paste(outputPath, "/newPhenoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
      isNewPhenoCreated <- TRUE
      phenoData <- phenoData_red
      P_gid <- rownames(phenoData)
    }
    write.csv(phenoData, file = paste(outputPath, "/phenoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
    
    ###reduce (if needed) rel matrix
    #   relData_red <- merge(relData, data.frame(P_gid), by.x = "row.names", by.y = "P_gid", sort = FALSE)
    #   rownames(relData_red) <- relData_red[,1]
    #   relData_red <- relData_red[,-1]
    isNewRelMatCreated <- FALSE
    if (length(RP_diff_gid) != 0) {
      relData_red <- merge(relData, data.frame(P_gid), by.x = "row.names", by.y = "P_gid", sort = FALSE)
      #delete column
      relData_red <- relData_red[,-(match(RP_diff_gid,colnames(relData_red)))]
      rownames(relData_red) <- relData_red[,1]
      relData_red <- relData_red[,-1]
#       write.csv(relData_red, file = paste(outputPath, "/newRelData.csv", sep = ""), quote = FALSE, row.names = TRUE)
      isNewRelMatCreated <- TRUE
      relData <- relData_red
      Rrow_gid <- rownames(relData)
      Rcol_gid <- colnames(relData)
      ###Check that relData is a valid symmetric matrix?
      #     if (sum(Rrow_gid != Rcol_gid) != 0 || (!isSymmetric(relData))) stop("Invalid relationship matrix.")
    }
    write.csv(relData, file = paste(outputPath, "/relData.csv", sep = ""), quote = FALSE, row.names = TRUE)
  }
  
  #match geno and map data
  if (!is.null(geno_file) && !is.null(map_file)) {
    ##genotypes in genoData w/c are not in mapData; display (if any)
    GM_diff_mid <- setdiff(G_mid, M_mid)
    
    ##genotypes in mapData w/c are not in genoData; display (if any)
    MG_diff_mid <- setdiff(M_mid, G_mid)
    
    ###reduce (if needed) genoData, sort markers as in mapData
    #   phenoData_red <- merge(phenoData, data.frame(G_gid), by.x = "row.names", by.y = "G_gid", sort = FALSE)  
    #   rownames(phenoData_red) <- phenoData_red[,1]
    #   phenoData_red <- phenoData_red[,-1]
    #   View(phenoData_red)
    if (is.null(isNewGenoCreated)) isNewGenoCreated <- FALSE
    if (length(GM_diff_mid) != 0) {
      genoDatat <- as.data.frame(t(genoData)) # t(G_dataRed) ###no -1 row?
      genoDatat_red <- merge(genoDatat, data.frame(M_mid), by.x = "row.names", by.y = "M_mid", sort = FALSE)  
      rownames(genoDatat_red) <- genoDatat_red[,1]
      genoDatat_red <- genoDatat_red[,-1]
      genoData_red <- as.data.frame(t(genoDatat_red))
#       write.csv(genoData_red, file = paste(outputPath, "/newGenoData.csv", sep = ""), quote = FALSE, row.names = TRUE)
      isNewGenoCreated <- TRUE
      genoData <- genoData_red
      G_mid <- colnames(genoData)
    }
    write.csv(genoData, file = paste(outputPath, "/genoData.csv", sep = ""), quote = FALSE, row.names = TRUE)

    ###reduce (if needed) mapData
    #   genoData_red <- merge(genoData, data.frame(P_gid), by.x = "row.names", by.y = "P_gid", sort = FALSE)
    #   rownames(genoData_red) <- genoData_red[,1]
    #   genoData_red <- genoData_red[,-1]
    isNewMapCreated <- FALSE
    if (length(MG_diff_mid) != 0) {
      mapData_red <- merge(mapData, data.frame(G_mid), by.x = "row.names", by.y = "G_mid", sort = FALSE)
      rownames(mapData_red) <- mapData_red[,1]
      mapData_red <- mapData_red[,-1]
#       write.csv(mapData_red, file = paste(outputPath, "/newMapData.csv", sep = ""), quote = FALSE, row.names = TRUE)
      isNewMapCreated <- TRUE
      mapData <- mapData_red
      M_mid <- rownames(mapData)
    }
    write.csv(mapData, file = paste(outputPath, "/mapData.csv", sep = ""), quote = FALSE, row.names = TRUE)
  }
  

  #
#   gsObject <- create.gpData(pheno = phenoData, geno = genoData, map = mapData) #, pedData)
#   gsObject$relData <- relData #, pedData)
##   gsData <- create.gpData(data.frame(phenoData[,-c(2:3,6)]), genoData, mapData) #, pedData)

  return(list(
    isNewPhenoCreated = isNewPhenoCreated,
    isNewGenoCreated = isNewGenoCreated,
    isNewRelMatCreated = isNewRelMatCreated,
    isNewMapCreated = isNewMapCreated,
#     isNewPedCreated = isNewPedCreated,
    phenoData = phenoData,
    genoData = genoData,
    mapData = mapData,
    relData = relData
# ,
#     gsObject = gsObject
##     ,
##     pedData = pedData
    ))

}