###########
#  createGSPlots
#' Function for creating graphs for genomic selection
#
#  Parameters:
#' @param outputPath path where output will be saved
#' @param gpData1 R gpData object
#' @param relData R data object of class relationship matrix
#' @param traitNames names of trait variables
# 
###########

createGSPlots <- function(outputPath, gpData1 = NULL, relData, traitNames) {
  
    phenoData <- na.omit(data.frame(gpData1$pheno))

  #1. plot phenotypic trait(s)
  png(paste(outputPath, "phenoPlot.png", sep = ""), res = NA)
  if (ncol(phenoData) == 1) { 
    #input trait
    print(hist(phenoData[,1], main = paste("Histogram of", traitNames[1]), xlab = traitNames[1]))
  } else {
    print(pairs(phenoData, diag.panel = panel.hist, lower.panel = panel.corCoef))
  }
  dev.off()

  #2. plot marker map - column 2- chr and column 3 - pos
  if (!is.null(gpData1$map)) { 
    png(paste(outputPath, "markerPlot.png", sep = ""), res = NA)
    mapData <- gpData1$map 
    print(plotGenMap(mapData[,c(2:3)], TRUE, FALSE, ylab = "pos [cM]"))
    dev.off()
  }
# } else {


  #3. plot LD
  if (!is.null(gpData1$map)) { 
    mygp <- gpData1 
    mygp$info$codeGeno <- TRUE #TODO: check with recode
    nchro <- nlevels(as.factor(mapData$chr))
#  if (!is.null(gpData1$map)) { 

# #     gpdataLDchN <- pairwiseLD(mygp,chr=NULL,type="data.frame") ###***1st chromosome only?
#     #one graph per file
#     for (i in 1:nchro) {
# #     for (i in 1:length(gpdataLDchN)) {
#       gpdataLD <- pairwiseLD(mygp,chr=i,type="data.frame") ###***for all chromosomes
#       png(paste(outputPath, "LDplot_chr", i,".png", sep = ""))
# #       plot(gpdataLDchN[1])
#       plot(gpdataLD)
#       dev.off()    
#     }

    #all graphs in one file
    numr <- nchro/3
    numr = if (nchro%%3 == 0) nchro/3 else round(nchro/3) + 1
    png(paste(outputPath, "LDplot.png", sep = ""), width = 960, height = 480*numr, res = NA)

    par(mfrow = c(numr,3))
    for (i in 1:nchro) {
      gpdataLD <- pairwiseLD(mygp,chr=i,type="data.frame") ###***for all chromosomes
      print(plot(gpdataLD))
    }
    dev.off()
  }

  #4. manhattan plot

  if (!is.null(gpData1$map)) { 
    #   if (!is.null(mapData)) {
    png(paste(outputPath, "manhattanPlot.png", sep = ""))  
    gpDataC <- codeGeno(gpData1)
    rrblup <- gpMod(gpDataC, mod = "BLUP", markerEffects = TRUE)
    print(manhattanPlot(rrblup$markerEffects, gpDataC, ylab = "marker effects"))
    dev.off()  
  }

  #5. heat map of relationship matrix
  #get matrix from file named rel_file if not in gsData1
  if (!is.null(relData)) {
    relMat = relData
    class(relMat) = "relationshipMatrix"
    png(paste(outputPath, "heatMap_relMat.png", sep = ""), res = NA)  
    print(plot(relMat))
    dev.off()
  }

#   #if not from file
#   str(gpData1)
#   class(gpData1$relMat) = "relationshipMatrix"
#   plot(gpData1$relMat)
}