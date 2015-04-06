# --------------------------------------------------------
# GraphHist: creates a histogram of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# nVar - name of the numeric variable
# mTitle - main title for the boxplot, if any
# yAxisLab - label to be used for the y-axis
# xAxisLab - label to be used for the x-axis 
# yMinValue - minimum value for the the y-axis
# yMaxValue - maximum value for the y-axis
# xMinValue - minimum value for the the x-axis
# xMaxValue - maximum value for the x-axis
# useFreq - logical; whether to use frequencies in the y-axis or not (probabilities)
# numBins - number of bins; by default, it is determined using Sturges' formula
# byVar - name of categorical variable for which separate graphs are created
# barColors - vector of rgb values for the colors of the bars
# dispCurve - logical; whether the kernel density estimates are displayed (in a curve) or not
# lineTypes - vector of numbers for types of lines used
# lineColors - vector of rgb values for the colors of the lines, if any
# lineWidths - vector of widths for the lines, if any
# showLeg - logical; whether the graph legend is displayed or not
# legPos - position of the legend (if showLeg is TRUE)
# legTitle - text for the title of the legend, if displayed
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphHist <- function(data, outputPath, nVar, mTitle = NULL, yAxisLab = NULL, xAxisLab = NULL, 
              yMinValue = NULL, yMaxValue = NULL, xMinValue = NULL, xMaxValue = NULL, 
              useFreq = TRUE, numBins = "Sturges", byVar = NULL, 
              barColors, dispCurve = FALSE, lineTypes, lineColors, lineWidths,
              showLeg = FALSE, legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", 
              "top", "topright", "right",  "center"), legTitle = NULL, 
              boxed = TRUE, multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
              orientGraphs = c("left-right", "top-bottom")) { 
  
  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #converts to factor the grouping variable(s)
  if (!is.null(byVar)) { data[,byVar] <- factor(data[,byVar]) }  
  
  if (!multGraphs) {
    numRowsGraphs = 1
    numColsGraphs = 1
  } 
  
  #determines number of cells allocated for graphs (esp. for multiple graphs)
  numCells = numRowsGraphs * numColsGraphs
  
  #determines number of graphs to be created
  if (!is.null(byVar)) {
    numGraphs = (nlevels(data[,byVar]))
  } else numGraphs = 1
  
  graphNum = 1
  #counts the number of files to save
  k = 1
  
  for (m in 1:numGraphs) {
    #creates device for saving graph(s)
    if (!multGraphs) {
      png(filename = paste(outputPath,"histogram", k,".png",sep=""))
      par(mfrow=c(numRowsGraphs, numColsGraphs))
        
    } else {
      if (graphNum == 1 || graphNum %% numCells == 1)	{
        png(filename = paste(outputPath,"histogram", k,".png",sep=""))
         
        if (orientGraphs == "top-bottom") {
          par(mfcol=c(numRowsGraphs, numColsGraphs))
        } else if (orientGraphs == "left-right") {
          par(mfrow=c(numRowsGraphs, numColsGraphs))
        }
      }
    }
    
    #creates data by subgroup, if a grouping variable is defined
    if (!is.null(byVar)) {
      tempData = data[which(data[,byVar] == levels(data[,byVar])[m]),]
      subTitle = paste(byVar,"=",levels(tempData[,byVar])[m], sep=" ")
    } else {
      tempData = data
      subTitle = NULL
    }

    #determines lower and upper limits for the x-axis
    xMinLim = if(!is.null(xMinValue)) { xMinValue
    } else min(0, min(tempData[,nVar])) 
    xMaxLim = if(!is.null(xMaxValue)) { xMaxValue
    } else max(tempData[,nVar]) + 0.25*max(tempData[,nVar]) 
      
    xVarLim = c(xMinLim,xMaxLim)
    
    histData = NULL  
    minCount = NULL
    maxCount = NULL
    minDensity = NULL
    maxDensity = NULL
    
    for (i in 1:length(nVar)) {
      histData[[i]] =  hist(tempData[,nVar[i]], breaks = numBins, plot = F)  
      minCount = min(0, minCount, min(histData[[i]]$counts)) 
      maxCount = max(maxCount, histData[[i]]$counts)
      minDensity = min(0, minDensity, min(histData[[i]]$density)) 
      maxDensity = max(maxDensity, histData[[i]]$density)
    }
    
    #determines lower and upper limits for the y-axis
    if (useFreq) {
      yMinLim = if(!is.null(yMinValue)) { yMinValue
      } else minCount
      yMaxLim = if(!is.null(yMaxValue)) { yMaxValue
      } else {
         maxCount
      }
    } else {
      yMinLim = if(!is.null(yMinValue)) { yMinValue
      } else minDensity
      yMaxLim = if(!is.null(yMaxValue)) { yMaxValue
      } else {
        if (dispCurve) maxDensity + maxDensity * 0.25
        else maxDensity
      }
    }
    yVarLim = c(yMinLim,yMaxLim)
    
    plot(histData[[1]], xlim = xVarLim, ylim = yVarLim, col = barColors[1], freq = useFreq)
    
    if (dispCurve && !useFreq)
      lines(density(tempData[,nVar[1]]), col = lineColors[1], lty = lineTypes[1], lwd = lineWidths[1])
    
    if (length(nVar) > 1) {
      for (i in 2:length(nVar)) {
        plot(histData[[i]], add = T, col = barColors[i], freq = useFreq)
        if (dispCurve && !useFreq)
          lines(density(tempData[,nVar[i]]), col = lineColors[i], lty = lineTypes[i], lwd = lineWidths[i])
      }
    }
        
    legLab = NULL
    if (showLeg) legLab = nVar
    
    if (showLeg) 
      legend(legPos, legend = legLab, title = legTitle, fill = barColors)
    #legend(legPos, legend = legLab, title = legTitle, col = barColors, pch = 15)
    
    #draws a box around the plot
    if (boxed) box(lty = 1)
    
    #increments file number
    if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
      dev.off()
      k = k + 1
    }
    graphNum = graphNum + 1
  }
	
}