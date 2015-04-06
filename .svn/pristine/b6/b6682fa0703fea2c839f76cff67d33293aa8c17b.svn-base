# --------------------------------------------------------
# GraphBoxplot: creates a boxplot of the data
# ARGUMENTS:
# data - name of R dataframe
# outputPath - folder where graph(s) will be saved
# nVar - name of the numeric variable
# cVar - name of the categorical variable
# mTitle - main title for the boxplot, if any
# nAxisLab - label to be used for the axis of the numeric variable
# cAxisLab - label to be used for axis of the categorical variable
# yMinValue - minimum value for the the y-axis
# yMaxValue - maximum value for the y-axis
# plotHoriz - logical; whether plots are drawn horizontally or not
# byVar - name of categorical variable for which separate graphs are created
# boxColor - vector of rgb values for the color of the box
# legPos - position of the legend (if showLeg is TRUE)
# legTitle - text for the title of the legend, if displayed
# boxed - logical; whether a box is drawn around the plot or not
# multGraphs - logical; whether multiple graphs will be displayed in a page or not
# numRowsGraphs - number of rows of graphs to allow to be displayed
# numColsGraphs - number of columns of graphs to allow to be displayed
# orientGraphs - whether multiple graphs are to be displayed from left-to-right or top-to-bottom
# --------------------------------------------------------

GraphBoxplot <- function(data, outputPath, nVar, cVar = NULL, mTitle = NULL, nAxisLab = NULL, cAxisLab = NULL,
                yMinValue = NULL, yMaxValue = NULL, plotHoriz = FALSE, byVar = NULL, 
                boxSize = 0.8, boxColor = NULL, boxFillColor = NULL,
                medLineType = 1, medLineWidth = 3, medColor = NULL, whiskLineType = 2, whiskLineWidth = 1,
                whiskColor = NULL, outChar = 1, outCharSize = 1, outColor = NULL,
                boxed = TRUE, multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, 
                orientGraphs = c("left-right", "top-bottom")) { 

  #reads data
  if (is.character(data)) { data <- eval(parse(text = data)) }
  
  #creates a grouping variable with ones for all rows if no grouping variable is declared
  if (is.null(cVar)) { 
    grp.Var = rep(1,nrow(data))
    cVar = "grp.Var"
    data = cbind(data,grp.Var)
  }  
  
  #converts to factor the grouping variable(s)
  data[,cVar] <- factor(data[,cVar]) 
  if (!is.null(byVar)) { data[,byVar] <- factor(data[,byVar]) }  
  
  if (!multGraphs) {
    numRowsGraphs = 1
    numColsGraphs = 1
  } 
  
  #determines number of cells allocated for graphs (esp. for multiple graphs)
  numCells = numRowsGraphs * numColsGraphs
  
  numGroups = 1
  #determines number of graphs to be created
  if (!is.null(byVar)) {
    numGroups = nlevels(data[,byVar])
    numGraphs = nlevels(data[,byVar]) * length(nVar)
  } else numGraphs = length(nVar)
  
  graphNum = 1
  #counts the number of files to save
  k = 1
  
  for (m in 1:numGroups) {
    #creates data by subgroup, if a grouping variable is defined
    if (!is.null(byVar)) {
      tempData = data[which(data[,byVar] == levels(data[,byVar])[m]),]
      subTitle = paste(byVar,"=",levels(tempData[,byVar])[m], sep=" ")
    } else {
      tempData = data
      subTitle = NULL
    }

    tempData = tempData[order(tempData[,cVar]),]
     
    namesX = levels(tempData[,cVar])
 
    for (i in 1:length(nVar)) {
      #creates device for saving graph(s)
      if (!multGraphs) {
        png(filename = paste(outputPath,"boxplot",k,".png",sep=""))
        par(mfrow=c(numRowsGraphs, numColsGraphs))
        
      } else {
        if (graphNum == 1 || graphNum %% numCells == 1)  {
          png(filename = paste(outputPath,"boxplot",k,".png",sep=""))
          
          if (orientGraphs == "top-bottom") {
            par(mfcol=c(numRowsGraphs, numColsGraphs))
          } else if (orientGraphs == "left-right") {
            par(mfrow=c(numRowsGraphs, numColsGraphs))
          }
        }
      }
      
      #determines lower and upper limits for the y-axis
      yMinLim = if(!is.na(yMinValue[i])) { yMinValue[i]
      } else min(0, min(tempData[,nVar[i]])) 
      
      yMaxLim = if(!is.na(yMaxValue[i])) { yMaxValue[i]
      } else max(tempData[,nVar[i]]) 
      
      yVarLim = c(yMinLim,yMaxLim)
      
      if (!is.null(cVar)) eqn = "tempData[,nVar[i]] ~ tempData[,cVar]" else eqn = "tempData[,nVar[i]]"
      eqn <- eval(parse(text = eqn))
      
      if (plotHoriz) {
        yAxisLab = cAxisLab
        xAxisLab = nAxisLab[i]
      } else {
        xAxisLab = cAxisLab
        yAxisLab = nAxisLab[i]
      }
      boxplot(eqn, ylab = yAxisLab, xlab = xAxisLab, main = mTitle, 
              sub = subTitle, ylim = yVarLim, las = 1,
              frame.plot = boxed, horizontal = plotHoriz, names.arg = namesX, 
              pars = list(boxwex = boxSize, boxcol = boxColor, 
              boxfill = boxFillColor, medlty = medLineType, medlwd = medLineWidth, medcol = medColor,
              whisklty = whiskLineType, whisklwd = whiskLineWidth, whiskcol = whiskColor, staplecol = whiskColor,
              outpch = outChar, outcex = outCharSize, outcol = outColor))

      #increments file number
      if ((graphNum %% numCells == 0) || graphNum == numGraphs) {
        dev.off()
        k = k + 1
      }
      graphNum = graphNum + 1
      
    }
    
  }
	
}