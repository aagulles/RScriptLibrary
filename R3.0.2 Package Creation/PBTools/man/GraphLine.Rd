\name{GraphLine}
\alias{GraphLine}
\alias{GraphLine.default}
\title{GraphLine}
\description{
constructs line graph
}
\usage{
GraphLine(data, outputPath, yVars, xVar = NULL, lineVars = NULL, mTitle = NULL, yAxisLab = NULL, 
                      xAxisLab = NULL, yMinValue = NULL, yMaxValue = NULL, axisLabelStyle = 1,
                      byVar = NULL, errBars = FALSE, typeErrBar = c(NULL, "stdErr", "confInt"), errMult = 1, confLev = NULL, 
                      plotCol = NULL, showLineLabels = FALSE, showLeg = FALSE, 
                      legPos = c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center"), 
                      legTitle = NULL, legCol= 1, boxed = TRUE, linePtTypes, lineTypes,  lineWidths, pointChars, pointCharSizes, 
                      multGraphs = FALSE, numRowsGraphs = 1, numColsGraphs = 1, orientGraphs = c("left-right", "top-bottom"))
}
\arguments{
  \item{data}{name of R dataframe}
  \item{outputPath}{folder where graph(s) will be saved}
  \item{yVars}{names of the variable(s) to be plotted against the factor variable}
  \item{xVar}{name of the factor variable}
  \item{lineVars}{name of grouping variable(s) for creating different lines}
  \item{mTitle}{main title for line graph, if any}
  \item{yAxisLab}{labels to be used for the y-axis variables}
  \item{xAxisLab}{label to be used for the x-axis}
  \item{yMinValue}{minimum value(s) for the the y-axis variables}
  \item{yMaxValue}{maximum value for the y-axis}
  \item{axisLabelStyle}{style for the axes labels}
  \item{byVar}{name of categorical variable for which separate graphs are created}
  \item{errBars}{logical; whether error bars are displayed or not}
  \item{typeErrBar}{whether error bars displayed (if any) are based on: standard error or confidence interval}
  \item{errMult}{multiplier to be used on standard error}
  \item{confLev}{confidence level, if confidence interval type error bars are to be displayed}
  \item{plotCol}{vector of rgb values for the color(s) of the lines (or names of colors)}
  \item{showLineLabels}{logical; whether the labels beside the lines are displayed or not}
  \item{showLeg}{logical; whether the graph legend is displayed or not}
  \item{legPos}{position of the legend (if showLeg is TRUE)}
  \item{legTitle}{text for the title of the legend, if displayed}
  \item{legCol}{number of columns used for legend text}
  \item{boxed}{logical; whether a box is drawn around the plot or not}
  \item{linePtTypes}{a}
  \item{lineTypes}{a}
  \item{lineWidths}{a}
  \item{pointChars}{a}
  \item{pointCharSizes}{a}
  \item{multGraphs}{logical; whether multiple graphs will be displayed in a page or not}
  \item{numRowsGraphs}{number of rows of graphs to allow to be displayed}
  \item{numColsGraphs}{number of columns of graphs to allow to be displayed}
  \item{orientGraphs}{whether multiple graphs are to be displayed from left-to-right or top-to-bottom}
}
\author{
Rose Imee Zhella A. Morantte, Nellwyn L. Sales
}
\keyword{linegraph}
