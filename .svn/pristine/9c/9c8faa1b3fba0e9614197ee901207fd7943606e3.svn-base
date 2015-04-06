# --------------------------------------------------------
# ARGUMENTS:
# data - a dataframe
# respvar - a string; variable name of the response variable
# env - a string; variable name of the environment variable
# is.random - logical; indicates whether genotype/treatment is random or not; default value is FALSE (FIXED factor)
# result - output of single environment analysis 
# --------------------------------------------------------

graph.mea2s.diagplots <- function(data, respvar, is.random = FALSE, result) {
	
  #dir.create("plots")
  
  #create diag plots
  for (i in (1:length(respvar))) {
    xlabel = respvar[i]
	cc = rgb(0.3, 0.3, 0.7, 0.2)

      if (is.random) {
        filename1 = paste(getwd(),"/diagPlotsMea2s_",respvar[i], "_random.png",sep="");
      } else filename1 = paste(getwd(),"/diagPlotsMea2s_",respvar[i], "_fixed.png",sep="");
      
	    png(filename = filename1); par(mfrow = c(2,2)) #par(mfrow = n2mfrow(length(respvar)*nlevels(data[,match(env, names(data))])));

      #scatterplot of residuals against fitted values
  		plot(result[[i]]$residuals~result[[i]]$fitted.values, xlab = "Predicted Values", pch = 15, cex = 0.7, col = cc, ylab = "Residuals", main = xlabel);

      #qqplot of residuals
		qqnorm(result[[i]]$residuals); qqline(result[[i]]$residuals, col=2, main=title,sub=xlabel)
	
	  #freq dist of residuals
		hist(result[[i]]$residuals, main = "Residuals",  col = cc, xlab = "Residual", ylab = 'Frequency' )

      dev.off();
	}	
}