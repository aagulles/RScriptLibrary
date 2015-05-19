# ------------------------------------------------------------------
# graphSEAdiagplot
# Description: diagnostic plots for single site analysis
# Created by: Alaine A. Gulles for IRRI
# ------------------------------------------------------------------

graphSEAdiagplot <- function(residualValues, fittedValues, respvar, is.random = FALSE, env = NULL, envLevel = NULL) UseMethod("graphSEAdiagplot")

graphSEAdiagplot.default <- function(residualValues, fittedValues, respvar, is.random = FALSE, env = NULL, envLevel = NULL) {
     par(mfrow = c(2,2))
     # scatterplot of the residuals
     plot(residualValues ~ fittedValues, xlab = "Predicted Values", pch = 15, cex = 0.7, col = rgb(0.3, 0.3, 0.7, 0.2),
          ylab = "Residuals", main = "Scatterplot of Residuals \nagainst Fitted Values")
     # qqplot of residuals
     qqnorm(residualValues)
     qqline(residualValues, col = 2, main = title, sub = respvar)
     # freq dist of residuals
     hist(residualValues, main = "Histogram of Residuals", 
          col = rgb(0.3, 0.3, 0.7, 0.2), xlab = "Residual", 
          ylab = "Frequency")
     # create a blank plot
     plot(seq(1:10) ~ seq(1:10), type = "n", axes = FALSE, xlab = "", ylab = "")
     noteString <- paste("Note: Residuals plotted are taken from the model where genotype is ", sep = "")
     if (is.random) { noteString <- paste(noteString, "random ", sep = "") 
     } else { noteString <- paste(noteString, "fixed ", sep = "") }
     noteString <- paste(noteString, "and response variable = ", respvar, sep = "")
     if (!is.null(env)) { noteString <- paste(noteString, " and ", env, " = ", envLevel, sep = "") }
     noteString <- paste(noteString, ".", sep = "")
     text(5,7, paste(strwrap(noteString, width = 30), sep = "", collapse = "\n"))
}

