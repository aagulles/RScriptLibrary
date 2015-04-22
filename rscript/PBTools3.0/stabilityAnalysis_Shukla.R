# --------------------------------------------------------------------------------------------
# Description: This function perform Stability Analysis using Shukla's minimum variance model.  
# R Script Created by: Violeta I Bartolome for International Rice Research Institute
# R Function Created by: Alaine A. Gulles 02.27.2015 for IRRI
# --------------------------------------------------------------------------------------------
# Parameters:
# data - a data frame
# respvar - a character string; name of the response variable
# env - a character string; name of the environment factor
# geno - a character string; name of the genotype factor
# outputPath - NULL or path where graphs will be save
# -------------------------------------------------------------------------------------

ShuklaStability <- function(data, respvar, geno, env, outputPath = NULL) UseMethod("ShuklaStability")

ShuklaStability.default <- function(data, respvar, geno, env, outputPath = NULL) {
     
     if (is.na(match(respvar, names(data))) ||  
              is.na(match(geno, names(data))) || 
              is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }
     
     if (length(respvar) > 1) { respvar <- respvar[1] }
     if (length(geno) > 1) { geno <- geno[1] }
     if (length(env) > 1) { env <- env[1] }
     
     data[,match(geno, names(data))] <- factor(data[,match(geno, names(data))])
     data[,match(env, names(data))] <- factor(data[,match(env, names(data))])
     
     g <- nlevels(data[,geno])
     e <- nlevels(data[,env])

     if (!is.na(match("package:lme4", search()))) { 
          is.pkgload <- TRUE
          detach("package:lme4")
     } else is.pkgload <- FALSE
     
     library(nlme)
     
     data <- subset(data, subset = (is.na(data[,match(respvar, names(data))]) == F))
     errormsg <- NULL
     #myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + ", env, sep = "")
     #myformula2 <- paste("varIdent(form = ~ 1 |", geno,")", sep = "")
     command <- paste("gls(formula(", respvar," ~ 1 + ", geno," + ", env,"), weights = varIdent(form = ~ 1 |", geno,"), data = data)", sep = "") 
     if (class(try(parse(text = command), TRUE)) == "try-error") {
          errormsg <- geterrmessage()
     } else {
          model1 <- eval(parse(text = command))
          if (class(try(parse(text = "intervals(model1)"), TRUE)) == "try-error"){ 
               errormsg <- geterrmessage()  
          } else {
               int <- intervals(model1)
               par <- as.data.frame(int$varStruct)
               par <- rbind(c(1,1,1), par)
               rownames(par) <- levels(data[,match(geno, names(data))])
               sigma <- as.data.frame(int$sigma)
               par$lower <- par$lower*sigma[rownames(sigma) == "lower",]
               par$est.   <- par$est.*sigma[rownames(sigma) == "est.",]
               par$upper <- par$upper*sigma[rownames(sigma) == "upper",]
               par[,geno] <- rownames(par)
               rownames(par) <- NULL
               stability <- par[,c(geno, "lower", "est.", "upper")]
               
               respvarMean <- aggregate(data[,respvar] ~ data[,geno], data = data, mean)
               colnames(respvarMean) <- c(geno, respvar)
               dataAll <- merge(par, respvarMean, by = geno)
               
               ymin <- min(dataAll$est.)
               ymax <- max(dataAll$est.)
               
               xmin <- min(dataAll[,respvar])
               xmax <- max(dataAll[,respvar])
               
               if (is.null(outputPath)) { outputPath <- getwd() }
               
               png(filename = paste(outputPath, "/",respvar, "_SA Shukla.png", sep = ""))
               par(cex = 0.8)
               plot(dataAll[,respvar], dataAll$est., cex = 0.8, pch = " ", ylim = c(ymin,ymax), 
                    xlim = c(xmin,xmax), main = "Stability Analysis based on Shukla's model",
                    frame = TRUE, xlab = paste("Mean", respvar), ylab = "Variance")
               
               text(dataAll[,respvar], dataAll$est., labels = dataAll[, geno], col = 1)
               abline(h = mean(dataAll$est.), v = mean(dataAll[,respvar]), lty = 2, col = "red")
               dev.off()
          }
               
     }
     
     detach("package:nlme")
     #if (is.pkgload) library(lme4)     
     return(list(stability = stability, errormsg = errormsg, method = "Stability Analysis using Shukla's model"))
}
