# -------------------------------------------------------------------------------------
# Description: This function perform Stability Analysis using Finlay Method.  
#              Computation based on Statistical and Biometrical Techniques in Plant 
#              Breeding by Jawahar R. Sharma
# R Script Created by: Violeta I Bartolome for International Rice Research Institute
# R Function Created by: Alaine A. Gulles 12.16.2014 for IRRI
# R Function Modified by: Alaine A. Gulles 05.29.2014
# Note: If data file contains vector of Residual Errors, get routine from 
#       agricolate/ammi on how to get average Pooled error
# -------------------------------------------------------------------------------------
# Parameters:
# data - a data frame
# respvar - a character string; name of the response variable
# env - a character string; name of the environment factor
# geno - a character string; name of the genotype factor
# outputPath - NULL or path where graphs will be save
# -------------------------------------------------------------------------------------

FinlayStability <- function(data, respvar, geno, env, outputPath = NULL) UseMethod("FinlayStability")

FinlayStability.default <- function(data, respvar, geno, env, outputPath = NULL) {
     
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
     
     # check if the dataset contain the residual variance and number of rep
     nameIndex <- match(c(paste(respvar, "_sigma2",sep = ""),paste(respvar, "_No.rep", sep = "")),colnames(data))
     
     if (any(is.na(nameIndex))) {
          Pooled.error <- 0
          r <- 1
     } else {
          numrepEnv <- tapply(data[,paste(respvar, "_No.rep", sep = "")], data[,env], mean)
          sigma2Env <- tapply(data[,paste(respvar, "_sigma2", sep = "")], data[,env], mean)
          Pooled.error <- sum(numrepEnv*sigma2Env)/sum(numrepEnv)
          r <- 1/mean(1/numrepEnv)
     }
     
     #### combined analysis over env
     
     model1 <- aov(formula(paste(respvar, "~", env, "+", geno)), data = data)
     out <- anova(model1)
     
     ### site.means
     Site.means <- data.frame(levels(data[, match(env, names(data))]), as.data.frame(tapply(data[, match(respvar, names(data))], data[, match(env, names(data))], mean, na.rm = TRUE)))
     colnames(Site.means) <- c(env, "Site.Index")
     
     ### genotype means
     Geno.means <- data.frame(levels(data[, match(geno, names(data))]), as.data.frame(tapply(data[, match(respvar, names(data))], data[, match(geno, names(data))], mean, na.rm = TRUE)))
     colnames(Geno.means) <- c(geno, "Means")
     
     ###### G x E means in parallel format 
     GE.means.wide <- reshape(data[,c(respvar, env, geno)], v.names = respvar, idvar = env, timevar = geno, direction = "wide")
     colnames(GE.means.wide) <- gsub(paste(respvar, ".", sep = ""),"",colnames(GE.means.wide))
     
     data.all <- merge(Site.means, GE.means.wide, by = env)
     
     ###### sum of squares due to regression
     temp <- as.matrix(c(1:(e*g)), nrow = e, ncol = g)
     dim(temp) <- c(e,g)
     for (i in 1:nlevels(data[,geno]))  { temp[,i] <- data.all[,i+2]/data.all[,2] }
     columnsum2 <- (colSums(temp))^2
     

     ### ANOVA Table
     newAOVTable <- rbind(out, out)
     rownames(newAOVTable)[3:6] <- c(paste(geno, "x", env), "Regression", "Dev.fr. Reg", "Error")
     newAOVTable[4:6,"Df"]<- c(g -1, (e - 1)*(g - 1) - (g - 1), e * (r - 1) * (g - 1))
     
     newAOVTable[1:4,"Sum Sq"] <- c(r * newAOVTable[1:3,"Sum Sq"],r^2*sum(columnsum2))
     newAOVTable[5:6,"Sum Sq"] <- c(newAOVTable[3,"Sum Sq"] - newAOVTable[4,"Sum Sq"], Pooled.error*newAOVTable[6,"Df"])
     
     newAOVTable[1:3,"Mean Sq"] <- c(r * newAOVTable[1:3,"Mean Sq"])
     newAOVTable[4:6,"Mean Sq"] <- c(newAOVTable[4:5,"Sum Sq"]/newAOVTable[4:5,"Df"], Pooled.error)
     
     newAOVTable[2:5,"F value"] <- newAOVTable[2:5,"Mean Sq"]/Pooled.error
     newAOVTable[1,c("F value", "Pr(>F)")]   <- NA
     newAOVTable[2:5,"Pr(>F)"] <- 1 - pf(newAOVTable[2:5,"F value"], newAOVTable[2:5,"Df"], newAOVTable[6,"Df"])
     
     cat("Analysis of Variance\n\n")
     printAOVTable(newAOVTable)
     cat("\n\n")
     
     #### regression for each Genotype with site index 
     slope <- as.matrix(c(1:(6*g)),nrow=g,ncol=6)
     dim(slope) <- c(g,6)
          
     ###### check whether variation among site means is not equal to zero  
     
     if (var(data.all$Site.Index) != 0)  {
          for(i in 1:g)   {
               #### check number of observations if < 3 give error message
               if (length(na.omit(data.all[[i+2]])) < 3)   {
                    slope[i,] <- c(NA,NA,NA,NA, NA, NA)
               } else {  
                    model <- lm(data.all[[i+2]]~data.all[[2]],data.all)
                    temp <- summary(model)
                    slope[i,] <- c(temp$coef[2,], anova(model)[1, "Mean Sq"], anova(model)[2, "Mean Sq"]) 
               }   
          }  
          rownames(slope) <- levels(factor(data$Geno))
          colnames(slope) <- c("Slope", "SE", "Tvalue", "Prob", "MSReg", "MSDev")
          #title <- "Stability analysis using Finlay-Wilkinson model"      
          #print(title, quote=F)
          #print(slope, quote=F)
          slopeNew <- data.frame(rownames(slope), slope)
          colnames(slopeNew)[1] <- " "
          
          printDataFrame(slopeNew, digits = 4)
          
          slope <- as.data.frame(slope)
          slope$Tvalue <- (slope$Slope-1)/slope$SE
          slope$Prob <- 2*pt(-abs(slope$Tvalue),e-2)
          slope[,geno] <- rownames(slope)
          slope$MSReg <- slope$MSDev <- rownames(slope)<- NULL
          
          slope <- merge (Geno.means, slope, by=geno)
          #### Note:  t-test is to test Ho: slope=1
          
          #### plot of predicted lines
          
          data.all <- data.all[order(data.all$Site.Index),]
          
          minY <- min(data.all[,c(3:(g+2))]) 
          maxY <- max(data.all[,c(3:(g+2))]) 
          
          if (is.null(outputPath)) { outputPath <- getwd() }
          
          png(filename = paste(outputPath, "/",respvar, "_SA Finlay_predictedLines.png", sep = ""))
          plot(data.all[[3]]~Site.Index,data.all, pch=" ", ylim=c(minY,maxY), xlab="Site Index", ylab=respvar)
          
          for (i in 1:g)  {
               model <- lm(data.all[[i+2]]~Site.Index,data.all)
               abline(model, col=i) 
               text(data.all$Site.Index[[e]], predict(model)[e], labels=levels(data[,geno])[i],col=i)
          }
          abline(0,1, col=(i+1), lty=2)
          text(data.all$Site.Index[[e]],data.all$Site.Index[[e]], labels="mean")
          dev.off()
          
          #### plot of means vs slope of each genotype
          
          png(filename = paste(outputPath,"/", respvar, "_SA Finlay_means.png", sep = ""))
          plot(slope$Mean, slope$Slope, pch=' ', xlab="Mean Yield", ylab="Slope")  
          text(slope$Mean, slope$Slope, labels=slope[,geno],col=1)
          
          abline(h=1,v=mean(slope$Mean),lty=1,col="Black")
          dev.off()
          
     }  else  { 
          print ("no variation in site means")  
          slopeNew <- NULL
     }
     return(list(aovtable = newAOVTable, stability = slopeNew, method = "Stabiliyt Analysis using Finlay-Wilkinson Model"))
}
