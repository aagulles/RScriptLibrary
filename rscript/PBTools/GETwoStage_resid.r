GETwoStage_resid <- function(ge2sout,is.genoRandom = FALSE,geno,env,is.single.out = TRUE,respvar,stderr=NULL,sigma2,numrep) {
    for (i in (1:length(ge2sout))) {
      if (is.genoRandom) {names.resid <- paste(ge2sout[[i]]$respvar,"resid_random",sep="_")
  	  } else names.resid <- paste(ge2sout[[i]]$respvar,"resid_fixed",sep="_")
        
	  #to create data frame for residuals
  		if (is.single.out) {
        if (i==1) {
			residuals <- as.data.frame(ge2sout[[i]]$residuals)
  		  	names(residuals) <- names.resid

          	residuals <-  cbind(ge2sout[[i]]$data[,match(env, names(ge2sout[[i]]$data))],ge2sout[[i]]$data[,match(geno, names(ge2sout[[i]]$data))], residuals)
          	colnames(residuals) <- c(env, geno, names.resid)

		} else {
		    resid2 <- as.data.frame(ge2sout[[i]]$residuals)
  		    names(resid2) <- names.resid
	        resid2 <- cbind(ge2sout[[i]]$data[,match(env, names(ge2sout[[i]]$data))],ge2sout[[i]]$data[,match(geno, names(ge2sout[[i]]$data))], resid2)
    	    colnames(resid2) <- c(env, geno, names.resid)

			byVariables <- c(env,geno)
		  	residuals <- merge(residuals,resid2,by= byVariables,sort = TRUE)
        }
  		} else {
        if (is.genoRandom) {
          if (i == 1) {
            residuals <- as.data.frame(ge2sout[[i]]$residuals)
      	    names(residuals) <- names.resid
            residuals <- cbind(ge2sout[[i]]$data[,match(env, names(ge2sout[[i]]$data))],ge2sout[[i]]$data[,match(geno, names(ge2sout[[i]]$data))], residuals)
            colnames(residuals) <- c(env, geno, names.resid)
          } else {
            resid2 <- as.data.frame(ge2sout[[i]]$residuals)
            names(resid2) <- names.resid
            residuals <- cbind(residuals,resid2)
          }
        }
  		}

    }
    return(residuals)
}
