GEOneStage_resid <- function(ge1sout, exptl.design = c("RCB", "Alpha", "RowCol"), respvar, geno, row, column = NULL, rep = NULL, env,is.genoRandom = FALSE) {
    for (i in (1:length(ge1sout$output))) {
      if (is.genoRandom) {names.resid <- paste(ge1sout$output[[i]]$respvar,"resid_random",sep="_")
      } else names.resid <- paste(ge1sout$output[[i]]$respvar,"resid_fixed",sep="_")
        if (i == 1) {
    	  	residuals <- as.data.frame(ge1sout$output[[i]]$residuals)
  			  names(residuals) <- names.resid
          if (exptl.design == "RCB"){
            residuals <-  cbind(ge1sout$output[[i]]$data[,match(env, names(ge1sout$output[[i]]$data))],ge1sout$output[[i]]$data[,match(geno, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(row, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(respvar[[i]], names(ge1sout$output[[i]]$data))], residuals)
            colnames(residuals) <- c(env, geno, row, respvar[i], names.resid)
          } else if (exptl.design == "Alpha") {
            residuals <-  cbind(ge1sout$output[[i]]$data[,match(env, names(ge1sout$output[[i]]$data))],ge1sout$output[[i]]$data[,match(geno, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(rep, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(row, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(respvar[[i]], names(ge1sout$output[[i]]$data))], residuals)
            colnames(residuals) <- c(env, geno, rep, row, respvar[i], names.resid)
          } else if (exptl.design == "RowCol") {
            residuals <-  cbind(ge1sout$output[[i]]$data[,match(env, names(ge1sout$output[[i]]$data))],ge1sout$output[[i]]$data[,match(geno, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(rep, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(row, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(col, names(ge1sout$output[[i]]$data))], ge1sout$output[[i]]$data[,match(respvar[[i]], names(ge1sout$output[[i]]$data))], residuals)
            colnames(residuals) <- c(env, geno, rep, row, col, respvar[i], names.resid)
          }
		    } else {
			    resid2 <- as.data.frame(ge1sout$output[[i]]$residuals)
  		    names(resid2) <- names.resid
          residuals <- cbind(residuals,resid2)		

        }
    }
    return(residuals)
}
