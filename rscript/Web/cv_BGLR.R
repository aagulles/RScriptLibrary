# cross validation for Bayessian methods in BGLR package

# input
## geno: recoded SNP markers (BGLR format)
## pheno: phenotype file (for one trait and one replicate for each trait)
## R2: trait heritability
## model: Bayessian methods used ("BayesA","BayesB","BayesC",or "RKHS")
## Seed: seed used for randomization whil sampling
## k: number of fold for CV
## Rep: times of repetition sampling for CV
## sampling: sampling method used ("random" or "within popStruc")
## popStruc: population or family information used (if sampling=="within popStruc")
## nIter 
## burnIn 



### genotype data
# rows: individual
# rownames: individual ID
# columns: marker
# columnnames: marker ID

### phenotype data
# rows: individual
# rownames: individual ID  

### NOTE: 
# individual IDs of phenotype file should be identical with that of genotype file 
# NAs are not permitted in both phenotype and genotype data

# output
# a cvData object (defined in synbreed package)
library(synbreed)
library(BGLR)
library(randomForest)
cv_BGLR<-function(geno,pheno,R2,model=c("BayesA","BayesB","BayesC","RKHS","randomForest"),Seed=NULL,k=5,Rep=1,sampling=c("random","within popStruc"),popStruc_file = NULL,nIter,burnIn,ntree=500){
  model=match.arg(model)
  sampling=match.arg(sampling)

  DF=5
  probIn=0.5
  counts=10
  shape=1.01

  MSx=sum(apply(FUN=var,MARGIN=2,X=geno))

  mtrys<-floor(sqrt(ncol(geno)))                 # Numbers of candidate splitting variables at each node for random Forest

  pheno=as.matrix(na.omit(pheno))
  colnames(pheno) <- "TRAIT"
  n=nrow(pheno)
  if (k < 2) 
    stop("folds should be equal or greater than 2")
  if (k > n) 
    stop("folds should be equal or less than the number of observations")
  if (sampling=="within popStruc" ) {
    popStruc=as.matrix(na.omit(read.csv(popStruc_file,sep=',',row.names=1,stringsAsFactors=F)))
    if(!sum(rownames(pheno)%in%rownames(popStruc))==nrow(pheno))
      stop("Some phenotyped individuals had not population structure information")
    popStruc=as.matrix(popStruc[rownames(popStruc)%in%rownames(pheno),])
    popStruc=as.matrix(popStruc[order(rownames(popStruc)),1])
#     identical(rownames(popStruc),rownames(pheno)) # check
  }
    
  if (!is.null(Seed)) 
    set.seed(Seed)
  
  seed2 <- round(runif(Rep, 1, 1e+05), 0)

  COR2 <- matrix(NA, nrow = k, ncol = Rep)
  colnames(COR2) <- paste("rep", 1:Rep, sep = "")
  rownames(COR2) <- paste("fold", 1:k, sep = "")

  rCOR2 <- matrix(NA, nrow = k, ncol = Rep)
  colnames(rCOR2) <- paste("rep", 1:Rep, sep = "")
  rownames(rCOR2) <- paste("fold", 1:k, sep = "")

  lm.coeff <- matrix(NA, nrow = k, ncol = Rep)
  colnames(lm.coeff) <- paste("rep", 1:Rep, sep = "")
  rownames(lm.coeff) <- paste("fold", 1:k, sep = "")

  y.TS2 <- NULL

  n.TS <- matrix(NA, nrow = k, ncol = Rep)
  colnames(n.TS) <- paste("rep", 1:Rep, sep = "")
  rownames(n.TS) <- paste("fold", 1:k, sep = "")

  n.DS <- matrix(NA, nrow = k, ncol = Rep)
  colnames(n.DS) <- paste("rep", 1:Rep, sep = "")
  rownames(n.DS) <- paste("fold", 1:k, sep = "")

  id.TS2 <- list()

  m10 <- matrix(NA, nrow = Rep, ncol = 1)
  colnames(m10) <- "m10"
  rownames(m10) <- paste("rep", 1:Rep, sep = "")
  
  mse <- matrix(NA, nrow = k, ncol = Rep)
  colnames(mse) <- paste("rep", 1:Rep, sep = "")
  rownames(mse) <- paste("fold", 1:k, sep = "")
  
  D=as.matrix(dist(geno,method="euclidean"))
  h=quantile(as.vector(D)^2,probs=.05)
  K1=exp(-5/h*(as.matrix(D)^2))
  K2=exp(-1/h*(as.matrix(D)^2))
  K3=exp(-1/5/h*(as.matrix(D)^2))
  
  for (i in 1:Rep) {
    #i=1
    y.u <- unique(rownames(pheno))
    if (sampling == "random") {
      set.seed(seed2[i])
      modu <- n%%k
      val.samp2 <- sample(c(rep(1:k, each = (n - modu)/k), 
                            sample(1:k, modu)), n, replace = FALSE)
      val.samp3 <- data.frame(y.u, val.samp2)
    }
    if (sampling == "within popStruc") {
      which.pop <- unique(popStruc)
      val.samp3 <- NULL
      for (j in 1:length(which.pop)) {
        y2 <- matrix(y.u[popStruc == which.pop[j]], ncol = 1)
        set.seed(seed2[i] + j)
        modu <- nrow(y2)%%k
        if (!modu == 0) 
          val.samp <- sample(c(rep(1:k, each = (nrow(y2) - modu)/k), sample(1:k, modu)), nrow(y2), replace = FALSE)
        if (modu == 0) 
          val.samp <- sample(rep(1:k, each = (nrow(y2))/k), 
                             nrow(y2), replace = FALSE)
        val.samp2 <- data.frame(y2, val.samp)
        val.samp3 <- as.data.frame(rbind(val.samp3, val.samp2))
      }
      val.samp3 <- val.samp3[order(as.character(val.samp3[, 1])), ]
    }
    
    y.TS <- NULL
    id.TS <- list()
    
    for (ii in 1:k) {
      #ii=1
      samp.es <- val.samp3[!(val.samp3[, 2] %in% ii), ]
      samp.ts <- val.samp3[!is.na(val.samp3[, 2]), ]
      samp.ts <- samp.ts[samp.ts[, 2] == ii, ]
      namesDS <- c(samp.es[, 1], samp.ts[, 1])
      pheno.samp <- pheno
      samp.kf <- val.samp3[, 2] == ii
      samp.kf[is.na(samp.kf)] <- FALSE
      pheno.samp[samp.kf == TRUE, "TRAIT"] <- NA
      
      y2 <- as.matrix(pheno[rownames(pheno) %in% samp.ts[, 1], "TRAIT"])
      X2 <- geno[(rownames(geno) %in% samp.ts[, 1]), ]
      if(model %in% c("BayesA","BayesB","BayesC")){
        if (.Platform$OS.type == "unix") {
          Test <- system(paste("ls"), intern = TRUE)
          if (!any(Test %in% model)) 
            system(paste("mkdir",model))
        }
        if (.Platform$OS.type == "windows") {
          Test <- shell(paste("dir /b"), intern = TRUE)
          if (!any(Test %in% model)) 
            shell(paste("md ",model))
        }
        
        Vy=var(pheno.samp,na.rm=TRUE)
        Sr=Vy*R2*(DF-2)/MSx
        rate=(shape-1)/Sr
        SrB=Vy*R2*(DF-2)/MSx/probIn
        rateB=(shape-1)/SrB
        Se=Vy*(1-R2)*(DF-2)
        
        ETA=switch(model,
                   BayesA=list(list(X=geno,model="BayesA",df0=DF,rate0=rate,shape0=shape,R2=R2)),
                   BayesB=list(list(X=geno,model="BayesB",df0=DF,rate0=rateB,shape0=shape,probIn=probIn,counts=counts,R2=R2)),
                   BayesC=list(list(X=geno,model="BayesC",df0=DF,S0=Se,probIn=probIn,counts=counts,R2=R2)))
        
        
        capture.output(mod50k <- BGLR(y=as.matrix(pheno.samp),ETA=ETA,nIter=nIter,burnIn=burnIn,
                                      saveAt = paste(model,"/rep", i, "_fold", ii, sep = "")), 
                       file = paste(model,"/",model,"out_rep",i, "_fold", ii, ".txt", sep = ""))
        bu <- as.numeric(mod50k$ETA[[1]]$b)
        y.dach <- as.matrix(X2 %*% bu+mod50k$mu)
      }
      if(model=="RKHS"){
        if (.Platform$OS.type == "unix") {
          RKHSTest <- system(paste("ls"), intern = TRUE)
          if (!any(RKHSTest %in% "RKHS")) 
            system(paste("mkdir RKHS"))
        }
        if (.Platform$OS.type == "windows") {
          RKHSTest <- shell(paste("dir /b"), intern = TRUE)
          if (!any(RKHSTest %in% "RKHS")) 
            shell(paste("md RKHS"))
        }
        Vy=var(pheno.samp,na.rm=TRUE)
        Se=Vy*(1-R2)*(DF-2)
        ETA=list(list(K=K1,model='RKHS',df0=DF,S0=Se),list(K=K2,model='RKHS',df0=DF,S0=Se),list(K=K3,model='RKHS',df0=DF,S0=Se))
        capture.output(mod50k <- BGLR(y=as.matrix(pheno.samp),ETA=ETA,nIter=nIter,burnIn=burnIn,
                                      saveAt = paste("RKHS/rep", i, "_fold", ii, sep = "")), 
                       file = paste("RKHS/RKHSout_rep",i, "_fold", ii, ".txt", sep = ""))
        y.dach <- as.matrix(mod50k$yHat[rownames(pheno) %in% samp.ts[, 1]])
      }
      if(model=="randomForest"){
        unphenod=which(is.na(pheno.samp[,1]))
        phenod=setdiff(1:nrow(pheno.samp),unphenod)  
                
        mod50k<-randomForest(x=geno[phenod,],y=pheno.samp[phenod,1],xtest=geno[unphenod,],ntree=ntree,mtry=mtrys)
        y.dach <- as.matrix(mod50k$test$predicted)
      }
      rownames(y.dach)=rownames(y2)
      n.TS[ii, i] <- nrow(y.dach)
      n.DS[ii, i] <- length(namesDS)
      y.TS <- rbind(y.TS, y.dach)
      
      COR <- round(cor(y2, y.dach), digits = 4)
      COR2[ii, i] <- COR
      rCOR <- round(cor(y2, y.dach, method = "spearman"), digits = 4)
      rCOR2[ii, i] <- rCOR
      lm1 <- lm(y2 ~ as.numeric(y.dach))
      lm.coeff[ii, i] <- lm1$coefficients[2]
      mse[ii, i] <- mean((y2 - as.numeric(y.dach))^2)
      id.TS[[ii]] <- as.character(unique(samp.ts[, 1]))
      names(id.TS)[[ii]] <- paste("fold", ii, sep = "")
    }
    y.TS <- y.TS[order(rownames(y.TS)), ]
    y.TS2 <- cbind(y.TS2, y.TS)
    colnames(y.TS2)[i] <- paste("rep", i, sep = "")
    id.TS2[[i]] <- id.TS
    names(id.TS2)[[i]] <- paste("rep", i, sep = "")

    TS10 <- y.TS[order(-y.TS)]
    n10 <- round(0.1 * length(y.TS))
    TS10sel <- TS10[1:n10]
    m10[i, 1] <- mean(pheno[rownames(pheno) %in% names(TS10sel), "TRAIT"])
  
  }

  obj <- list(n.SNP = ncol(geno), n.TS = n.TS, n.DS = n.DS, 
              id.TS = id.TS2, y.TS = y.TS2, PredAbi = COR2, 
              rankCor = rCOR2, bias = lm.coeff, k = k, Rep = Rep, sampling = sampling, 
              Seed = Seed, rep.seed = seed2, VC.est.method = model,m10 = m10, mse = mse)
  class(obj) <- "cvData"
  return(obj)
}