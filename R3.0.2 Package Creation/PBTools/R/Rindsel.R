#-------------------------------------------------
# These functions are created for R Analytical Pipeline (RAP)
# Authors: S. Perez-Elizalde, J. Crossa, J. Ceron-Rojas, and G. Alvarado
#-------------------------------------------------

#Rindsel scripts

Read <- 
function (file.dat = NULL, file.wgt = NULL) 
{   
   fcol <<- 4
#    if(!is.null(traits)&&!is.null(weights)){
##    if (is.null(file.dat)) 
##        fname <- choose.files(, "Select data file")
##    else 
    fname <- file.dat
    traits1 <- read.csv(file = fname, header = T)
    traits1[, 1] <- as.factor(traits1[, 1])
    traits1[, 2] <- as.factor(traits1[, 2])
    traits1[, 3] <- as.factor(traits1[, 3])
    nom <- colnames(traits1[, -(1:(fcol - 1))])
    n <- length(nom)
##    if (is.null(file.wgt)) 
##        fname <- choose.files(, "Select weigths file")
##    else 
    fname <- file.wgt
    wgts <<- read.csv(file = fname, header = T)#}
#    else{traits1 <- traits
#         wgts <<- weights} 
    traits <- TraitsSelect(traits1, wgts[, 2])
    ntraits <<- ncol(traits) - (fcol - 1)
    return(traits)
}

Readwgt <-
function () 
{
    n <- ntraits
    nom <- wgts[as.logical(wgts[, 2]), 1]
    ewgt <- wgts[1:n, 3]
    names(ewgt) <- nom
    return(ewgt)
}

Readsgn <-
function () 
{
    n <- ntraits
    nom <- wgts[as.logical(wgts[, 2]), 1]
    sgn <- wgts[1:n, 4]
    names(sgn) <- nom
    return(sgn)
}

Readrest <-
function (file.wgt = NULL) 
{
    n <- ntraits
    nom <- wgts[as.logical(wgts[, 2]), 1]
    indrest <- wgts[1:n, 5]
    W <- diag(n)[, which(indrest == 1)]
    return(W)
}

TraitsSelect <-
function (data, ind) 
{
    traits <- data[, -(1:(fcol - 1))][, as.logical(ind)]
    cat("SELECTED TRAITS", "\n")
    cat(colnames(traits), "\n")
    traits <- cbind(data[, (1:(fcol - 1))], traits)
    return(traits)
}

adjmeans <-
function (i, data, design = "lattice") 
{
    mu <- switch(design, lattice = lmer(data[, i] ~ ENTRY + (1 | 
        REP/Block) - 1, data = data, REML = TRUE)@mu, rcb = lmer(data[, 
        i] ~ ENTRY + (1 | REP) - 1, data = data, REML = TRUE)@mu)
    matmu <- matrix(mu, ncol = length(unique(data$REP)))
    means <- apply(matmu, 1, mean)
    names(means) <- paste("Entry ", unique(data$ENTRY), sep = "")
    return(means)
}

Means <-
function (data, design = "lattice") 
{
    am <- sapply(fcol:ncol(data), adjmeans, data, design)
    am <- as.data.frame(am)
    colnames(am) <- colnames(data[, -(1:(fcol - 1))])
    return(am)
}

Adjcovs <-
function (i, data, design = "lattice") 
{
    if (design == "lattice") {
        varcorr <- VarCorr(lmer(data[, i] ~ 1 + (1 | REP:Block) + 
            (1 | ENTRY), data = data, REML = TRUE))
        cov.names <- c("block(rep)", "ENTRY", "Residual")
    }
    else if (design == "rcb") {
        varcorr <- VarCorr(lmer(data[, i] ~ 1 + (1 | REP) + (1 | 
            ENTRY), data = data, REML = TRUE))
        cov.names <- c("rep", "ENTRY", "Residual")
    }
    covs <- c(varcorr$REP, varcorr$ENTRY, attr(varcorr, "sc")^2)
    names(covs) <- cov.names
    return(covs)
}

gCovariances <-
function (data, design = "lattice") 
{
    ncols <- ntraits
    is <- unlist(sapply(1:ncols, function(i) rep(i, ncols - i + 
        1)))
    js <- unlist(sapply(1:ncols, function(i) i:ncols))
    nom <- colnames(data[, -(1:(fcol - 1))])
    adjcovij <- function(i, j, data, design) {
        coupij <- function(i, j, data) {
            data[, i] + data[, j]
        }
        coup12 <- data.frame(data[, 1:(fcol - 1)], coupij(i, 
            j, data))
        cov1 <- Adjcovs(i, data, design)[2]
        if (i != j) {
            cov2 <- Adjcovs(j, data, design)[2]
            cov3 <- Adjcovs(4, coup12, design)[2]
            covar <- (cov3 - cov1 - cov2) * 0.5
        }
        else {
            covar <- cov2 <- cov1
        }
        covij <- c(i - fcol + 1, j - fcol + 1, covar)
        names(covij) <- c("row", "col", "cov")
        return(covij)
    }
    matij <- function(i, j, data) {
        if (i == 1 & j == 1) {
            ncols <- max(data[, 1])
            mat <<- matrix(rep(0, ncols^2), ncol = ncols)
        }
        ind <- which((data[, 1] == i) & (data[, 2] == j))
        mat[i, j] <<- ifelse(i != j, data[ind, 3], data[ind, 
            3] * 0.5)
    }
    covdata <- t(mapply(adjcovij, fcol - 1 + is, fcol - 1 + js, 
        MoreArgs = list(data = data, design = design)))
    mapply(matij, is, js, MoreArgs = list(data = covdata))
    matcov <- mat + t(mat)
    dimnames(matcov) <- list(nom, nom)
    return(matcov)
}

fCovariances <-
function (data) 
{
    Cov <- cov(data)
    nom <- colnames(data)
    dimnames(Cov) <- list(nom, nom)
    return(Cov)
}

ReadMarkers <-
function(markers=NULL)
{ 
##	if(is.null(markers)) markers <- choose.files(,"Select markers data file")
	markers <- read.csv(file=markers,na.strings = ".")
	nmark <- ncol(markers)
	colnames(markers) <- paste("M",1:nmark,sep="")
	
	recode <- function(i)
	{
		cod <- ifelse(markers[,i]==0,-1, 
				ifelse(markers[,i]==1,0,
				ifelse(markers[,i]==2,1,NA)))
    	cod[which(is.na(cod))] <- 0
    	return(cod)
	}
	
	if(any(markers==2)) codes <- sapply(1:nmark,recode)  else codes <- markers
	return(as.matrix(codes))
}

ReadQTL <-
function(qtl=NULL)
{
##    if(is.null(qtl)) qtl0 <- choose.files(,"QTL file")
##    else 
    qtl0 <- qtl #added 021612
    qtl0<- read.csv(file=qtl0,skip=1)
    qtl <- qtl0[,seq(1,ncol(qtl0)-1,by=2)]
    link <- qtl0[,seq(2,ncol(qtl0),by=2)]
    return(list(QTL=qtl,Link=link))
}

CovMol <-
function(codes,QTLS)
{
    qtl <- as.matrix(QTLS$QTL)
    link <- as.matrix(QTLS$Link)
    scores <- matrix(0,nrow(codes),ncol(qtl))
    for(i in 1:ncol(qtl)){
    scores[,i] <- codes[,link[!is.na(link[,i]),i]]%*%qtl[!is.na(qtl[,i]),i]
    }
    covar <- cov(scores)
    return(list(S=covar,Scores=scores))
}

Stan <-
function (data, design = "lattice") 
{
    fstd <- function(i) {
        stdres(lm(data[, i] ~ 1))
    }
    std <- sapply(fcol:ncol(data), fstd)
    colnames(std) <- paste("V", 1:ncol(std), sep = "")
    std <- data.frame(data[, 1:(fcol - 1)], std)
    return(std)
}




SmithIndex_2 <-
function (file.dat = NULL, weights = NULL, selval = 5, design = "lattice", 
    #corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE) 
    corr = FALSE, rawdata = TRUE)
{
    library(lme4)
    library(Hmisc)
    if (rawdata == TRUE) {
        traits <<- NULL
        covG <<- NULL
        covP <<- NULL
        allmean <<- NULL
    }
    if (is.null(traits)) {
        traits <<- Read(file.dat, weights)
    }
    theta <- Readwgt()
    #cat("\n", "   Wait ...", "\n\n")
    nom.traits <- colnames(traits)[-(1:fcol - 1)]
    if (is.null(covG) & is.null(covP) & is.null(allmean)) {
        covG <<- gCovariances(traits, design = design)
        dimnames(covG) <<- list(nom.traits, nom.traits)
        allmean <<- Means(traits)
        colnames(allmean) <<- nom.traits
        covP <<- fCovariances(allmean)
        dimnames(covP) <<- dimnames(covG)
    }
    if (corr == TRUE) {
        mat.name <- "CORRELATION"
        MVGS <<- cov2cor(covG)
        MVG <- MVGS
        MVPS <<- cov2cor(covP)
        MVP <- MVPS
    }
    else {
        mat.name <- "COVARIANCE"
        MVG <- covG
        MVP <- covP
    }
    BSmith1 <- as.vector(solve(MVP) %*% MVG %*% theta)
    BSmith2 <- as.vector(sqrt(t(BSmith1) %*% BSmith1))
    BSmith <- BSmith1/BSmith2
    YSmith <<- scale(allmean) %*% BSmith
    rnames <- rownames(YSmith)
    vGb <- t(theta) %*% MVG %*% BSmith
    VIDS <- t(BSmith) %*% MVP %*% BSmith
    VDCG <- t(theta) %*% MVG %*% theta
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
    corrSMITH <- min(0.9999, vGb/(vGv %*% bPb))
    selected <- YSmith >= quantile(YSmith, 1 - selval/100)
    selYSmith <- YSmith[selected]
    ord <- order(selYSmith, decreasing = TRUE)
    selYSmith <- selYSmith[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MsmithSI <- apply(selentry, 2, mean)
    Msmithall <- apply(allmean, 2, mean)
    GAINbyMEANS = MsmithSI - Msmithall
    ks <- if (selval == 5) 
        2.063
    else if (selval == 10) 
        1.755
    else 1.4
    SMITHGain = as.vector(ks * (covG %*% BSmith)/as.numeric(sqrt(t(BSmith) %*% 
        covP %*% BSmith)))
    selentry2 <- rbind(selentry, MsmithSI, Msmithall, GAINbyMEANS, 
        SMITHGain)
    selSMITH <- data.frame(selentry2, c(selYSmith, rep(NA, 4)))
    colnames(selSMITH)[ncol(selSMITH)] <- "Smith index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selSMITH)[(length(selYSmith) + 1):(length(selYSmith) + 
        4)] <- rnames
    allSMITH <- data.frame(allmean, YSmith)
    colnames(allSMITH)[ncol(selSMITH)] <- "Smith index"
    #cat("SMITH SELECTION INDEX METHOD", "\n", file = out)
    #cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVG, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVP, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "COVARIANCE BETWEEN THE SMITH SELECTION INDEX AND THE BREEDING VALUE: ", 
    #    vGb, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE SMITH SELECTION INDEX:                               ", 
    #    VIDS, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE BREEDING VALUE:                                      ", 
    #    VDCG, "\n", file = out, append = T)
    #cat("\n", "CORRELATION BETWEEN THE SMITH SELECTION INDEX AND THE BREEDING VALUE:", 
    #    corrSMITH, "\n", file = out, append = T)
    #cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE SMITH SELECTION, MEANS AND GAINS FOR", 
    #    paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
    #    append = T)
    #print.char.matrix(round(selSMITH, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE SMITH SELECTION INDEX", 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(allSMITH, 2), file = out, col.names = T, 
    #    append = T)
    #file.show(out)
    #write.csv(selSMITH, na = "", file = outcsv)
    #write.csv(allSMITH, file = paste("all", outcsv, sep = ""))
    #cat("\n", paste("Output saved in the files:", paste(out, 
    #    ",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""), 
    #    sep = " "), "\n")
    rm(fcol, ntraits, mat, pos = ".GlobalEnv")

    detach("package:lme4")
    detach("package:Hmisc")
    return(list(
          MVGeno = round(MVG, 2),
          MVPheno = round(MVP, 2),
          CovIndBV = vGb,
          VarInd = VIDS,
          VarBV = VDCG,
          CorrIndBV = corrSMITH,
          SelGen = round(selSMITH,2),
          AllGen = round(allSMITH,2)
      ))
}



ESIMIndex_2 <-
function (file.dat = NULL, weights = NULL, selval = 5, design = "lattice", 
#    corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE) 
     corr = FALSE, rawdata = TRUE)
{
    if (rawdata == TRUE) {
        traits <<- NULL
        covG <<- NULL
        covP <<- NULL
        allmean <<- NULL
    }
    if (is.null(traits)) {
        traits <<- Read(file.dat, weights)
    }
    theta <- Readsgn()
    nom.traits <- colnames(traits)[-(1:fcol - 1)]
    #cat("\n", "   Wait ...", "\n\n")
    if (is.null(covG) & is.null(covP) & is.null(allmean)) {
        covG <<- gCovariances(traits, design = design)
        dimnames(covG) <<- list(nom.traits, nom.traits)
        allmean <<- Means(traits)
        colnames(allmean) <<- nom.traits
        covP <<- fCovariances(allmean)
        dimnames(covP) <<- dimnames(covG)
    }
    if (corr == TRUE) {
        mat.name <- "CORRELATION"
	  MVGS <<- cov2cor(covG)
        MVG <- MVGS
        MVPS <<- cov2cor(covP)
        MVP <- MVPS   
 	}
    else {
        mat.name <- "COVARIANCE"
        MVG <- covG
        MVP <- covP
    }
    svdIPG <- svd(solve(MVP) %*% MVG)
    EigVecESIM <- as.matrix(abs(svdIPG$u[, 1])) * theta
    EigVecESIM <- EigVecESIM/as.vector(sqrt(t(EigVecESIM)%*%EigVecESIM))
    YESIM <<- scale(allmean) %*% EigVecESIM
    rnames <- rownames(YESIM)
    thetaESIM <- solve(MVG) %*% (MVP) %*% EigVecESIM
    thetaNESIM <- thetaESIM %*% (1/sqrt(t(thetaESIM) %*% thetaESIM))
    vGb <- t(thetaNESIM) %*% MVG %*% EigVecESIM
    VIDS <- t(EigVecESIM) %*% MVP %*% EigVecESIM
    VDCG <- t(thetaNESIM) %*% MVG %*% thetaNESIM
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
    corrESIM <- min(0.9999, vGb/(vGv %*% bPb))
    selected <- YESIM >= quantile(YESIM, 1 - selval/100)
    selYESIM <- YESIM[selected]
    ord <- order(selYESIM, decreasing = TRUE)
    selYESIM <- selYESIM[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MESIMSI <- apply(selentry, 2, mean)
    MESIMall <- apply(allmean, 2, mean)
    GAINbyMEANS = MESIMSI - MESIMall
    ks <- if (selval == 5) 
        2.063
    else if (selval == 10) 
        1.755
    else 1.4
    ESIMGain <- as.vector(ks * (covG %*% EigVecESIM)/as.numeric(sqrt(t(EigVecESIM) %*% 
        covP %*% EigVecESIM)))
    selentry2 <- rbind(selentry, MESIMSI, MESIMall, GAINbyMEANS, 
        ESIMGain)
    selESIM <- data.frame(selentry2, c(selYESIM, rep(NA, 4)))
    colnames(selESIM)[ncol(selESIM)] <- "ESIM index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selESIM)[(length(selYESIM) + 1):(length(selYESIM) + 
        4)] <- rnames
    allESIM <- data.frame(allmean, YESIM)
    colnames(allESIM)[ncol(allESIM)] <- "ESIM index"
    #cat("ESIM SELECTION INDEX METHOD", "\n", file = out)
    #cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVG, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVP, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "COVARIANCE BETWEEN THE ESIM SELECTION INDEX AND THE BREEDING VALUE: ", 
    #    vGb, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE ESIM SELECTION INDEX:                               ", 
    #    VIDS, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE BREEDING VALUE:                                      ", 
    #    VDCG, "\n", file = out, append = T)
    #cat("\n", "CORRELATION BETWEEN THE ESIM SELECTION INDEX AND THE BREEDING VALUE:", 
    #    corrESIM, "\n", file = out, append = T)
    #cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE ESIM SELECTION, MEANS AND GAINS FOR", 
    #    paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
    #    append = T)
    #print.char.matrix(round(selESIM, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE ESIM SELECTION INDEX", 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(allESIM, 2), file = out, col.names = T, 
    #    append = T)
    #file.show(out)
    #write.csv(selESIM, na = "", file = outcsv)
    #write.csv(allESIM, file = paste("all", outcsv, sep = ""))
    #cat("\n", paste("Output saved in the files:", paste(out, 
    #    ",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""), 
    #    sep = " "), "\n")
    rm(fcol, ntraits, mat, pos = ".GlobalEnv")
    
    return(list(
          MVGeno = round(MVG,2),
          MVPheno = round(MVP,2),
          CovIndBV = vGb,
          VarInd = VIDS,
          VarBV = VDCG,
          CorrIndBV = corrESIM,
          SelGen = round(selESIM,2),
          AllGen = round(allESIM,2)
      ))
}


KNIndex_2 <-
function (file.dat = NULL, weights = NULL, selval = 5, design = "lattice", 
    #corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE) 
     corr = FALSE, rawdata = TRUE) 
{
    if (rawdata == TRUE) {
        traits <<- NULL
        covG <<- NULL
        covP <<- NULL
        allmean <<- NULL
    }
    if (is.null(traits)) {
        traits <<- Read(file.dat, weights)
    }
    thetaKN <- Readwgt()
    nom.traits <- colnames(traits)[-(1:fcol - 1)]
    W <- Readrest(weights)
    #cat("\n", "   Wait ...", "\n\n")
    if (is.null(covG) & is.null(covP) & is.null(allmean)) {
        covG <<- gCovariances(traits, design = design)
        dimnames(covG) <<- list(nom.traits, nom.traits)
        allmean <<- Means(traits)
        colnames(allmean) <<- nom.traits
        covP <<- fCovariances(allmean)
        dimnames(covP) <<- dimnames(covG)
    }

    if (corr == TRUE) {
        mat.name <- "CORRELATION"
	  MVGS <<- cov2cor(covG)
        MVG <- MVGS
        MVPS <<- cov2cor(covP)
        MVP <- MVPS   
 	}

    else {
        mat.name <- "COVARIANCE"
        MVG <- covG
        MVP <- covP
    }
    IcovP <- solve(MVP)
    covGV <- MVG %*% thetaKN
    Bs <- IcovP %*% covGV
    C <- t(MVG) %*% W
    Port <- IcovP %*% C %*% solve(t(C) %*% IcovP %*% C) %*% t(C)
    NR <- nrow(Port)
    A2 <- diag(NR)
    A <- A2 - Port
    BKNnn <- A %*% Bs
    BKN <- BKNnn/sqrt(sum(BKNnn^2))
    YKN <<- scale(allmean) %*% BKN
    rnames <- rownames(YKN)
    vGb <- t(thetaKN) %*% MVG %*% BKN
    VIDS <- t(BKN) %*% MVP %*% BKN
    VDCG <- t(thetaKN) %*% MVG %*% thetaKN
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
    corrKN <- min(0.9999, vGb/(vGv %*% bPb))
    selected <- YKN >= quantile(YKN, 1 - selval/100)
    selYKN <- YKN[selected]
    ord <- order(selYKN, decreasing = TRUE)
    selYKN <- selYKN[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MKNSI <- apply(selentry, 2, mean)
    MKNall <- apply(allmean, 2, mean)
    GAINbyMEANS = MKNSI - MKNall
    ks <- if (selval == 5) 
        2.063
    else if (selval == 10) 
        1.755
    else 1.4
    KNGain <- as.vector(ks * (covG %*% BKN)/as.numeric(sqrt(t(BKN) %*% 
        covP %*% BKN)))
    selentry2 <- rbind(selentry, MKNSI, MKNall, GAINbyMEANS, 
        KNGain)
    selKN <- data.frame(selentry2, c(selYKN, rep(NA, 4)))
    colnames(selKN)[ncol(selKN)] <- "KN index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selKN)[(length(selYKN) + 1):(length(selYKN) + 4)] <- rnames
    allKN <- data.frame(allmean, YKN)
    colnames(allKN)[ncol(allKN)] <- "KN index"
    #cat("THE KEMPTHORNE AND NORDSKOG (KN) RESTRICTIVE SELECTION INDEX", 
    #    "\n", file = out)
    #cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVG, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVP, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "COVARIANCE BETWEEN THE KN SELECTION INDEX AND THE BREEDING VALUE: ", 
    #    vGb, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE KN SELECTION INDEX:                               ", 
    #    VIDS, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE BREEDING VALUE:                                      ", 
    #    VDCG, "\n", file = out, append = T)
    #cat("\n", "CORRELATION BETWEEN THE KN SELECTION INDEX AND THE BREEDING VALUE:", 
    #    corrKN, "\n", file = out, append = T)
    #cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE KN SELECTION, MEANS AND GAINS FOR", 
    #    paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
    #    append = T)
    #print.char.matrix(round(selKN, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE KN SELECTION INDEX", 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(allKN, 2), file = out, col.names = T, 
    #    append = T)
    #file.show(out)
    #write.csv(selKN, na = "", file = outcsv)
    #write.csv(allKN, file = paste("all", outcsv, sep = ""))
    #cat("\n", paste("Output saved in the files:", paste(out, 
    #    ",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""), 
    #    sep = " "), "\n")
    rm(fcol, ntraits, mat, pos = ".GlobalEnv")

    return(list(
          MVGeno = round(MVG,2),
          MVPheno = round(MVP, 2),
          CovIndBV = vGb,
          VarInd = VIDS,
          VarBV = VDCG,
          CorrIndBV = corrKN,
          SelGen = round(selKN,2),
          AllGen = round(allKN,2)
      ))    
}


RESIMIndex_2 <-
function (file.dat = NULL, weights = NULL, selval = 5, design = "lattice", 
    #corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE) 
    corr = FALSE, rawdata = TRUE)
{
    if (rawdata == TRUE) {
        traits <<- NULL
        covG <<- NULL
        covP <<- NULL
        allmean <<- NULL
    }
    if (is.null(traits)) {
        traits <<- Read(file.dat, weights)
    }
    theta <- Readsgn()
    nom.traits <- colnames(traits)[-(1:fcol - 1)]
    W <- Readrest(weights)
    #cat("\n", "   Wait ...", "\n\n")
    if (is.null(covG) & is.null(covP) & is.null(allmean)) {
        covG <<- gCovariances(traits, design = design)
        dimnames(covG) <<- list(nom.traits, nom.traits)
        allmean <<- Means(traits)
        colnames(allmean) <<- nom.traits
        covP <<- fCovariances(allmean)
        dimnames(covP) <<- dimnames(covG)
    }

    if (corr == TRUE) {
        mat.name <- "CORRELATION"
	  MVGS <<- cov2cor(covG)
        MVG <- MVGS
        MVPS <<- cov2cor(covP)
        MVP <- MVPS   
 	}

    else {
        mat.name <- "COVARIANCE"
        MVG <- covG
        MVP <- covP
    }
    IcovP <- solve(MVP)
    C <- t(MVG) %*% W
    Port <- IcovP %*% C %*% solve(t(C) %*% IcovP %*% C) %*% t(C)
    NR <- nrow(Port)
    QI <- diag(NR)
    A2 <- QI - Port
    A <- A2 %*% IcovP %*% MVG
    BRESIM <- abs(svd(A)$u[, 1]) * theta
    thetaRESIM <- solve(MVG) %*% (MVP) %*% BRESIM
    thetanRESIM <- thetaRESIM/sqrt(sum(thetaRESIM^2))
    YRESIM <<- scale(allmean) %*% BRESIM
    rnames <- rownames(YRESIM)
    vGb <- t(thetanRESIM) %*% MVG %*% BRESIM
    VIDS <- t(BRESIM) %*% MVP %*% BRESIM
    VDCG <- t(thetanRESIM) %*% MVG %*% thetanRESIM
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
    corrRESIM <- min(0.9999, vGb/(vGv %*% bPb))
    selected <- YRESIM >= quantile(YRESIM, 1 - selval/100)
    selYRESIM <- YRESIM[selected]
    ord <- order(selYRESIM, decreasing = TRUE)
    selYRESIM <- selYRESIM[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MRESIMSI <- apply(selentry, 2, mean)
    MRESIMall <- apply(allmean, 2, mean)
    GAINbyMEANS = MRESIMSI - MRESIMall
    ks <- if (selval == 5) 
        2.063
    else if (selval == 10) 
        1.755
    else 1.4
    RESIMGain <- as.vector(ks * (covG %*% BRESIM)/as.numeric(sqrt(t(BRESIM) %*% 
        covP %*% BRESIM)))
    selentry2 <- rbind(selentry, MRESIMSI, MRESIMall, GAINbyMEANS, 
        RESIMGain)
    selRESIM <- data.frame(selentry2, c(selYRESIM, rep(NA, 4)))
    colnames(selRESIM)[ncol(selRESIM)] <- "RESIM index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selRESIM)[(length(selYRESIM) + 1):(length(selYRESIM) + 
        4)] <- rnames
    allRESIM <- data.frame(allmean, YRESIM)
    colnames(allRESIM)[ncol(allRESIM)] <- "RESIM index"
    #cat("THE RESIM SELECTION INDEX", "\n", file = out)
    #cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVG, 2), file = out, col.names = T, 
    #   append = T)
    #cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(MVP, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "COVARIANCE BETWEEN THE RESIM SELECTION INDEX AND THE BREEDING VALUE: ", 
    #    vGb, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE RESIM SELECTION INDEX:                               ", 
    #    VIDS, "\n", file = out, append = T)
    #cat("\n", "VARIANCE OF THE BREEDING VALUE:                                      ", 
    #    VDCG, "\n", file = out, append = T)
    #cat("\n", "CORRELATION BETWEEN THE RESIM SELECTION INDEX AND THE BREEDING VALUE:", 
    #    corrRESIM, "\n", file = out, append = T)
    #cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE RESIM SELECTION, MEANS AND GAINS FOR", 
    #    paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
    #    append = T)
    #print.char.matrix(round(selRESIM, 2), file = out, col.names = T, 
    #    append = T)
    #cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE RESIM SELECTION INDEX", 
    #    "\n", file = out, append = T)
    #print.char.matrix(round(allRESIM, 2), file = out, col.names = T, 
    #    append = T)
    #file.show(out)
    #write.csv(selRESIM, na = "", file = outcsv)
    #write.csv(allRESIM, file = paste("all", outcsv, sep = ""))
    #cat("\n", paste("Output saved in the files:", paste(out, 
    #    ",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""), 
    #    sep = " "), "\n")
    rm(fcol, ntraits, mat, pos = ".GlobalEnv")

    return(list(
          MVGeno = round(MVG, 2),
          MVPheno = round(MVP, 2),
          CovIndBV = vGb,
          VarInd = VIDS,
          VarBV = VDCG,
          CorrIndBV = corrRESIM,
          SelGen = round(selRESIM,2),
          AllGen = round(allRESIM,2)
      ))    
}


LTIndex_2 <-
#function(file.dat=NULL,weights=NULL,selval=5,design="lattice", corr=FALSE,out="out.txt",outcsv="out.csv",rawdata=TRUE)
function(file.dat=NULL,weights=NULL,selval=5,design="lattice",
    #corr=FALSE,out="out.txt",outcsv="out.csv",rawdata=TRUE,markers,qtls)  
    corr=FALSE,rawdata=TRUE,markers,qtls)
{
  
  if(rawdata==TRUE){
    traits <<- NULL
##    traits2 <<- NULL
    covG <<- NULL
    covP <<- NULL
    allmean <<- NULL
  }
  if(is.null(traits)){
  	traits <<- Read(file.dat,weights)
  }
  theta <- rep(Readsgn(),2)
  ##nom.traits <- colnames(traits1[,-(1:fcol-1)])[as.logical(traits2[,2])]
  #cat("\n","   Wait ...","\n\n")
  nom.traits <- colnames(traits)[-(1:fcol - 1)]
  if(is.null(covG)&is.null(covP)&is.null(allmean)){
    covG <<- gCovariances(traits,design=design)
    dimnames(covG) <<- list(nom.traits,nom.traits)
    allmean <<- Means(traits)
    colnames(allmean) <<- nom.traits
    covP <<- fCovariances(allmean)
    dimnames(covP) <<- dimnames(covG)
  }
 
  codes <- ReadMarkers(markers)   #codes <- ReadMarkers()
  QTLS <- ReadQTL(qtls)           #QTLS <- ReadQTL()
    
  cov.scores <- CovMol(codes,QTLS)
  MVMNS <<- as.matrix(cov.scores$S[1:ncol(allmean),1:ncol(allmean)])
  dimnames(MVMNS) <- dimnames(covG)
  
  if(corr==TRUE) {
    mat.name <- "CORRELATION"
    MVGS <<- cov2cor(covG)
    MVG <- MVGS
    MVPS <<- cov2cor(covP)
    MVP <- MVPS
    MVMS <<- cov2cor(MVMNS)
    MVM <- MVMS
    
  } else{
    mat.name <- "COVARIANCE"
    MVG <- covG
    MVP <- covP
    MVM <- MVMNS
  }
  
  
  MMVP <<- rbind(cbind(MVP,MVM),cbind(MVM,MVM))
  MMVG <<- rbind(cbind(MVG,MVM),cbind(MVM,MVM))
  scores <- cov.scores$Scores[,1:ncol(allmean)]
  BLT1 <<- as.vector(solve(MMVP)%*%MMVG%*%theta)
  BLT2 <<- as.vector(sqrt(t(BLT1)%*%BLT1))
  BLT <<- BLT1/BLT2
  YLT <<- cbind(scale(allmean),scale(scores))%*%BLT
  rnames <- rownames(YLT)
  thetaLT <- solve(MMVG)%*%(MMVP)%*%BLT
  thetaNLT <- thetaLT%*%(1/sqrt(t(thetaLT)%*%thetaLT))
  vGb <- t(thetaNLT)%*%MMVG%*%BLT
  ##vGb <- t(theta)%*%MMVG%*%BLT
  VIDS <- t(BLT)%*%MMVP%*%BLT
  VDCG <- t(thetaNLT)%*%MMVG%*%thetaNLT
  ##VDCG <- t(theta)%*%MMVG%*%theta
  vGv <- sqrt(abs(VDCG))
  bPb <- sqrt(abs(VIDS))
  corrLT <- min(0.9999,vGb/(vGv%*%bPb))  
  selected <- YLT>=quantile(YLT,1-selval/100)
  selYLT <- YLT[selected]
  ord <- order(selYLT,decreasing=TRUE)
  selYLT <- selYLT[ord]
  selentry <-  allmean[selected,]
  nom.selentry <- rownames(allmean)[selected]
  selentry <- data.frame(as.matrix(selentry)[ord,])
  dimnames(selentry) <- list(nom.selentry[ord],nom.traits)
  MLTSI <- apply(selentry,2,mean)
  MLTall <- apply(allmean,2,mean)
  GAINbyMEANS=MLTSI-MLTall
  ks <- if(selval==5) 2.063 else if(selval==10) 1.755 else 1.400
  LTGain <- as.vector(ks*(MMVG%*%BLT)/as.numeric(sqrt(t(BLT)%*%MMVP%*%BLT)))
  selentry2 <- rbind(selentry,MLTSI,MLTall,GAINbyMEANS,LTGain)
  selLT <- data.frame(selentry2,c(selYLT,rep(NA,4)))
  colnames(selLT)[ncol(selLT)] <- "LT index"
  rnames <- c("Mean of Selected Individuals", "Mean of all Individuals",
            'Selection Differential',paste("Expected Genetic Gain for",paste(selval,"%",sep=""),sep=" "))
  rownames(selLT)[(length(selYLT)+1):(length(selYLT)+4)] <- rnames
  allLT <- data.frame(allmean,YLT)
   colnames(allLT)[ncol(allLT)] <- "LT index"
  #cat("LANDE & THOMPSON SELECTION INDEX METHOD","\n",file=out)
  #cat("\n",paste("GENETIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(MVG,2),file=out,col.names=T,append=T)
  #cat("\n",paste("PHENOTYPIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(MVP,2),file=out,col.names=T,append=T)
  #cat("\n",paste("MOLECULAR",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(MVM,2),file=out,col.names=T,append=T)
  #cat("\n\n","COVARIANCE BETWEEN THE LT SELECTION INDEX AND THE BREEDING VALUE: ", vGb,"\n",file=out,append=T)
  #cat("\n","VARIANCE OF THE LT SELECTION INDEX:                               ", VIDS,"\n",file=out,append=T)
  #cat("\n","VARIANCE OF THE BREEDING VALUE:                                      ", VDCG,"\n",file=out,append=T)
  #cat("\n","CORRELATION BETWEEN THE LT SELECTION INDEX AND THE BREEDING VALUE:", corrLT,"\n",file=out,append=T)
  #cat("\n\n",
  #    paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE LT SELECTION, MEANS AND GAINS FOR", paste(selval,"%",
  #          sep =""),sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(selLT,2),file=out,col.names=T,append=T)
  #cat("\n\n","VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE LT SELECTION INDEX","\n",file=out,append=T)
  #print.char.matrix(round(allLT,2),file=out,col.names=T,append=T)
  #file.show(out)
  #write.csv(selLT,na="",file=outcsv)
  #write.csv(allLT,file=paste("all",outcsv,sep=""))
  #cat("\n",paste("Output saved in the files:",paste(out,",",sep=""),outcsv,"and",paste("all",outcsv,sep=""),sep=" "),"\n")
  rm(fcol, ntraits, mat, pos = ".GlobalEnv")
  
  return(list(MVGeno = round(MVG,2), 
              MVPheno = round(MVP,2),
              MVMole = round(MVM,2),
              CovIndBV = vGb,
              VarInd = VIDS,
              VarBV = VDCG,
              CorrIndBV = corrLT,
              SelGen = round(selLT,2),
              AllGen = round(allLT,2)
              ))
}


MESIMIndex_2 <-
#function(file.dat=NULL,weights=NULL,selval=5,design="lattice",corr=FALSE,out="out.txt",outcsv="out.csv",rawdata=TRUE){
#function(file.dat=NULL,weights=NULL,selval=5,design="lattice",corr=FALSE,out="out.txt",outcsv="out.csv",rawdata=TRUE,
#           markers, qtls)
function(file.dat=NULL,weights=NULL,selval=5,design="lattice",corr=FALSE,rawdata=TRUE,markers, qtls)
{ 
  if(rawdata==TRUE){
    traits <<- NULL
    traits2 <<- NULL
    covG <<- NULL
    covP <<- NULL
    allmean <<- NULL
  }
  if(is.null(traits)){
  traits <<- Read(file.dat,weights)
  }
  theta <- rep(Readsgn(),2)
  nom.traits <- colnames(traits)[-(1:fcol - 1)]
  #cat("\n","   Wait ...","\n\n")
  if(is.null(covG)&is.null(covP)&is.null(allmean)){
    covG <<- gCovariances(traits,design=design)
    dimnames(covG) <<- list(nom.traits,nom.traits)
    allmean <<- Means(traits)
    colnames(allmean) <<- nom.traits
    covP <<- fCovariances(allmean)
    dimnames(covP) <<- dimnames(covG)
  }
 
  codes <- ReadMarkers(markers)   #codes <- ReadMarkers()
  QTLS <- ReadQTL(qtls)               #QTLS <- ReadQTL()
    
  cov.scores <- CovMol(codes,QTLS)
  MVMNS <<- as.matrix(cov.scores$S[1:ncol(allmean),1:ncol(allmean)])
  dimnames(MVMNS) <- dimnames(covG)
  
  if(corr==TRUE) {
    mat.name <- "CORRELATION"
    MVGS <<- cov2cor(covG)
    MVG <- MVGS
    MVPS <<- cov2cor(covP)
    MVP <- MVPS
    MVMS <<- cov2cor(MVMNS)
    MVM <- MVMS
    
  } else{
    mat.name <- "COVARIANCE"
    MVG <- covG
    MVP <- covP
    MVM <- MVMNS
  }
  
 
  MMVP <<- rbind(cbind(MVP,MVM),cbind(MVM,MVM))
  MMVG <<- rbind(cbind(MVG,MVM),cbind(MVM,MVM))
  svdIPG <- svd(solve(MMVP)%*%MMVG)
  EigVecMESIM <- as.matrix(abs(svdIPG$u[,1]))*theta
  EigVecMESIM <- EigVecMESIM/as.vector(sqrt(t(EigVecMESIM)%*%EigVecMESIM))
  scores <- cov.scores$Scores[,1:ncol(allmean)]
  YMESIM <<- cbind(scale(allmean),scale(scores))%*%EigVecMESIM
  rnames <- rownames(YMESIM)
  thetaMESIM <- solve(MMVG)%*%(MMVP)%*%EigVecMESIM
  thetaNMESIM <- thetaMESIM%*%(1/sqrt(t(thetaMESIM)%*%thetaMESIM))
  vGb <- t(thetaNMESIM)%*%MMVG%*%EigVecMESIM
  VIDS <- t(EigVecMESIM)%*%MMVP%*%EigVecMESIM
  VDCG <- t(thetaNMESIM)%*%MMVG%*%thetaNMESIM
  vGv <- sqrt(abs(VDCG))
  bPb <- sqrt(abs(VIDS))
  corrMESIM <- min(0.9999,vGb/(vGv%*%bPb))  
  selected <- YMESIM>=quantile(YMESIM,1-selval/100)
  selYMESIM <- YMESIM[selected]
  ord <- order(selYMESIM,decreasing=TRUE)
  selYMESIM <- selYMESIM[ord]
  selentry <-  allmean[selected,]
  nom.selentry <- rownames(allmean)[selected]
  selentry <- data.frame(as.matrix(selentry)[ord,])
  dimnames(selentry) <- list(nom.selentry[ord],nom.traits)
  MMESIMSI <- apply(selentry,2,mean)
  MMESIMall <- apply(allmean,2,mean)
  GAINbyMEANS=MMESIMSI-MMESIMall
  ks <- if(selval==5) 2.063 else if(selval==10) 1.755 else 1.400
  MESIMGain <- as.vector(ks*(MMVG%*%EigVecMESIM)/as.numeric(sqrt(t(EigVecMESIM)%*%MMVP%*%EigVecMESIM)))
  selentry2 <- rbind(selentry,MMESIMSI,MMESIMall,GAINbyMEANS,MESIMGain)
  selMESIM <- data.frame(selentry2,c(selYMESIM,rep(NA,4)))
  colnames(selMESIM)[ncol(selMESIM)] <- "MESIM index"
  rnames <- c("Mean of Selected Individuals", "Mean of all Individuals",
            'Selection Differential',paste("Expected Genetic Gain for",paste(selval,"%",sep=""),sep=" "))
  rownames(selMESIM)[(length(selYMESIM)+1):(length(selYMESIM)+4)] <- rnames
  allMESIM <- data.frame(allmean,YMESIM)
   colnames(allMESIM)[ncol(allMESIM)] <- "MESIM index"
  #cat("MESIM SELECTION INDEX METHOD","\n",file=out)
  #cat("\n",paste("GENETIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(MVG,2),file=out,col.names=T,append=T)
  #cat("\n",paste("PHENOTYPIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(MVP,2),file=out,col.names=T,append=T)
  #cat("\n","MOLECULAR COVARIANCE MATRIX","\n",file=out,append=T)
  #print.char.matrix(round(MVM,2),file=out,col.names=T,append=T)
  #cat("\n\n","COVARIANCE BETWEEN THE MESIM SELECTION INDEX AND THE BREEDING VALUE: ", vGb,"\n",file=out,append=T)
  #cat("\n","VARIANCE OF THE MESIM SELECTION INDEX:                               ", VIDS,"\n",file=out,append=T)
  #cat("\n","VARIANCE OF THE BREEDING VALUE:                                      ", VDCG,"\n",file=out,append=T)
  #cat("\n","CORRELATION BETWEEN THE MESIM SELECTION INDEX AND THE BREEDING VALUE:", corrMESIM,"\n",file=out,append=T)
  #cat("\n\n",
  #    paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE MESIM SELECTION, MEANS AND GAINS FOR", paste(selval,"%",
  #          sep =""),sep=" "),"\n",file=out,append=T)
  #print.char.matrix(round(selMESIM,2),file=out,col.names=T,append=T)
  #cat("\n\n","VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE MESIM SELECTION INDEX","\n",file=out,append=T)
  #print.char.matrix(round(allMESIM,2),file=out,col.names=T,append=T)
  #file.show(out)
  #write.csv(selMESIM,na="",file=outcsv)
  #write.csv(allMESIM,file=paste("all",outcsv,sep=""))
  #cat("\n",paste("Output saved in the files:",paste(out,",",sep=""),outcsv,"and",paste("all",outcsv,sep=""),sep=" "),"\n")
  rm(fcol, ntraits, mat, pos = ".GlobalEnv")

  return(list(
      MVGeno = round(MVG,2),
      MVPheno = round(MVP,2),
      MVMole = round(MVM,2),
      CovIndBV = vGb,
      VarInd = VIDS,
      VarBV = VDCG,
      CorrIndBV = corrMESIM,
      SelGen = round(selMESIM,2),
      AllGen = round(allMESIM,2)
    ))
    
}