load.cross.data<-function(P.data, G.data, map.data, cross, heterozygotes=TRUE, genotype, env.label = NULL, env=NULL, ST=NULL){
  library(qtl)  

  #Import the genotypic data file
  G.data<-read.table(G.data, header=TRUE,as.is=TRUE,check.names=FALSE)
  
  #Import the map data file
  map.data<-read.table(map.data, header=FALSE,as.is=TRUE)

  if ((is.null(env.label)==FALSE) & (is.null(env)==FALSE)) {P.data<-P.data[which(P.data[match(env.label,names(P.data))]==env),]}
  
  ####remove parents from genotype data and  match geno and pheno data
  if(cross!="am"){G.data<-G.data[3:nrow(G.data),]}
  
  #   P.data<-P.data[match(rownames(G.data),P.data$genotype),]
  P.data<-P.data[match(rownames(G.data),as.character(P.data[[match(genotype,names(P.data))]])),]
  a<-matrix(rep("", (ncol(P.data)*2)), 2, ncol(P.data))
  colnames(a)<-colnames(P.data)
  names(a)<-names(P.data)
  b<-rbind(G.data[,0], t(map.data[,2:3]))
  names(b)<-as.character(map.data[,1])
  c<-data.frame(a,b, check.names=FALSE, stringsAsFactors=FALSE)
  rownames(c) <- c(rownames(a), rownames(b))
  
  #extract id information for geno and pheno
  if(cross!="am"){
    G.id<-row.names(G.data)
    P.id <- as.character(P.data[[match(genotype,names(P.data))]])
    if (sum(G.id!=P.id)>0) simpleError("IDs don't match")
    geno.data<-G.data
    geno.data[geno.data=="1"]<-"AA"
    geno.data[geno.data=="1/1"]<-"AA"
    geno.data[geno.data=="2"]<-"BB"
    geno.data[geno.data=="2/2"]<-"BB"
    geno.data[geno.data=="1/2"]<-"AB"
    geno.data[geno.data=="1/-"]<-"C"
    geno.data[geno.data=="2/-"]<-"D"
    geno.data[geno.data=="-"]<-NA
    e<-cbind(P.data,  geno.data[1:nrow(geno.data),])
    f<-rbind(c,e)
    names(f)[1]="id"
    zz <- file("cross.data.file.csv", "w") 
    write.csv(f,zz, quote=FALSE, row.names=FALSE)
    close(zz)
  }
  
  if(cross=="am"){
    G.id<-row.names(G.data)
    P.id <- as.character(P.data[[match(genotype,names(P.data))]])
    if (sum(G.id!=P.id)>0) simpleError("IDs don't match")
    geno.data<-G.data
    geno.data[geno.data=="1"]<-"BB"
    geno.data[geno.data=="0"]<-"AA"
    geno.data[geno.data=="0.5"]<-"AB"
    geno.data[geno.data=="-"]<-NA
    d<-rbind(G.data[0,],geno.data)
    names(d)<-names(G.data)
    e<-cbind(P.data, d)
    for(i in 2:dim(P.data)[2]){
      e[,i]<-as.factor(e[,i])
    }
    f<-rbind(c,e)
    names(f)[1]="id"
    zz <-file("cross.data.file.csv", "w") 
    write.csv(f,zz, quote=FALSE, row.names=FALSE)
    close(zz)
  }
  
  if(cross=="am" & heterozygotes=="FALSE") cross<-"dh"
  if(cross=="am" & heterozygotes=="TRUE") cross<-"f2"
  
  #this is to import the cross
  if(cross=="dh"){data.prueba<-read.cross("csv", file="cross.data.file.csv", genotypes=c("AA", "BB"))} #check if more
  if(cross=="bc"){data.prueba<-read.cross("csv", file="cross.data.file.csv", genotypes=c("AA", "AB"))} #check if more
  if(cross=="f2"){data.prueba<-read.cross("csv", file="cross.data.file.csv", genotypes=c("AA", "AB", "BB", "C", "D"))}
  if(cross=="riself" | cross=="ri4self" | cross=="ri8self"){data.prueba<-read.cross("csv", file="cross.data.file.csv", genotypes=c("AA", "BB", "C", "D"))}
  if(cross=="risib" | cross=="ri4sib" | cross=="ri8sib"){data.prueba<-read.cross("csv", file="cross.data.file.csv", genotypes=c("AA", "BB", "C", "D"))}
  
  class(data.prueba)[1]=paste(cross)
  
  #cross.data<-NULL
  #cross.data<-data.prueba
  
  detach("package:qtl")
  return(data.prueba)
  
} #end of load.cross.data function 

MQ.marker.diag<-function(crossobj, QTL.path, estmarker=FALSE, thresholdMQ=0.1,quantile=FALSE,p.val=0.01,na.cutoff=0.1){
  library(qtl)
  #proj.path<-paste(getwd(),sep="")
  #QTL.path <- paste(proj.path, "/", sep="")
  
  MQ.cf.file <- paste(QTL.path, "MQ_cfplot", sep="")
  MQ.markermap.file <- paste(QTL.path, "MQ_markermapplot", sep="") 
  MQ.estmarkermap.file <- paste(QTL.path, "MQ_estmarkermapplot", sep="")
  MQ.genotype.file <- paste(QTL.path, "MQ_genotypeplot", sep="")
  MQ.missinggenotype.file <- paste(QTL.path, "MQ_missingegnotypeplot", sep="")
  MQ.comparegenotypes.file <- paste(QTL.path, "MQ_comparegenotypesplot", sep="")
  MQ.identical.genotypes<- paste(QTL.path, "MQ_identical.genotypes.plot", sep="")
  MQ.genotypic.distortion<- paste(QTL.path, "MQ_genotypic.distortion.plot", sep="")
  MQ.summary.markers.file<- paste(QTL.path, "MQ_summary.markers.report", sep="")
  MQ.problems.markers.file<- paste(QTL.path, "MQ_problems.markers.report", sep="")
  
  MQ.markermap.plot <- function(crossobj, out=MQ.markermap.file){
    #To plot markermap 
    #pdf(file = paste(out,".pdf",sep=""), onefile=FALSE)
    png(file = paste(out,".png",sep=""))
    par(mfrow = c(1,1))
    plot.map(crossobj)
    dev.off()
  }
  
  MQ.estmarkermap.plot <- function(crossobj, out=MQ.estmarkermap.file){
    #To estimate and plot new markermap 
    newmap <- est.map(crossobj)
    #pdf(file = paste(out,".pdf",sep=""), onefile=FALSE)
    png(file = paste(out,".png",sep=""))
    par(mfrow = c(1,1))
    plot.map(newmap, crossobj)
    dev.off()
  }
  
  MQ.genotype.plot <- function(crossobj, out=MQ.genotype.file){
    #pdf(file = paste(out,".pdf",sep=""), onefile=FALSE)
    png(file = paste(out,".png",sep=""))
    par(mfrow = c(1,1))
    geno.image(crossobj)
    dev.off()
  }
  
  MQ.missinggenotype.plot <- function(crossobj, out=MQ.missinggenotype.file){
    #pdf(file = paste(out,".pdf",sep=""), onefile=FALSE)
    png(file = paste(out,".png",sep=""))
    par(mfrow = c(1,1))
    plot.missing(crossobj)
    dev.off()
  }
  
  MQ.comparegenotypes.plot <- function(crossobj, out=MQ.comparegenotypes.file){
    #To compare pairs of genotypes
    #pdf(file = paste(out,".pdf",sep=""), onefile=FALSE)
    png(file = paste(out,".png",sep=""))
    par(mfrow = c(1,1))
    output <- comparegeno(crossobj)
    n.ind <- nind(crossobj)
    image(1:n.ind, 1:n.ind, output, col=gray((0:99)/99),
          breaks=seq(0,1,len=101))
    dev.off()
  }
  
  MQ.cf.plot <- function(crossobj, out=MQ.cf.file){
    #To compare pairs of genotypes
    #pdf(file = paste(out,".pdf",sep=""), onefile=FALSE)
    png(file = paste(out,".png",sep=""))
    par(mfrow = c(1,1))
    plot.rf(crossobj, alternate.chrid=FALSE)
    dev.off()
  }
  
  ######create reports on problems etc
  
  ######For MISSING data exploration  

  MQ.summary.markers<-function(crossobj, p.val=0.01,na.cutoff=0.1, out=MQ.genotypic.distortion){
    
    #To add chromosome and position to markers
    names.lg<-names(pull.map(crossobj))
    num.lg<-length(names.lg)
    a<-t(t(unlist(pull.map(crossobj))))
    library(stringr)
    row.names(a)<-str_sub(row.names(a),3)
    ch<-NULL
    for(i in 1:num.lg){
      ch<-rbind(ch, t(t(rep(i,summary.map(crossobj)[i,1]))))
    } 
    b<-cbind(ch,a)
    
    #To list markers with more than 10% missing data 
    n.missing<-nmissing(crossobj, what="mar")[(nmissing(crossobj, what="mar"))/sum(summary(crossobj)$n.ind)>na.cutoff]  
    p.missing<-((nmissing(crossobj, what="mar"))/sum(summary(crossobj)$n.ind))[(nmissing(crossobj, what="mar"))/sum(summary(crossobj)$n.ind)>na.cutoff]
    
    if(length(n.missing)==0){n.missing<-NA; names(n.missing)<-"NA"} 
    if(length(p.missing)==0){p.missing<-NA; names(p.missing)<-"NA"} 
    
    Chr<-NA
    Pos<-NA
    for(i in 1: dim(b)[1]){
      for(j in 1:length(n.missing)){
        if(row.names(b)[i]==names(n.missing)[j]) {
          c<-b[i,1] 
          p<-b[i,2]
        } else {
          c<-NULL 
          p<-NULL
        }
        Chr<-rbind(Chr, c)
        Pos<-rbind(Pos, p)
      }
    }
    missing.markers<-data.frame(Chr, Pos, n.missing, p.missing, check.names=FALSE, row.names=names(n.missing))
    names(missing.markers)<-c("Chr.", "Pos.", "Num. missing", "Frec. missing")
    
    #To list individuals with more than cutoff% missing data
    n.missing.i<-nmissing(crossobj, what="ind")[(nmissing(crossobj, what="ind"))/sum(summary(crossobj)$n.mar)>na.cutoff]  
    p.missing.i<-((nmissing(crossobj, what="ind"))/sum(summary(crossobj)$n.mar))[(nmissing(crossobj, what="ind"))/
      sum(summary(crossobj)$n.mar)>na.cutoff]
    
    if(length(n.missing.i)==0){n.missing.i<-NA; names(n.missing.i)<-"NA"} 
    if(length(p.missing.i)==0){p.missing.i<-NA; names(p.missing.i)<-"NA"} 
    
    missing.individuals<-data.frame(n.missing.i, p.missing.i)
    names(missing.individuals)<-c("Num. missing", "Frec. missing")
    
    #To check for segregation distortion
    gt<-geno.table(crossobj)
    
    segregation.distortion<-gt[gt$P.value < p.val,]
    gt2<-data.frame(gt, ch, a, -log10(gt$P.value))
    
    #pdf(file=paste(out,".pdf",sep=""), onefile=FALSE)
    png(file=paste(out,".png",sep=""))
    if(nchr(crossobj)<4) par(mfrow=c(1,nchr(crossobj)))
    if(nchr(crossobj)>4 & nchr(crossobj)<8) par(mfrow=c(2,((round(nchr(crossobj))/2)+1)))
    if(nchr(crossobj)>8 ) par(mfrow=c(3,((round(nchr(crossobj))/3)+1)))
    for(i in 1:nchr(crossobj)){
      c<-gt2[,7][gt2[,1]==i]
      d<-gt2[,8][gt2[,1]==i]
      e<-cbind(c, d)
      plot(e[,1], e[,2], type="h", xlab="position", ylab="-logP", main=paste("Chr", i), ylim=c(0,max(gt2[,8])))
    }
    dev.off()
    
    out<-NULL
    out$missing.markers<-missing.markers
    out$missing.individuals<-missing.individuals
    out$segregation.distortion<-segregation.distortion
    
    filename<-MQ.summary.markers.file
    
    write("",file=filename)
    
    write("-----------------------------------------------", append=TRUE,file=filename)
    
    sink(file=filename,append=TRUE)
    print(out)
    sink()
  }

  MQ.list.problems<-function(crossobj, thresholdMQ=FALSE, quant=FALSE,out=MQ.identical.genotypes){
    #pdf(file=paste(out,".pdf",sep=""), onefile=FALSE)
    png(file=paste(out,".png",sep=""))
    par(mfrow=c(1,1))
    cg<- comparegeno(crossobj)
    hist(cg, breaks=200, xlab="Proportion of shared alleles", col="red", main="Identical genotypes")
    rug(cg, col="blue")
    dev.off()
    
    if(quant!="FALSE") {
      a<-which(cg>quantile(cg,(1-quant), na.rm=TRUE), arr.ind=TRUE)
      c<-NULL
      if(length(a)>2){
        for(i in 1:dim(a)[1]){
          b<-cg[a[i,1], a[i,2]]
          c<-rbind(c,b)
        }
        identical.genotypes<-cbind(a,c)
      }
      if(length(a)==0){
        identical.genotypes<-"No unusually identical genotypes"
      }
      
      a<-which(cg<quantile(cg, (quant), na.rm=TRUE), arr.ind=TRUE)
      c<-NULL
      if(length(a)>2){
        for(i in 1:dim(a)[1]){
          b<-cg[a[i,1], a[i,2]]
          c<-rbind(c,b)
        }
        different.genotypes<-cbind(a,c)
      }
      if(length(a)==0){
        different.genotypes<-"No unusually different genotypes"
      }
    }
    if(thresholdMQ!="FALSE") {
      a<-which(cg>(1-thresholdMQ), arr.ind=TRUE)
      c<-NULL
      if(length(a)>2){
        for(i in 1:dim(a)[1]){
          b<-cg[a[i,1], a[i,2]]
          c<-rbind(c,b)
        }
        identical.genotypes<-cbind(a,c)
      }
      if(length(a)==0){
        identical.genotypes<-"No unusually identical genotypes"
      }
      
      a<-which(cg<thresholdMQ, arr.ind=TRUE)
      c<-NULL
      if(length(a)>2){
        for(i in 1:dim(a)[1]){
          b<-cg[a[i,1], a[i,2]]
          c<-rbind(c,b)
        }
        different.genotypes<-cbind(a,c)
      }
      if(length(a)==0){
        different.genotypes<-"No unusually different genotypes"
      }
    }
    
    #Check marker order
    crossobj<-est.rf(crossobj)
    marker.order<-checkAlleles(crossobj)
    out<-NULL
    out$identical.genotypes<-identical.genotypes
    out$different.genotypes<-different.genotypes
    out$marker.order<-marker.order
    #out
    
    filename<-MQ.problems.markers.file
    
    write("",file=filename)
    
    write("-----------------------------------------------", append=TRUE,file=filename)
    sink(file=filename,append=TRUE)
    print(out)
    sink()
    
  }
  
  MQ.markermap.plot(cross.data)
  
  #reestimate a new marker map and plot it -- WARNING: computationally time consuming
  if (estmarker==TRUE){MQ.estmarkermap.plot(cross.data)}
  
  #plot alleles of the genotypes
  MQ.genotype.plot(cross.data)
  
  #plot about the missing genotypes
  MQ.missinggenotype.plot(cross.data)
  
  #comparison of genotypes
  MQ.comparegenotypes.plot(cross.data)
  
  #pairwise comparion of genotypes in terms of recombination fraction
  suppressWarnings(MQ.cf.plot(cross.data))  
  
  MQ.list.problems(cross.data,thresholdMQ=thresholdMQ)
  
  MQ.summary.markers(cross.data,p.val=p.val,na.cutoff=na.cutoff)
  
  detach("package:qtl")
  detach("package:stringr")
} #end of MQ.marker.diag function 

# QTL Functions 

QTL.analysis<-function(crossobj=cross.data, QTL.path, env.label=NULL, env = NULL, trait="pred", step, method, threshold, distance, cofactors, window.size){
  library(qtl)
  library(lattice)
  #proj.path<-paste(getwd(),sep="")
  #QTL.path <- paste(proj.path, "/", sep="")
  
  u.e.=FALSE #obsolete
  file<-crossobj  
  
  #gen.predictors 
  file<-calc.genoprob(file, step=step)
  
  #replace pseudomarkernames with proper ones
  mark.names<-c(); for(i in 1:length(file$geno)){mark.names<-c(mark.names,colnames(file$geno[[i]]$data));}#make list of markernames 
  
  #extract genotypic predictors
  probs1<-NULL
  probs2<-NULL 
  for(i in 1:nchr(file)){
    names<-dimnames(file$geno[[i]]$prob)[[2]]  
    sel.mark<-which(names%in%mark.names==FALSE); names[sel.mark]<-paste(i,names[sel.mark],sep="_")
    dimnames(file$geno[[i]]$prob)[[2]]<-names	 
    names(attributes(file$geno[[i]]$prob)$map)<-names
    p1<-file$geno[[i]]$prob
    names(p1)<-dimnames(p1)[[2]]
    probs1<-cbind(probs1, p1[,,1])
    if(class(file)[1]=="f2"|class(file)[1]=="4way"){ 
      probs2<-cbind(probs2, p1[,,3])
    }
    if(class(file)[1]=="dh"| class(file)[1]=="bc" | class(file)[1]=="riself" | class(file)[1]=="ri4self" | class(file)[1]=="ri8self" | class(file)[1]=="risib" | class(file)[1]=="ri4sib" | class(file)[1]=="ri8sib" ){ 
      probs2<-cbind(probs2, p1[,,2])  #cbind(probs,, p1[,,2]) in original , check! JvH
    }
  }
  additive<-probs2-probs1
  
  ###
  
  new.additive=additive   #JvH single phenotype check if these are supposed to be reps?
  
  #Import the phenotypic data file
  P.data<-file$pheno
  
  #use some labels...
  Gen<-"id"
  
  Trait<-paste(trait)
  
  if (is.null(env.label) == FALSE) {
    if(is.null(P.data[[match(env.label,names(P.data))]])==FALSE){ENV<-as.factor(P.data[[match(env.label,names(P.data))]])}
    if(is.null(P.data[[match(env.label,names(P.data))]])==TRUE){ENV<-as.factor(rep("env",nrow(P.data)))}
  } else {
    ENV<-as.factor(rep("env",nrow(P.data)))
  }
  
  GEN<-as.factor(P.data[, Gen])  
  MEAN<-as.numeric(as.matrix(P.data[,Trait])) 
  if(u.e.=="FALSE") {
    all.means<-data.frame(ENV,GEN,MEAN, stringsAsFactors=FALSE)      
  }
  if(u.e.!="FALSE") {
    ue<-paste(u.e.)
    UE<-as.numeric(as.matrix(P.data[,ue])) 
    all.means<-data.frame(ENV,GEN,MEAN,UE, stringsAsFactors=FALSE)      
  }
  
  b<-matrix(,0,2) 
  for(i in 1:nchr(file)){
    a<-paste("file$geno$'", i, "'$prob", sep="")
    mp<-attributes(eval(parse(text=a)))$map
    mp<-cbind(rep(i,length(mp)),as.matrix(mp))	
    b<-rbind(b,mp) 	
  }
  
  #################For MB and SIM
  if(method=="SIM"){
    p.values<-NULL
    fixeff<-NULL
    for( i in 1:dim(new.additive)[2]){
      marker<-new.additive[,i]
      GE.data<-data.frame(all.means, marker)
      if(length(levels(ENV))>1) {CS.Model<-lmer(MEAN~ ENV+ marker:ENV + (1|GEN), data=GE.data)}
      if(length(levels(ENV))==1) {CS.Model<-lm(MEAN~ marker, data=GE.data);fstats<-as.vector(summary(CS.Model)$fstatistic)}
      
      if(length(levels(ENV))>1){p.value<-1-pf(anova(CS.Model)[2,4], anova(CS.Model)[2,1], (dim(GE.data)[1]-anova(CS.Model)[1,1]-anova(CS.Model)[2,1]-1))}
      if(length(levels(ENV))==1){p.value<-1-pf(fstats[1],fstats[2],fstats[3])}   
      
      p.value<-data.frame(rownames(b)[i], b[i,1], b[i,2], p.value, stringsAsFactors=FALSE)
      p.values<-rbind(p.values, p.value)
      if(length(levels(ENV))>1){ f<-t(t(CS.Model@fixef[(length(unique(all.means$ENV))+1):length(CS.Model@fixef)]))}
      if(length(levels(ENV))==1){f<-CS.Model$coefficients[2]}
      fixeff<-cbind(fixeff, f)
    }
  }
  
  #################For CIM
  if(method=="CIM"){   ###CIM
    cofactor.list<-NULL
    cofactor.pos<-NULL
    for(i in 1:length(cofactors)){
      c<-new.additive[,which(colnames(new.additive)==cofactors[i])]
      cofactor.list<-cbind(cofactor.list, c)
      cof.pos<-cbind(cofactors[i], t(b[row.names(b)==cofactors[i]]))
      cofactor.pos<-rbind(cofactor.pos, cof.pos)
    }
    
    cofactor.win.f<-rep(0, dim(b)[1])
    for(j in 1:length(cofactors)){
      c<-b[,1]==as.numeric(as.character(cofactor.pos[j,2]))
      c[c==FALSE]=0
      c[c==TRUE]=1
      win<-(b[,2]>(as.numeric(as.character(cofactor.pos[j,3]))[1]-0.5*window.size)& b[,2]<(as.numeric(as.character(cofactor.pos[j,3]))[1]+0.5*window.size))
      win[win==FALSE]=0
      win[win==TRUE]=1
      cofactor.win<-c*win
      cofactor.win[cofactor.win==1]=j
      cofactor.win.f<-cofactor.win.f+cofactor.win
    }
    
    p.values<-NULL
    fixeff<-NULL
    for( i in 1:dim(new.additive)[2]){
      marker<-new.additive[,i]
      if(unlist(cofactor.win.f[i])>0){
        new.list<-cofactor.list[,-c(unlist(cofactor.win.f[i]))]
        number.cofactors<-length(cofactors)-1
      }
      if(unlist(cofactor.win.f[i])==0){ 
        new.list<-cofactor.list
        number.cofactors<-length(cofactors)
      }
      GE.data<-data.frame(all.means, marker)
      if(length(levels(ENV))>1) {CS.Model<-lmer(MEAN~ ENV+ marker:ENV + new.list + (1|GEN), data=GE.data)}
      if(length(levels(ENV))==1) {CS.Model<-lm(MEAN~ new.list+marker, data=GE.data);fstats<-c(anova(CS.Model)[2,4],anova(CS.Model)[2,1],anova(CS.Model)[3,1])}
      
      if(length(levels(ENV))>1) {p.value<-1-pf(anova(CS.Model)[dim(anova(CS.Model))[1],4], anova(CS.Model)[dim(anova(CS.Model))[1],1], (dim(GE.data)[1]-dim(GE.data)[2]-number.cofactors-anova(CS.Model)[dim(anova(CS.Model))[1],1]))}
      if(length(levels(ENV))==1){p.value<-1-pf(fstats[1],fstats[2],fstats[3])}   
      p.value<-data.frame(rownames(b)[i], b[i,1], b[i,2], p.value, stringsAsFactors=FALSE)
      
      p.values<-rbind(p.values, p.value)
      if(length(levels(ENV))>1){ f<-t(t(CS.Model@fixef[(length(unique(all.means$ENV))+(length(new.list)/length(all.means$ENV))+1):length(CS.Model@fixef)]))}
      if(length(levels(ENV))==1){f<-CS.Model$coefficients[ncol(new.list)+2]}
      fixeff<-cbind(fixeff, f)
    }
  }  #CIM
  
  names(p.values)=c("marker", "Chr", "Pos", "LOD")
  outem<-p.values
  
  #Threshold options
  if(threshold>0) threshold.f<-threshold
  if(threshold=="Li&Ji"){
    a<-eigen(cor(probs1), only.values=TRUE)
    a<-a$values
    for(i in 1:length(a)){
      if(a[i]>1) {a[i]<-a[i]
      } else {a[i]<-0 }
    }
    c<-NULL
    for(i in 1:length(a)){
      bt<-((a[i]-1)^2)/(totmar(file)-1)
      c<-rbind(c,bt)
    }
    v.lambda<-sum(c)
    M.eff<-1+((totmar(file)-1)*(1-(v.lambda/totmar(file))))
    alpha.e<-0.05
    alpha.p<-1-((1-alpha.e)^(1/M.eff))
    #-log(alpha.p)
    if(class(file)[1]=="bc") dff<-2
    if(class(file)[1]=="dh") dff<-2
    if(class(file)[1]=="f2") dff<-3
    if(class(file)[1]=="ril") dff<-3
    #threshold.f<-qchisq((1-alpha.p), df=dff)/4.6
    threshold.f<-alpha.p
  }
  
  ####New Looping function JvH
  pot.qtl<-outem[which(outem$'LOD'<threshold.f),] #potential qtl under threshold
  res.qtl<-c()   #result frame with selected markers
  
  if(nrow(pot.qtl)>0){ #condition to deal with 0 qtl
    pot.qtl$select<-1   #select flag to be modified in loop
    pot.qtl$eval<-0   #evaluated flag to be modified in loop
    
    #loop through chromosomes
    for(chr in unique(pot.qtl$Chr)){
      t.pot.qtl<-pot.qtl[pot.qtl$Chr==chr,]	
      while(sum(t.pot.qtl$eval)<nrow(t.pot.qtl)){
        min.p<-min(t.pot.qtl$'LOD'[which(t.pot.qtl$eval==0)])
        sel.row<-which(t.pot.qtl$'LOD'==min.p&t.pot.qtl$eval==0)[1]
        d<-abs(t.pot.qtl$Pos-t.pot.qtl$Pos[sel.row])
        t.pot.qtl$select[d<=distance]<-0
        t.pot.qtl$eval[d<=distance]<-1
        t.pot.qtl$select[sel.row]<-1
      }
      t.pot.qtl<-t.pot.qtl[which(t.pot.qtl$select==1),]
      res.qtl<-rbind(res.qtl,t.pot.qtl)	
    }
    
    res.qtl$select<-NULL
    res.qtl$eval<-NULL	
  }
  
  ##### New Looping function JvH
  if(length(levels(ENV))>1)   {  ####test for interactions if multiple environments
    ##############for second part model
    if(method=="SIM"){
      new.additive2<-new.additive[,outem[,4]<threshold.f]
      p.values.main<-NULL
      p.values.inter<-NULL
      fixeff2<-NULL
      for( i in 1:dim(new.additive2)[2]){
        marker<-new.additive2[,i]
        GE.data<-data.frame(all.means, marker)
        CS.Model<-lmer(MEAN~ ENV + marker + marker:ENV + (1|GEN), data=GE.data)
        
        p.value.inter<-1-pf(anova(CS.Model)[3,4], anova(CS.Model)[3,1], (dim(GE.data)[1]-anova(CS.Model)[1,1]-anova(CS.Model)[2,1]-2))
        p.values.inter<-rbind(p.values.inter, p.value.inter)
        p.value.main<-1-pf(anova(CS.Model)[2,4], anova(CS.Model)[2,1], (dim(GE.data)[1]-anova(CS.Model)[1,1]-anova(CS.Model)[2,1]-2))
        p.values.main<-rbind(p.values.main, p.value.main)
        f<-CS.Model@fixef[length(unique(all.means$ENV))+1]+CS.Model@fixef[1]
        
        fixeff2<-cbind(fixeff2, f)                                                               
      }
    }
    
    if(method=="CIM"){
      new.additive2<-new.additive[,outem[,4]<threshold.f]
      p.values.main<-NULL
      p.values.inter<-NULL
      fixeff2<-NULL
      
      #I will not do it for the reduced set now because of the windows size
      for( i in 1:dim(new.additive)[2]){
        marker<-new.additive[,i]
        if(unlist(cofactor.win.f[i])>0){
          new.list<-cofactor.list[,-c(unlist(cofactor.win.f[i]))]
          number.cofactors<-length(cofactors)-1
        }
        if(unlist(cofactor.win.f[i])==0){ 
          new.list<-cofactor.list
          number.cofactors<-length(cofactors)
        }
        GE.data<-data.frame(all.means, marker)
        
        CS.Model<-lmer(MEAN~ ENV + marker + marker:ENV + new.list + (1|GEN), data=GE.data)
        p.value.inter<-1-pf(anova(CS.Model)[3,4], anova(CS.Model)[3,1], (dim(GE.data)[1]-anova(CS.Model)[1,1]-anova(CS.Model)[2,1]-2))
        
        p.values.inter<-rbind(p.values.inter, p.value.inter)
        p.value.main<-1-pf(anova(CS.Model)[2,4], anova(CS.Model)[2,1], (dim(GE.data)[1]-anova(CS.Model)[1,1]-anova(CS.Model)[2,1]-2))
        
        p.values.main<-rbind(p.values.main, p.value.main)
        f<-CS.Model@fixef[length(unique(all.means$ENV))+1]+CS.Model@fixef[1]
        
        fixeff2<-cbind(fixeff2, f)                                                               
      }
      
      #this because of new.additive instead of new.additive2
      fixeff2<-fixeff2[,outem[,4]<threshold.f]
      p.values.inter<-p.values.inter[outem[,4]<threshold.f,]
      p.values.main<-p.values.main[outem[,4]<threshold.f,]
    }
    
    
    p.values.inter<-data.frame(row.names(b)[outem[,4]<threshold.f], b[outem[,4]<threshold.f,], p.values.inter)
    p.values.main<-data.frame(row.names(b)[outem[,4]<threshold.f], b[outem[,4]<threshold.f,], p.values.main)
    colnames(p.values.main)=c("marker", "Chr", "Pos", "LOD")
    colnames(p.values.inter)=c("marker", "Chr", "Pos", "LOD")
    
  } ####test for interactions if multiple environments
  
  #####To plot profile
  if(max(-log10(outem[,4]),na.rm=TRUE)=="Inf") {max=10} 
  if(max(-log10(outem[,4]),na.rm=TRUE)!="Inf") {max=(max(-log10(outem[,4]))+0.05)} 
  
  #############for heatmaps...
  colnames(fixeff)<-colnames(new.additive)
  rescale<-fixeff[, outem[,4]<threshold.f]
  
  if(length(levels(ENV))>1){
     #pdf(file=paste(QTL.path,"QTLxE heatmap.pdf",sep=""), onefile=TRUE)                        
     png(file=paste(QTL.path,"QTLxE heatmap.png",sep=""))
     par(mfrow = c(1,1))
     heatmap(rescale, Rowv=NA, Colv=NA, scale="none", col=colorRampPalette(c("yellow","blue"))(50)) 
     dev.off()
                             
     #####QTL effect by location and marker plots
     rownames(rescale)<-unique(all.means$ENV)
     
     #pdf(file=paste(QTL.path,"QTLxE effects by location and marker.pdf",sep=""), onefile=TRUE)
     png(file=paste(QTL.path,"QTLxE effects by location and marker.png",sep=""))
     par(mfrow=c(2,1))
     barplot(rescale, beside=TRUE, legend.text=TRUE, xlim=c(0,dim(rescale)[1]*dim(rescale)[2]+30), col=heat.colors(dim(rescale)[1],alpha=1), main="QTL specific effects sorted by Marker and Location")
     barplot(t(rescale), beside=TRUE, legend.text=TRUE, xlim=c(0,dim(rescale)[1]*dim(rescale)[2]+30), col=colorRampPalette(c("dark blue","light blue"))(dim(rescale)[2]))
     dev.off()
  }   

  ######For reporting final model choose only significant markers calculate effects and R squared
  bw.sel<-function(y,X,alpha=0){ #bw function
    XXX<-X
    pval<-rep(alpha,ncol(XXX))
    names(pval)<-colnames(XXX)
    ##
    test=1
    Rsq<-0
    Rvec<-c()
    toss.names<-c()
    while(test==1){
      toss<-which(pval==max(pval)&pval>alpha)
      if(length(pval)>1){toss.name<-names(pval)[toss]
      } else {toss.name<-setdiff(colnames(X),toss.names) }
      toss.names<-c(toss.names,toss.name)
      sel<-setdiff(c(1:ncol(XXX)),toss)
      if (length(sel)==0) {break}
      XXX<-as.matrix(XXX[,sel])	
      CS.Model<-lm(y~XXX);
      r<-summary(CS.Model)$r.squared
      Rvec<-c(Rvec,Rsq-r)
      Rsq<-r
      pval<-summary(CS.Model)$coefficients[,4]; pval <-pval[2:length(pval)]
      names(pval)<-gsub("XXX","",names(pval));
      test<-((sum(pval>=alpha)>0)*1)	
    }
    Rvec<-c(Rvec[2:length(Rvec)],r)
    names(Rvec)<-toss.names
    coef<-coefficients(CS.Model); coef<-coef[2:length(coef)]		
    qtl.names<-names(coef); qtl.names<-gsub("XXX","",qtl.names); 
    
    return(list(qtl.names,Rvec,coef))
  } #bw function
  
  if(nrow(pot.qtl)>0){  
    if(length(levels(ENV))==1) {
      
      #backward select final model
      new.additive.final<-new.additive[,colnames(new.additive)%in%res.qtl$marker]
      BW<-bw.sel(MEAN,new.additive.final,alpha=0.1)
      qtl.names<-BW[[1]]
      res.qtl<-res.qtl[match(qtl.names,res.qtl$marker),]
      m.eff<-BW[[3]]
      
      #backward selection of markers to set R squared
      new.additive.final<-new.additive.final[,match(qtl.names,colnames(new.additive.final))]
      Rsq<-bw.sel(MEAN,new.additive.final,alpha=0)[[2]]
    } 
  }
  
  if(nrow(pot.qtl)==0){  m.eff<-NA; Rsq<-NA;} #output NA if no qtl found
  QTL.result<-NULL
  QTL.result$all<-outem
  QTL.result$selected<-cbind(res.qtl,m.eff,Rsq)
  
  m.eff<-NULL
  
  #convert p tot lod 
  QTL.result$all$'LOD'<-(-log10(QTL.result$all$'LOD'))
  if(nrow(pot.qtl)>0){QTL.result$selected$'LOD'<-(-log10(QTL.result$selected$'LOD'));}
  
  if(length(levels(ENV))>1) {QTL.result$interaction<-p.values.inter} 
  if(length(levels(ENV))==1) {QTL.result$interaction<-NULL}
  if(length(levels(ENV))>1) {QTL.result$main<-p.values.main}
  if(length(levels(ENV))>1) {QTL.result$main<-outem}
  
  row.names(QTL.result$all) <- NULL
  row.names(QTL.result$selected) <- NULL
  
  ####write data to file 
  if (is.null(env.label)==FALSE) {
    fnam = paste(QTL.path,"QTL_", tolower(trait), '_',  env.label, "=", env, "_", tolower(method), '.Rdata', sep = '')
    pngFile = paste(QTL.path,"QTLmap_", trait, '_', env.label, "=", env, "_", method,".png",sep="")
  } else {
    fnam = paste(QTL.path,"QTL_", tolower(trait), '_',  tolower(method), '.Rdata', sep = '')
    pngFile = paste(QTL.path,"QTLmap_", trait, '_', method,".png",sep="")
  }
  
  save(QTL.result, file = fnam)
  #pdf(file = pdfFile,onefile=TRUE);
  png(file = pngFile);
  print(xyplot(-log10(outem[,4])~ outem[,3] | factor(outem[,2]), type="l", layout=c(nchr(file),1), col="red", xlab="Chromosome position", ylab="-log10(P)", main=paste("QTL mapping", method, sep=""), scales = list(x = "free"), ylim=c(0,max), lwd=3,panel = function(x,y,...) {panel.abline(h =-log10(threshold.f),lty=2);llines(x,y,col="red",lwd=2)}));
  
  dev.off()
  
  detach("package:qtl")
  return(QTL.result)
  
} #end of QTL analysis

