Estimate_coefy<-function(Phenotype,inputorders){
  Vecy = as.matrix(read.table(paste(Phenotype),head=T))
  fixeffect=NULL
  if(ncol(Vecy)!=3){
    ncolVecy=ncol(Vecy)
    fixeffect = (1:ncolVecy)[-c(1,ncolVecy-1,ncolVecy)]
    fix = as.matrix(cbind(Vecy[,1],Vecy[,fixeffect]))
    Vecy= Vecy[,-fixeffect]
    colnames(fix)=c("id",1:length(fixeffect))
    fix=unique(fix)
  }
  Vecy=cbind(Vecy[,3],Vecy[,2],Vecy[,1])
  colnames(Vecy)=c("time","trait","id")
  y = reshape(as.data.frame(Vecy),direction = "wide")
  id = y[,1]
  y = y[,-1]
  MinnumInd = min(apply(apply(y,2,na.omit),2,length))
  nt <- inputorders+1
  nobs=nrow(y)
  coefpheno <- matrix(,nobs,nt)
  if(MinnumInd>nt){
    X = cbind(1,CLegendre(c(1:ncol(y)),inputorders))
    for(ij in 1:nobs){
      fit <- fastLmPure(y=as.numeric(y[ij,]),X=X)
      coefpheno[ij,]= t(fit$coef)
    }  
  }else{
    nobs=length(unique(Vecy$id))
    L <-Legendre (rep(c(1:length(colname)),each=nobs),orders)
    f0 <- "L[,1]"
    for (i in 2:orders){    
      f0 <- paste(f0,"+L[,",i,"]",sep="")   
    }
    f1 <- paste("trait ~ ",f0,sep="")
    f2 <- paste(" ~ 1+",f0,"|id",sep="")   
    f1 <- as.formula(f1)
    f2 <- as.formula(f2)
    control <-list()
    control$maxIter=1000
    control$msMaxIter=1000
    str(lCtr <- lmeControl(maxIter=1000,msMaxIter = 1000,returnObject=TRUE))
    do.call(lmeControl, lCtr)
    modI <- lme(f1, data = Vecy,random = f2,control=lCtr)
    coefpheno = coef(modI)
    id=rownames(coefpheno)
  }
  if(!is.null(fixeffect)){
    coefpheno=cbind(id,coefpheno)
    alldata=merge(coefpheno,fix,by.x="id",sort = F,all.y=F)
    ncolalldata=ncol(alldata)
    fixcol=(ncolalldata-length(fixeffect)+1):ncolalldata
    coefpheno = lm(as.matrix(a[,-fixcol]) ~ as.factor(alldata[,fixcol]))$resi
  }
  return(coefpheno)
}


Hi_RRM<-function(plinkfilename,y){
  y<-as.matrix(y)
  nt=ncol(y)
  nobs=nrow(y)
  Vg<-matrix(,nt,nt)
  Ve<-matrix(,nt,nt)
  fwrite(y,"mvpheno",col.names=F,row.names=F,quote=F)
  g0<-BEDMatrix(plinkfilename)
  nobs=nrow(g0)
  nmar=ncol(g0)
  ntraits <- paste(c(1:nt),collapse=" ")
  
  system(paste("./Hi_RRM_linux -bfile ",plinkfilename," -gk 2 -o kinship",sep=""))   
  system(paste("./Hi_RRM_linux -bfile ",plinkfilename,
               " -p mvpheno -k ./output/kinship.sXX.txt -lmm 3 -n ",ntraits," -o vc",sep=""))
  Vc=as.matrix(read.table("vcout"))
  Vgnull<-Vc[1:(length(Vc)/2)]
  Venull<-Vc[-(1:(length(Vc)/2))]
  count=1
  for (i in 1:nt) {
    for(j in 1:i){
      Vg[i,j]  = as.numeric(Vgnull[count])
      Ve[i,j] = as.numeric(Venull[count])
      count = count+1
    }
  }
  Vg[upper.tri(Vg)] <- t(Vg)[upper.tri(Vg)]
  Ve[upper.tri(Ve)] <- t(Ve)[upper.tri(Ve)]
  G<-fread("./output/kinship.sXX.txt")
  eigen<-eigen(G)
  ug<-eigen$vectors
  sg<-eigen$values
  ynew=t(ug)%*%y
  vl2<-c()
  for(j in 1:nrow(y)){
    deltv= eigen(Vg * sg[j] + Ve)
    vl2<-rbind(vl2,t(deltv$vectors)/sqrt(deltv$values))
  }
  Ynew<-calynew(ynew,vl2)
  blocksize=20000
  zone = floor(nmar/blocksize)+1
  thread = detectCores(,logical=F)
  cl <- makeCluster(thread) 
  zonenumb=c(1:zone)
  result<-matrix(,nmar,2*nt+1) 
  for(ii in zonenumb){
    cat(floor(ii*100/zone),"%")
    st=blocksize*(ii-1)+1
    end=blocksize*ii
    if(ii==zone) end=nmar
    gnew = MatMult(as.matrix(t(ug)),as.matrix(scale(g0[,st:end],center=T,scale=F))) 
    range = ncol(gnew)
    sth = c(0,floor(range/thread)*(1:(thread-1)),range)
    a = parLapply(cl,c(1:thread),fast_mvlmmthreads,Y=Ynew,gnew=gnew,vl=vl2,sth=sth)
    result1<-do.call("rbind",a)
    result[st:end,]=result1[order(result1[,ncol(result1)]),-ncol(result1)]
  }
  stopCluster(cl)
  print("Done!")
  pv=pchisq(result[,ncol(result)],nt,lower.tail = F)
  result=cbind(result,pv)
 
  colnames(result)=c(paste("BETA",1:nt,sep=" "),paste("VBETA",1:nt,sep=" "),"chisq","p-value")
  bim=fread(paste(plinkfilename,".bim",sep=""))
  result=cbind(bim[,c(1,4,2)],result)
  colnames(result)[1:3]=c("CHR","POS","SNP")
  fwrite(result,"HiRRM.txt",col.names=T,row.names=F,quote=F)
  return(result[which(pv<(0.05/nmar)),])
}
