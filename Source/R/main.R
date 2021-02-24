
Fit_Curve<-function(Vecy,inputorders){
  y = reshape(Vecy,direction = "wide")
  MinnumInd = min(apply(apply(y,2,na.omit),2,length))
  nt <- inputorders+1
  nobs=nrow(y)
  coefpheno <- matrix(,nobs,nt)
  if(MinnumInd>nt){
    for(ij in 1:nobs){
      fit <- fastLmPure(y=as.numeric(y[ij,]),X=cbind(1,CLegendre(c(1:ncol(y)),inputorders)))
      coefpheno[ij,]= t(fit$coef)
    }
  }else{
    Vecy=reshape(y, varying = colname,
                 timevar = "time", v.names="trait",
                 direction = 'long')
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
  
  system(paste("./Hi_RRM -bfile ",plinkfilename," -gk 2 -o kinship",sep=""))   
  system(paste("./Hi_RRM -bfile ",plinkfilename,
               " -p mvpheno -k ./output/kinship.sXX.txt -lmm 3 -n ",ntraits," -o vc",sep=""))
  Vc=read.table("vcout")
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
  Ve[upper.tri(Ve)] <- t(Ve)[upper.tri(Vg)]
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
  pv=pchisq(result[,ncol(result)],nt,lower.tail = F)
  result=cbind(result,pv)
  colnames(result)=c(paste("BETA",1:nt,sep=" "),paste("VBETA",1:nt,sep=" "),"chisq","p-value")
  bim=read.table(paste(plinkfilename,".bim",sep=""))
  result=cbind(bim[,c(1,4,2)],result)
  colnames(result)[1:3]=c("CHR","POS","SNP")
  return(result)
}
