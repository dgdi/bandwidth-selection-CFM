#Data are simulated according to the DGPs described in:
#"Two-stage DEA: caveat emptor"
#Leopold Simar Paul W. Wilson

library(truncnorm)
library(R.matlab)
library(matlab)
library(gtools)
source("Ker_LSCV_OUT.R")

n=100
r=100

for(l in 1:r){
        seed<-sample(1:500000,1)        
        set.seed(seed)
        p=3
        sdX=1
        X<-matrix(nrow=n, ncol=p,rtruncnorm(n*p, a=0, b=1, mean = 0.5, sd = sdX))
        d=2
        sdZ=5
        Z<-matrix(nrow=n, ncol=d,rtruncnorm(n*d, a=0, b=4, mean = 2, sd = sdZ))
        q=3
        inef<-matrix(nrow=n, ncol=q,0)
        meanU=sample(seq(0.25,1,by=0.25), q, replace=TRUE)
        sdU=sample(seq(0.25,2.5,by=0.25), q, replace=TRUE)
        inef<-sapply(1:q,function(x){rtruncnorm(n, a=0, b=2, mean = meanU[x], sd = sdU[x])})
        #We impose this structure of dependence
        #Y1=f(X1,Z1), Y2=f(X2,Z1), Y3=f(X3,Z2)
        #ad similarly
        #Yexp1=f(X1,Z1), Yexp2=f(X2,Z1), Yexp3=f(X3,Z2)
        colx=c(1,2,3)
        colz=c(1,1,2)
        Y<-sapply(1:q,function(h){sqrt(1L-((X[,colx[h]]-1L)^2L))*exp(-(Z[,colz[h]]-2L)^2L) * inef[,h]})
        Yexp<-sapply(1:q,function(h){sqrt(1L-((X[,colx[h]]-1L)^2L))*exp(-(Z[,colz][h]-2L)^2L) * exp(inef[,h])})
       
         #Bandwidths        
        
        h0one=t(c(seq(from=0.3, to=2.1, by=0.3),seq(from=2.5, to=3, by=0.5)))
        hPerm<-permutations(n=9,r=3,v=h0one,repeats.allowed=T)
  
        
        sdY<-apply(Y,2,sd)
        meansdY<-mean(sdY)
        sdYexp<-apply(Yexp,2,sd)
        meansdYexp<-mean(sdYexp)
        sdZ<-apply(Z,2,sd)
        sdTot<-c(meansdY, sdZ)
        sdTotexp<-c(meansdYexp, sdZ)
        h<-apply(hPerm,1,function(d){d*sdTot})
        hexp<-apply(hPerm,1,function(d){d*sdTotexp})
        
        band    <- sapply(1:ncol(h),function(g){sapply(1:n, function(m){Ker_LSCV_OUT(h=h[, g],x=sapply(X[m,], as.numeric),X=X,Y=Y,Z=Z, n=n, q=q, d=d)})})
        bandexp <- sapply(1:ncol(hexp),function(g){sapply(1:n, function(m){Ker_LSCV_OUT(h=hexp[, g], x=sapply(X[m,], as.numeric), X=X, Y=Yexp, Z=Z, n=n, q=q, d=d)})})    
        
        nameMat=paste0("validation_ce_", l,".mat")
        writeMat(con=nameMat, h=h, hexp=hexp, X=X, Y=Y, Yexp=Yexp, Z=Z, n=n, q=q, d=d)
        nameBand=paste0("CV_ce_R_", l,".csv")
        nameBandexp=paste0("CV_ce_exp_R_", l,".csv")
        write.table(seed,file="seeds_.csv", row.names=FALSE, col.names=FALSE, append=TRUE)
        write.table(band,file=nameBand, row.names=FALSE, col.names=FALSE)
        write.table(bandexp,file=nameBandexp, row.names=FALSE, col.names=FALSE)
        print(l)
}