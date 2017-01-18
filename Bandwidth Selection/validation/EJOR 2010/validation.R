#Data are simulated according to the DGPs described in:
#"Optimal bandwidth selection for conditional efficiency measures: A data-driven approach"
#Luiza Badin, Cinzia Daraio, Léopold Simar
library(matlab)
library(R.matlab)
library(truncnorm)
library(gtools)
source("Ker_LSCV_OUT.R")

n=100
r=100

for(l in 1:r){
        seed<-sample(1:500000,1)        
        set.seed(seed)
        X=matrix(ncol=2, runif(n*2, min=1, max=2))
        y=matrix(ncol=2, runif(n*2, min=0.2, max=5))
        #efficiency frontier Y[i,2] = 1.0845* ((X[i,1])^0.3)* ((x[i,2])^0.4)* - Y[i,1]
        Y<-matrix(nrow=n,ncol=2)
        Y[,1] = ((1.0845*((X[,1])^0.3)) * ((X[,2])^0.4))/((y[,2]/y[,1])+1)
        Y[,2] = ((1.0845*((X[,1])^0.3)) * ((X[,2])^0.4)) - Y[,1]
        #unidimensional Z
        db<-matrix(ncol=13,nrow=n,0)
        db[,1]=X[,1]
        db[,2]=X[,2]
        db[,3]<-Y[,1]
        db[,4]<-Y[,2]
        #U[i] Expo(mean=1/3) for univariate Z
        db[,5]<-rexp(n, rate = 3)
        db[,6]<-db[,3]*exp(-db[,5])
        db[,7]<-db[,4]*exp(-db[,5])
        db[,8:9]<-matrix(nrow=n,ncol=2, runif(n*2, min=1, max=4))
        #U[i] Expo(mean=1/2) for multidimensional Z
        db[,10]<-rexp(n, rate = 2)
        db[,11]<-(1+(2*((abs(db[, 9]-2.5))^3))) * db[,3] * exp(-db[,10])
        db[,12]<-(1+(2*((abs(db[, 9]-2.5))^3))) * db[,4] * exp(-db[,10])
        colnames(db)<-c("x1","x2","y1_","y2_","u1", "y1uni","y2uni","zind", "zrel", "u1mul","u2mul","y1mul","y2mul")
        h0one=t(c(seq(from=0.3, to=2.1, by=0.3),seq(from=2.5, to=3, by=0.5)))
        hPermUni<-permutations(n=9,r=2,v=h0one,repeats.allowed=T)        
        hPermMul<-permutations(n=9,r=3,v=h0one,repeats.allowed=T)
        sdY1uni<-sd(db[,6])
        sdY2uni<-sd(db[,7])
        meansduni<-mean(sdY1uni,sdY2uni)
        sdY1mul<-sd(db[,11])
        sdY2mul<-sd(db[,12])
        meansdmul<-mean(sdY1mul,sdY2mul)
        sdZind<-sd(db[,8])
        sdZrel<-sd(db[,9])
        sdTotUni<-c(meansduni, sdZind)
        sdTotMul<-c(meansdmul, sdZind, sdZrel)
        hUni<-apply(hPermUni,1,function(d){d*sdTotUni})
        hMul<-apply(hPermMul,1,function(d){d*sdTotMul})
        bandUni<- sapply(1:ncol(hUni),function(g){sapply(1:n, function(d){Ker_LSCV_OUT(h=hUni[, g],x=sapply(db[d,1:2], as.numeric),X=db[, 1:2],Y=db[, 3:4],Z=db[, 8, drop=FALSE],n=n,q=2,d=1)})})
        bandMul<- sapply(1:ncol(hMul),function(g){sapply(1:n, function(d){Ker_LSCV_OUT(h=hMul[, g],x=sapply(db[d,1:2], as.numeric),X=db[, 1:2],Y=db[, 11:12],Z=db[, 8:9],n=n,q=2,d=2)})})
        nameBandUni=paste0("CV_ban_Uni_R_", l,".csv")
        nameBandMul=paste0("CV_ban_Mul_R_", l,".csv")
        nameMat=paste0("validation_ban_", l,".mat")
        write.table(seed,file="seeds_ban.csv", row.names=FALSE, col.names=FALSE, append=TRUE)
        writeMat(con=nameMat, hUni=hUni, hMul=hMul, X=db[, 1:2], YUni=db[, 3:4], YMul=db[, 11:12], Z=db[, 8:9], n=n)
        write.table(bandUni,file=nameBandUni, row.names=FALSE, col.names=FALSE)
        write.table(bandMul,file=nameBandMul, row.names=FALSE, col.names=FALSE)
        print(l)
}