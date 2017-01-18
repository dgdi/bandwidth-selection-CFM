
Ker_LSCV_OUT <- function(h,x,X,Y,Z,n,q,d){
        
        
        hby=h[1]
        hy=hby * apply(Y,2,sd)
        wz=rep(1,n)
        Q2x=rep(0,n)
        Q1x=rep(0,n)
        Cst1=(4*pi)^(q/2)*prod(hy)
        Cst=1/Cst1
        Yh=scale(Y,center=FALSE, scale=hy)
        attr(Yh, "scaled:scale")<-NULL
        YY=Yh %*% t(Yh)
        DYY=diag(YY)
        DYY1= DYY %x% rep(1,n)
        DYY1= t(matrix(ncol=n, DYY1))
        Convol=Cst*exp(-(1/4) * (DYY1 + t(DYY1) - 2*YY))
        Xv=matrix(rep(1,n))%*%t(as.matrix(x))
        FlagXv=(X<=Xv)
        FlagX=apply(FlagXv,1, function(x) all(x==TRUE))
        xv=rep(1,n-1)%*%t(as.matrix(x))
        for(i in 1:n){
                zi=Z[i, , drop=FALSE]
                yi=Y[i, , drop=FALSE]
                Xi=X[-i, , drop=FALSE]
                Yi=Y[-i, , drop=FALSE]
                Zi=Z[-i, , drop=FALSE]
                flagx=(Xi<=xv)
                flagx1=apply(flagx,1, function(x) all(x==TRUE))
                tempz=(Zi-repmat(as.numeric(zi),n-1,1)) # this is a (n-1) x d matrix
                tempy=(Yi-repmat(as.numeric(yi),n-1,1)) # this is a (n-1) x q matrix
                tempyh=scale(tempy,center=FALSE, scale=hy)
                attr(tempyh, "scaled:scale")<-NULL
                keryi=(exp(-(tempyh^2)/2)/sqrt(2*pi))%*%diag(x=rep(1, q)/hy, nrow=q) # Individual Gaussian Kernel
                kery=apply(keryi,1,prod) # Gaussian product kernel: kery is a (n-1) x 1 vector
                hz=h[c(2:(d+1))]
                tempzh=scale(tempz,center=FALSE, scale=hz)
                attr(tempzh, "scaled:scale")<-NULL
                kerzi=(15/16)*((matrix(nrow=dim(Zi)[1],rep(1, prod(dim(Zi))))-tempzh^2)^2)*(abs(tempzh)<= 1)%*%diag(x = rep(1,d)/hz, nrow=d, ncol=d) # 
                kerz=apply(kerzi,1,prod) # Quartic product kernel: (n- 1) x 1 vector
                mxzi=mean(flagx1*kerz) 
                fyxzi=mean(flagx1*kerz*kery)
                if(mxzi<=1.0e-06){
                        Q2x[i]=0
                        Q1x[i]=0
                }else{
                        Q2x[i]=fyxzi*wz[i]/mxzi
                        integ=Convol[c(1:n)[-i], c(1:n)[-i]]
                        Bigi1i2=(flagx1 %x% t(flagx1)) * (kerz %x% t(kerz)) * integ
                        Gxzi=mean(Bigi1i2)
                        Q1x[i]=Gxzi*wz[i]/(mxzi^2)
                } 
        }
        I1x=mean(Q1x * FlagX)
        I2x=mean(Q2x * FlagX)
        return( CV=I1x - 2*I2x)
}
