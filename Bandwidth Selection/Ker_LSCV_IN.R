Ker_LSCV_IN <- function(h,y,X,Y,Z,n,p,d){
        
        
        hbx=h[1]
        hx=hbx * apply(X,2,sd)
        wz=rep(1,n)
        Q2y=rep(0,n)
        Q1y=rep(0,n)
        Cst1=(4*pi)^(p/2)*prod(hx)
        Cst=1/Cst1
        Xh=scale(X,center=FALSE, scale=hx)
        attr(Xh, "scaled:scale")<-NULL
        XX=Xh %*% t(Xh)
        DXX=diag(XX)
        DXX1= DXX %x% rep(1,n)
        DXX1= t(matrix(ncol=n, DXX1))
        Convol=Cst*exp(-(1/4) * (DXX1 + t(DXX1) - 2*XX))
        Yv=matrix(rep(1,n))%*%t(as.matrix(y))
        FlagYv=(Y>=Yv)
        FlagY=apply(FlagYv,1, function(x) all(x==TRUE))
        yv=rep(1,n-1)%*%t(as.matrix(y))
        for(i in 1:n){
                zi=Z[i, , drop=FALSE]
                xi=X[i, , drop=FALSE]
                Xi=X[-i, , drop=FALSE]
                Yi=Y[-i, , drop=FALSE]
                Zi=Z[-i, , drop=FALSE]
                flagy=(Yi>=yv)
                flagy1=apply(flagy,1, function(x) all(x==TRUE))
                tempz=(Zi-repmat(as.numeric(zi),n-1,1)) # this is a (n-1) x d matrix
                tempx=(Xi-repmat(as.numeric(xi),n-1,1)) # this is a (n-1) x p matrix
                tempxh=scale(tempx,center=FALSE, scale=hx)
                attr(tempxh, "scaled:scale")<-NULL
                kerxi=(exp(-(tempxh^2)/2)/sqrt(2*pi))%*%diag(x=rep(1, p)/hx, nrow=p) # Individual Gaussian Kernel
                kerx=apply(kerxi,1,prod) # Gaussian product kernel: kery is a (n-1) x 1 vector
                hz=h[c(2:(d+1))]
                tempzh=scale(tempz,center=FALSE, scale=hz)
                attr(tempzh, "scaled:scale")<-NULL
                kerzi=(15/16)*((matrix(nrow=dim(Zi)[1],rep(1, prod(dim(Zi))))-tempzh^2)^2)*(abs(tempzh)<= 1)%*%diag(x = rep(1,d)/hz, nrow=d, ncol=d) # 
                kerz=apply(kerzi,1,prod) # Quartic product kernel: (n- 1) x 1 vector
                mxzi=mean(flagy1*kerz) 
                fyxzi=mean(flagy1*kerz*kerx)
                if(mxzi<=1.0e-06){
                        Q2y[i]=0
                        Q1y[i]=0
                }else{
                        Q2y[i]=fyxzi*wz[i]/mxzi
                        integ=Convol[c(1:n)[-i], c(1:n)[-i]]
                        Bigi1i2=(flagy1 %x% t(flagy1)) * (kerz %x% t(kerz)) * integ
                        Gxzi=mean(Bigi1i2)
                        Q1y[i]=Gxzi*wz[i]/(mxzi^2)
                } 
        }
        I1x=mean(Q1y * FlagY)
        I2x=mean(Q2y * FlagY)
        return( CV=I1x - 2*I2x)
}
