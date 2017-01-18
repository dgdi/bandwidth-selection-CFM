function CV=Ker_LSCV_OUT(h,x,X,Y,Z,n,q,d)
%
% Evaluate the LSCV criterion for estimating a conditional pdf
% of (y |X<=x, Z=z) for baseline bandwidths hby and vector hz
% OUTPUT ORIENTATION
%
hby=h(1);
hy=hby*std(Y); % this is a (1 x q) vector but hby is one-dimensional
wz=ones(n,1);
Q2x=zeros(n,1);
Q1x=zeros(n,1);
Cst1=(4*pi)^(q/2)*prod(hy);
Cst=(1/Cst1) ;
Yh=Y*diag(ones(1,q)./hy);
YY=Yh*Yh';
DYY=diag(YY);
DYY1=kron(DYY,ones(1,n)); %  Kronecker tensor product 
Convol=Cst*exp(-(1/4) * (DYY1 + DYY1' - 2*YY)); 
Xv=ones(n,1)*x;
FlagXv=(X <= Xv);
FlagX=all(FlagXv,2);
xv=ones(n-1,1)*x;
for i=1:n
    zi=Z(i,:);
    yi=Y(i,:);
    Xi=X([(1:i-1)'; (i+1:n)'],:);
    Yi=Y([(1:i-1)'; (i+1:n)'],:);
    Zi=Z([(1:i-1)'; (i+1:n)'],:); 
    flagx=(Xi<=xv);flagx1=all(flagx,2); 
    tempz=(Zi-repmat(zi,n-1,1)); % this is a (n-1) x d matrix
    tempy=(Yi-repmat(yi,n-1,1)); % this is a (n-1) x q matrix
    tempyh=tempy*diag(ones(1,q)./hy); 
    keryi=(exp(-(tempyh .^2)/2)/sqrt(2*pi))*diag (ones(1,q)./hy); % Individual Gaussian Kernel 
    kery=prod(keryi,2); % Gaussian product kernel: kery is a (n-1) x 1 vector
    hz=h(2:d+1);
    tempzh=tempz*diag(ones(1,d)./(hz')); 
    kerzi=(15/16)*((ones(size(Zi))-tempzh .^2) .^2).*(abs(tempzh)<= 1)*diag(ones(1,d)./(hz')); % 
    kerz = prod(kerzi,2); % Quartic product kernel: (n- 1) x 1 vector
    mxzi=mean(flagx1.*kerz); 
    fyxzi=mean(flagx1.*kerz.*kery);
    if(mxzi <=1.0e-06);
       Q2x(i)=0; Q1x(i)=0; 
    else
        Q2x(i)=fyxzi*wz(i)/mxzi; 
    integ=Convol([(1:i-1)'; (i+1:n)'],[(1:i-1)'; (i+1:n)']);
    Bigi1i2=kron(flagx1,flagx1').*kron(kerz,kerz').* integ;
    Gxzi=mean(mean(Bigi1i2)); Q1x(i)=Gxzi*wz(i)/(mxzi^2);
    end
end
I1x=mean(Q1x.*FlagX);
I2x=mean(Q2x.*FlagX); 
CV=I1x - 2*I2x;
end