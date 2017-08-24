function [eA2,A] = unmix_RCTLS_IV(En,Xn,Ag,K,N,mu,lam,ni,alpha,IDDT,beta,nr,nc,gammaW)

% initialization
A=zeros(K,N);
V1=A; V2h=A; V2v=A; V3=A;
G1=0; G2h=0; G2v=0; G3=0;
eA2=zeros(ni,1);

% BCD iterations
for ii=1:ni
    % wrt E
    E=(Xn*A'+alpha*En)/(A*A'+alpha*eye(K)+gammaW);
    
    % wrt A
    A=(E'*En+(lam+mu)*eye(K))\(E'*Xn+lam*(V1+G1)+mu*(V3+G3));
    
    % wrt V1
    N1=A-G1;
    V1=ConvC(N1+diffchT(V2h+G2h,nr,nc,N,K)+diffcvT(V2v+G2v,nr,nc,N,K),IDDT,nr,nc,N,K);
    
    % wrt V2
    Nh=diffch(V1,nr,nc,N,K)-G2h;
    Nv=diffcv(V1,nr,nc,N,K)-G2v;
    [V2h,V2v]=vec_soft_col_iso(Nh,Nv,beta/lam);
    
    % wrt V3
    N3=A-G3;
    V3=N3;
    V3(V3<0)=0;
    V3(V3>1)=1;
    
    % update Gs
    G1 =V1-N1;
    G2h=V2h-Nh;
    G2v=V2v-Nv;
    G3 =V3-N3;
    
    % square abundance error
    eA2(ii)=norm(A(:)-Ag(:))^2;
end

end