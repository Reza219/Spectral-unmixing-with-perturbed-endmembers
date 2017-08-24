function [eA2,A] = unmix_CLS(En,Xn,Ag,K,N,mu,ni)
% CLS-ADMM

% initialization
%A=zeros(Ne,Np);
V=zeros(K,N);
G=0;
eA2=zeros(ni,1);

% constants
MAI=(En'*En+mu*eye(K))^-1;
ETX=En'*Xn;

% ADMM iterations
for ii=1:ni
    % wrt A
    A=MAI*(ETX+mu*(V+G));
    
    % wrt V
    V=A-G;
    V(V<0)=0;
    V(V>1)=1;
    
    % update G
    G=G-(A-V);
    
    % square abundance error
    eA2(ii)=norm(A(:)-Ag(:))^2;
end

end

