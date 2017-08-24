function [eA2,A,E] = unmix_CTLS_IV(En,Xn,Ag,K,N,mu,ni,alpha)

% initialization
A=zeros(K,N);
V=A;
G=0;
eA2=zeros(ni,1);

% BCD iterations
for ii=1:ni
    % wrt E
    E=(Xn*A'+alpha*En)/(A*A'+alpha*eye(K));
    
    % wrt A
    A=(E'*En+mu*eye(K))\(E'*Xn+mu*(V+G));
    
    % wrt V
    Nu=A-G;
    V=Nu;
    V(V<0)=0;
    V(V>1)=1;
    
    % update G
    G=V-Nu;
    
    % square abundance error
    eA2(ii)=norm(A(:)-Ag(:))^2;
end

end