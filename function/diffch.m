function XD=diffch(X,nr,nc,N,K)
% circular horizontal gradient

X=reshape(X',[nr,nc,K]);
XD=diff(X,1,2);
XD=cat(2,XD,X(:,1,:)-X(:,end,:));
XD=reshape(XD,[N,K])';

end