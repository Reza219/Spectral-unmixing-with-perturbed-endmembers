function XD=diffcv(X2,nr,nc,N,K)
% circular vertical gradient

X3=reshape(X2',[nr,nc,K]);
XD=diff(X3);
XD=cat(1,XD,X3(1,:,:)-X3(end,:,:));
XD=reshape(XD,[N,K])';

end

