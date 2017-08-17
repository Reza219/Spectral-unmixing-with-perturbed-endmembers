function XF=ConvC(X,F,nr,nc,N,K)
%ConvC - defines a circular convolution (the same for all bands) accepting
% a matrix and returnig a matrix. FK is the fft2 of a one-band filter

X=reshape(X',[nr,nc,K]);
F=repmat(F,[1,1,K]);
XF=real(ifft2(fft2(X).*F));
XF=reshape(XF,[N,K])';

end