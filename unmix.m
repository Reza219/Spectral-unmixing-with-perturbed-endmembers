clear
addpath('function')

load('C:\Data\MATLAB\HSIData\Botswana\BWs1.mat')
alpha=500;                              %regularization parameter
beta =1e-5;                             %TV regularization parameter
gamma=500;                              %regularization parameter
mu   =1e-3;                             %penalty parameter
lam  =1e-5;                             %penalty parameter
% load('C:\Data\MATLAB\HSIData\IP\IPu1.mat')
% alpha=600;                              %regularization parameter
% beta =1e-5;                             %TV regularization parameter
% gamma=600;                              %regularization parameter
% mu   =1e-3;                             %penalty parameter
% lam  =1e-5;                             %penalty parameter
% load('C:\Data\MATLAB\HSIData\DC\DCs1.mat')
% alpha=500;                              %regularization parameter
% beta=5e-5;                              %TV regularization parameter
% gamma=100;                              %regularization parameter
% mu  =1e-2;                              %penalty parameter
% lam =5e-5;                              %penalty parameter

% paramaters
nt=100;                                 %number of trials
ni=200;                                 %number of iterations
snrx=30;                                %SNR (dB)
snre=30;                                %SNR (dB)
sx=mean(Xg(:).^2)/10^(snrx/10);         %noise variance
se=mean(Eg(:).^2)/10^(snre/10);         %noise variance

dh=zeros(nr,nc); dh(1,1)=-1; dh(1,nc)=1;
dv=zeros(nr,nc); dv(1,1)=-1; dv(nr,1)=1;
IDDT=1./(abs(fft2(dh)).^2+abs(fft2(dv)).^2+1);
W=eye(K)-ones(K)/K;

% initialization
eAa1=0; eAa2=0; eAa3=0; eAa4=0;
eXa1=0; eXa2=0; eXa3=0; eXa4=0;
eaa1=zeros(K,1); eaa2=zeros(K,1); eaa3=zeros(K,1); eaa4=zeros(K,1);
epa1=zeros(N,1); epa2=zeros(N,1); epa3=zeros(N,1); epa4=zeros(N,1);
BAa1=0; BAa2=0; BAa3=0; BAa4=0;
t1=0; t2=0; t3=0; t4=0;

for tt=1:nt
    tt
    % add noise
    Xn=Xg+sqrt(sx)*randn(L,N);
    En=Eg+sqrt(se)*randn(L,K);
    
    tic
    [eA1,A1]=unmix_CLS(En,Xn,Ag,K,N,mu,ni);
    t1=t1+toc; 
    eAa1=eAa1+eA1;
    BAa1=BAa1+(A1-Ag);
    
    tic
    [eA2,A2]=unmix_CTLS(En,Xn,Ag,K,N,mu,ni,alpha);
    t2=t2+toc;
    eAa2=eAa2+eA2;
    BAa2=BAa2+(A2-Ag);
    
    tic
    [eA3,A3]=unmix_CTLS_IV(En,Xn,Ag,K,N,mu,ni,alpha);
    t3=t3+toc;
    eAa3=eAa3+eA3;
    BAa3=BAa3+(A3-Ag);

    tic
    [eA4,A4]=unmix_RCTLS_IV(En,Xn,Ag,K,N,mu,lam,ni,alpha,IDDT,beta,nr,nc,gamma*W);
    t4=t4+toc;
    eAa4=eAa4+eA4;
    BAa4=BAa4+(A4-Ag);
    
    for ii=1:K
        eaa1(ii)=eaa1(ii)+norm(A1(ii,:)-Ag(ii,:))^2;
        eaa2(ii)=eaa2(ii)+norm(A2(ii,:)-Ag(ii,:))^2;
        eaa3(ii)=eaa3(ii)+norm(A3(ii,:)-Ag(ii,:))^2;
        eaa4(ii)=eaa4(ii)+norm(A4(ii,:)-Ag(ii,:))^2;
    end
    for jj=1:N
        epa1(jj)=epa1(jj)+norm(A1(:,jj)-Ag(:,jj))^2;
        epa2(jj)=epa2(jj)+norm(A2(:,jj)-Ag(:,jj))^2;
        epa3(jj)=epa3(jj)+norm(A3(:,jj)-Ag(:,jj))^2;
        epa4(jj)=epa4(jj)+norm(A4(:,jj)-Ag(:,jj))^2;
    end
end

%% results
figure
plot(sqrt(eAa1/nt/K/N),'k:','linewidth',3)
hold on
plot(sqrt(eAa2/nt/K/N),'r--','linewidth',3)
plot(sqrt(eAa3/nt/K/N),'b-.','linewidth',3)
plot(sqrt(eAa4/nt/K/N),'g-','linewidth',3)
legend('CLS','CTLS','CTLS-IV','RCTLS-IV')
xlabel('iteration number')
ylabel('normalized RMSE')

figure
plot(sqrt(eaa1/nt/N),'k:o','linewidth',3)
hold on
plot(sqrt(eaa2/nt/N),'r-.o','linewidth',3)
plot(sqrt(eaa3/nt/N),'b--o','linewidth',3)
plot(sqrt(eaa4/nt/N),'g-o','linewidth',3)
legend('CLS','CTLS','CTLS-IV','RCTLS-IV')
xlabel('endmember number')
ylabel('normalized RMSE')

figure
plot(sort(sqrt(epa1/nt/K)),'k:','linewidth',3)
hold on
plot(sort(sqrt(epa2/nt/K)),'r--','linewidth',3)
plot(sort(sqrt(epa3/nt/K)),'b-.','linewidth',3)
plot(sort(sqrt(epa4/nt/K)),'g-','linewidth',3)
legend('CLS','CTLS','CTLS-IV','RCTLS-IV')
xlabel('pixel number')
ylabel('normalized RMSE')

bias_CLS     =norm(BAa1(:)/nt)/sqrt(K*N)
bias_CTLS    =norm(BAa2(:)/nt)/sqrt(K*N)
bias_CTLS_IV =norm(BAa3(:)/nt)/sqrt(K*N)
bias_RCTLS_IV=norm(BAa4(:)/nt)/sqrt(K*N)

time_CLS     =t1/nt
time_CTLS    =t2/nt
time_CTLS_IV =t3/nt
time_RCTLS_IV=t4/nt