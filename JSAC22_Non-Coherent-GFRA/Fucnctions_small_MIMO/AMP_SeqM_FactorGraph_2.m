%% amp-factorGraph+structured for solving media-modulation or Seq-modulation with small BS antennas
% in EM, N_seq=2, can be esaily extended
% Written by Li Qiao, qiaoli@bit.edu.cn, Beijing Institute of Technology
%% Algorithm 1 in the paper
function [Xhat, Afactor,iter,NMSE_re,xvar] = AMP_SeqM_FactorGraph_2(Y, HH,Phi, A_UE_index,Pn,damp,niter,tol,N_seq,sparsity,varH)
% save debug1214.mat Y  HH Phi  A_UE_index Pn damp niter tol N_seq sparsity varH
% AMP 结合 矩阵形式的去噪器 去噪部分的推导类似唐子涵的ACCESS文章
% Initializations
% MatMAX=log(1e200);
% MatMAX=log(1e200);
% MatMIN=1e-20;
[M,Q] = size(Y);
[~,N,~] = size(Phi);
alpha = M/N;
normal_cdf = @(x) 1/2.*(1+erf(x/sqrt(2)));
normal_pdf = @(x) 1/sqrt(2*pi).*exp(-x.^2/2);
alpha_grid = linspace(0,10,1024);
rho_SE = (1 - 2/alpha.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)))...
    ./ (1 + alpha_grid.^2 - 2.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)));
Afactor = alpha*max(rho_SE).*ones(N/N_seq,1);
nvar = Pn;
% snr0 = 100; 
% nvar = norm(Y(:))^2/(snr0+1)/(M*Q);
% mean and variance of target signal
xmean = eps;
Phi_fro = sum(sum(abs(Phi).^2, 1), 2);
Phi_fro = Phi_fro(:, ones(1,Q), :);
xvar = (sum(abs(Y).^2, 1) - M*nvar) ./ Phi_fro;
xvar = sum(xvar(:))/Q/(alpha*max(rho_SE));
% xvar=varH;
% other parameters initialization
Xhat = zeros(N,Q);
v = ones(N,Q);
V = ones(M,Q);
Z = Y;
NMSE=zeros(300,1);
%% AMP iteration
for iter = 1:200
    Xhat_pre = Xhat;
    V_pre = V;
    % factor node update
    V= damp.*V_pre + (1-damp).*abs(Phi).^2*v;
    Z = damp.*Z + (1-damp).*...
        (Phi*Xhat - (Y-Z)./(nvar+V_pre).*V);
    % variable node update
    Vi = 1 ./ ((abs(Phi).^2).'*(1./(nvar+V)));
    Ri = Xhat + Vi.*(Phi'*((Y-Z)./(nvar+V)));
    
    L_0= log(Vi./(Vi+xvar)) + abs(Ri).^2./Vi - abs(Ri-xmean).^2./(Vi+xvar);
    ExpL0=exp(L_0);
    SumExpL=zeros(N/N_seq,Q);
    ProdSEL=zeros(N/N_seq,1);
    for q=1:Q
       SumExpL(:,q)=sum(reshape(ExpL0(:,q),N_seq,N/N_seq),1); 
    end
    R_SumExpL=1./(1e-20+SumExpL);
    ProdSEL=prod(R_SumExpL,2);
    SUM=zeros(N,Q);
    PROD=repmat(reshape(repmat(ProdSEL.',N_seq,1),N,1),1,Q);
    for q=1:Q
        SUM(:,q)=reshape(repmat(SumExpL(:,q).',N_seq,1),N,1);
    end
    Lambda=repmat(reshape(repmat(Afactor.',N_seq,1),N,1),1,Q);
    fenmu=Lambda+N_seq^M*(1-Lambda).*PROD;
    PAI=max(1e-20,min(1-1e-20,(1e-20+Lambda.*ExpL0)./(1e-20+SUM.*fenmu)));
    PAI0=PAI./(1-PAI);
    PAI01=zeros(N/N_seq,Q);
    for q=1:Q
        PAI01(:,q)=sum(reshape(PAI0(:,q),N_seq,N/N_seq),1); 
    end
%         PAI1=zeros(Q,N/N_seq);
%     PAI2=zeros(Q,N/N_seq);
%     for q=1:Q
%         TEMMM=reshape(PAI(:,q),N_seq,N/N_seq);
%         PAI1(q,:)=TEMMM(1,:);
%         PAI2(q,:)=TEMMM(2,:);
%     end
%     THETA=exp(-abs(Ri).^2./Vi)/pi./Vi;
%     GAMMA=exp(-abs(Ri-xmean).^2./(xvar+Vi))/pi./(xvar+Vi);
%     THprod=zeros(N/N_seq,Q);
%     TH1GA2=zeros(N/N_seq,Q);
%     TH2GA1=zeros(N/N_seq,Q);
%     for q=1:Q
%         tempT=THETA(:,q);
%         tempG=GAMMA(:,q);
%         tempT2=reshape(tempT,N_seq,N/N_seq);
%         THprod(:,q)=prod(tempT2,1);
%         tempG2=reshape(tempG,N_seq,N/N_seq);
%         TH1GA2(:,q)=tempT2(1,:).*tempG2(2,:);
%         TH2GA1(:,q)=tempG2(1,:).*tempT2(2,:);
%     end
%     NorZ=(1-Afactor).*prod(THprod,2)+Afactor.*prod((TH1GA2+TH2GA1)/2,2);
%     PAI1=zeros(Q,N/N_seq);
%     PAI2=zeros(Q,N/N_seq);
%     PAI=zeros(N,Q);
%     for q=1:Q
%        PAI1(q,:)=Afactor.*prod((TH1GA2+TH2GA1)/2,2).*TH2GA1(:,q)./((TH1GA2(:,q)+TH2GA1(:,q))/2)/2./NorZ;
%        PAI2(q,:)=Afactor.*prod((TH1GA2+TH2GA1)/2,2).*TH1GA2(:,q)./((TH1GA2(:,q)+TH2GA1(:,q))/2)/2./NorZ;
%        PAI(:,q)=reshape([PAI1(q,:);PAI2(q,:)],N,1);
%     end
    meanBar=(Vi.*xmean+Ri.*xvar)./(xvar+Vi);
    varBar=xvar.*Vi./(xvar+Vi);
    meanBar1=zeros(N/N_seq,Q);
    meanBar2=zeros(N/N_seq,Q);
    varBar1=zeros(N/N_seq,Q);
    varBar2=zeros(N/N_seq,Q);
    for q=1:Q
        tempMB=reshape(meanBar(:,q),N_seq,N/N_seq);
        tempVB=reshape(varBar(:,q),N_seq,N/N_seq);
        meanBar1(:,q)=tempMB(1,:);
        meanBar2(:,q)=tempMB(2,:);
        varBar1(:,q)=tempVB(1,:);
        varBar2(:,q)=tempVB(2,:);
    end
    % posterior mean
    Xhat = PAI.*meanBar; 
    % posterior variance
    v = PAI.*(abs(meanBar).^2+varBar) - abs(Xhat).^2;
    Yhat = Phi*Xhat;
    %% EM update
%     % Noise variance
%     nvar = mean(sum(abs(Y-Z).^2./abs(1+V./nvar).^2 + V./(1+V./nvar))/M);
%     % Gaussian prior mean
%----------------------------------------------------------%
TTTTPPP=PAI.*meanBar;
TTTPPP2=PAI.*(abs(xmean-meanBar).^2+varBar);
xmean = sum(TTTTPPP(:))./sum(PAI(:));
    xvar = sum(TTTPPP2(:))./sum(PAI(:)); 
    %---------------------------------------------------------%
% xmean = mean(sum(PAI.*meanBar, 1)./sum(PAI, 1));
%     xvar = mean(sum(PAI.*(abs(xmean-meanBar).^2+varBar), 1)./sum(PAI, 1)); 
%     numerM=(1-PAI1.').*PAI2.'.*meanBar2+(1-PAI2.').*PAI1.'.*meanBar1;
%     DenoM=(1-PAI1.').*PAI2.'+(1-PAI2.').*PAI1.';
%     xmean=sum(numerM(:))./sum(DenoM(:));
% % %     % Gaussian prior variance
%     numerV=(1-PAI1.').*PAI2.'.*(abs(xmean-meanBar2).^2+varBar2)+(1-PAI2.').*PAI1.'.*(abs(xmean-meanBar1).^2+varBar1);
%     DenoV=DenoM;
%     xvar=sum(numerV(:))./sum(DenoV(:));
    % activity indicators
%     Afactor=max(1e-20,min(1-1e-20,mean(DenoM./(1-PAI1.'.*PAI2.'),2)));%只适用于N_seq=2;
Afactor=max(1e-20,min(1-1e-20,mean(1./(1+1./PAI01),2)));

    
    % stopping critrria
    NMSE_iter =  norm(Xhat(:)-Xhat_pre(:))^2 / norm(Xhat_pre(:))^2;
    NMSE_re = norm(Y(:)-Yhat(:))^2 / norm(Y(:))^2;
    NMSE(iter)=(norm((Xhat-HH),'fro'))^2/(norm(HH,'fro'))^2;
    if NMSE_iter < tol || NMSE_re < tol
        break;
    end
end
end