function [Xhat, A_UE_index,pai_tmp,iter,NMSE_re] = GMMV_AMP(Y, Phi, Pn,damp,niter,tol,K,N_seq,varH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GMMV-AMP algorithm for GMMV CS problem (estimating 3D matrix)
% where incremental EM algorithm is used to learn unknown hyper-parameters.
% An extended version of AMP-NNSPL algorithm in the paper:
% X. Meng et al, "Approximate message passing with nearest neighbor sparsity pattern learning".

% Inputs:
%   Y: received signal
%   Phi: measurement matrix
%   niter: number of AMP iterations
%   tol: termination threshold
%   nns_sel: select nearest neighbor set for NNSPL 0: strctured sparsity 1: clustered sparsity

% Outputs:
%   Xhat: the estimated matrix
%   lambda: belief indicators
%   iter: number of iteration

% Written by Malong Ke (kemalong@bit.edu.cn), Beijing Institute of Technology
% version: 2019.12.03
% xmean, xvar, and nvar are scalars 
% 适用于目标信号的元素独立同分布，即均值方差相等情况下；若均值方差不相等，性能恶化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initializations
[M,Q,P] = size(Y);
[~,N,~] = size(Phi);

% sparsity ratio
alpha = M/N;
normal_cdf = @(x) 1/2.*(1+erf(x/sqrt(2)));
normal_pdf = @(x) 1/sqrt(2*pi).*exp(-x.^2/2);
alpha_grid = linspace(0,10,1024);
rho_SE = (1 - 2/alpha.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)))...
         ./ (1 + alpha_grid.^2 - 2.*((1+alpha_grid.^2).*normal_cdf(-alpha_grid) - alpha_grid.*normal_pdf(alpha_grid)));
lambda = alpha*max(rho_SE).*ones(N,Q,P); 

% noise variance
% snr0 = 100; 
% nvar = norm(Y(:))^2/(snr0+1)/(M*Q*P);
nvar = Pn;
% nvar(ones(1,M),:,:);
% Pn;
% norm(Y(:))^2/(snr0+1)/(M*Q*P);

% mean and variance of target signal
xmean = eps;
Phi_fro = sum(sum(abs(Phi).^2, 1), 2);
Phi_fro = Phi_fro(:, ones(1,Q), :);
xvar = (sum(abs(Y).^2, 1) - M*nvar) ./ Phi_fro;
xvar = sum(xvar(:))/Q/P/(alpha*max(rho_SE));
% xvar=varH;
% other parameters initialization
Xhat = xmean.*ones(N,Q,P); 
v = xvar.*ones(N,Q,P);
V = ones(M,Q,P);
Z = Y;

% allocate memory
D = zeros(N,Q,P);
C = zeros(N,Q,P);
L_cal = zeros(N,Q,P);
pai = zeros(N,Q,P);
A = zeros(N,Q,P);
B = zeros(N,Q,P);
Yhat = zeros(M,Q,P);

% damp=opt.damp;
% niter=opt.niter;
% tol=opt.tol;

%% AMP iteration
for iter = 1:niter
    Xhat_pre = Xhat;
    V_pre = V;
    for p = 1:P
        % factor node update 
        V(:,:,p) = damp.*V_pre(:,:,p) + (1-damp).*abs(Phi(:,:,p)).^2*v(:,:,p);
        Z(:,:,p) = damp.*Z(:,:,p) + (1-damp).*...
                   (Phi(:,:,p)*Xhat(:,:,p) - (Y(:,:,p)-Z(:,:,p))./(nvar+V_pre(:,:,p)).*V(:,:,p));
                                          
        % variable node update 
        D(:,:,p) = 1 ./ ((abs(Phi(:,:,p)).^2).'*(1./(nvar+V(:,:,p))));
        C(:,:,p) = Xhat(:,:,p) + D(:,:,p).*(Phi(:,:,p)'*((Y(:,:,p)-Z(:,:,p))./(nvar+V(:,:,p))));
        
        % compute posterior mean and variance
        L_cal(:,:,p) = (1/2).*(log(D(:,:,p)./(D(:,:,p)+xvar)) + abs(C(:,:,p)).^2./D(:,:,p) - ...
                                   abs(C(:,:,p)-xmean).^2./(D(:,:,p)+xvar)); 
        pai(:,:,p) = lambda(:,:,p) ./ (lambda(:,:,p)+(1-lambda(:,:,p)).*exp(-L_cal(:,:,p))); 
        A(:,:,p) = (xvar.*C(:,:,p)+xmean.*D(:,:,p)) ./ (D(:,:,p)+xvar);
        B(:,:,p) = (xvar.*D(:,:,p)) ./ (xvar+D(:,:,p));
        Xhat(:,:,p) = pai(:,:,p).*A(:,:,p);
        v(:,:,p) = pai(:,:,p).*(abs(A(:,:,p)).^2+B(:,:,p)) - abs(Xhat(:,:,p)).^2;
        
        % reconstruct received signal
        Yhat(:,:,p) = Phi(:,:,p)*Xhat(:,:,p);
    end
    
    % hyper-parameter update based on EM learning
    xmeantemp = mean(sum(pai.*A, 1)./sum(pai, 1));
    xmean = repmat(xmeantemp,N,Q);
    xvartemp = mean(sum(pai.*(abs(xmean-A).^2+B), 1)./sum(pai, 1));
    xvar = repmat(xvartemp,N,Q);
%     nvartemp = mean(sum(abs(Y-Z).^2./abs(1+V./nvar).^2 + V./(1+V./nvar))/M);
%     nvar = repmat(nvartemp,M,Q);
%         xmean = sum(pai.*A, 1)./sum(pai, 1);
%     xmean = sum(xmean(:))/Q/P;
%     xvar = sum(pai.*(abs(xmean-A).^2+B), 1)./sum(pai, 1); 
%     xvar = sum(xvar(:))/Q/P;
%     nvar = sum(abs(Y-Z).^2./abs(1+V./nvar).^2 + V./(1+V./nvar))/M;
%     nvar = sum(nvar(:))/Q/P;
    
    % refine the update rule for sparsity ratio
%     switch nns_sel
%     case 0 % structured sparsity(common support)

        pai_tmp = sum(sum(pai,3),2)./Q./P;
        lambda = pai_tmp(:, ones(1,Q), ones(1,P)); 
%     end  
    
    % stopping critrria
    NMSE_iter =  norm(Xhat(:)-Xhat_pre(:))^2 / norm(Xhat_pre(:))^2;
    NMSE_re = norm(Y(:)-Yhat(:))^2 / norm(Y(:))^2;
    if NMSE_iter < tol || NMSE_re < tol
       break;
    end
    
%     % print status
%     NMSE = norm(Xhat(:)-X(:))^2 / norm(X(:))^2;
%     fprintf('iter = %d, NMSE = %7.9f, NMSE_ietr = %7.9f, NMSE_re = %7.9f\n',...
%              iter, NMSE, NMSE_iter, NMSE_re);
end
temppp=pai_tmp;
temppp(temppp>0.5)=1;
temppp(temppp<0.5)=0;
% 实际利用了结构性
Ind_temp=reshape(temppp,[N_seq,K]);
A_UE_index=sum(Ind_temp);
% fprintf('KML, ');
% fprintf('xmean=%f,',xmeantemp);
% fprintf('xvar=%f,',xvartemp);
% fprintf('|Pn-nvar|=%f\n',abs(Pn-nvartemp));
end