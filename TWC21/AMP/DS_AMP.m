function [LLR_MBM,LLR_Constell,x_hat,alpha_new]=DS_AMP(y,A,opt,X)
[Nr,N]=size(A);
K=opt.K; % number of users
Nt=N/K; % number of antennas per user
% sigma2=opt.sigma2; % noise variance of channel
alphabet_set=opt.Set; % the modulation symbol set
M = length(alphabet_set);
Iter_num=opt.max_iter;
sigma2=100;
damp=opt.damp;
% Initialization
T=opt.T;
% T=opt.T_coded;
y=y.';
Z = y; %
V = ones(T,Nr);
alpha=0.5*ones(1,K);
alpha_repmat=reshape(repmat(alpha,Nt,1),[1,K*Nt]);
% beta=1/Nt*ones(1,Nt*K);
x_hat=alpha_repmat.*sum(alphabet_set)/M.* ones(T,N);
var_hat=alpha_repmat.*sum(abs(alphabet_set).^2)/M/Nt.* ones(T,N)-x_hat.^2 ;
t=1;
spars=opt.sparsity;
% pluss=Nt*M*(1/spars-1);
while (1)
    
    V_new = (abs(A).^2 * var_hat.').';
    Z_new = (A * x_hat.').' - ((y - Z) ./ (sigma2 + V) ) .* V_new;
    
    Z_new=damp*Z+(1-damp)*Z_new;
    V_new=damp*V+(1-damp)*V_new;
    
    var_1 = (1 ./ (sigma2 + V_new) ) * (abs(A).^2);
    var_2 = ((y - Z_new) ./ (sigma2 + V_new) ) * conj(A);
    
    Ri = var_2 ./ var_1 + x_hat;
    Vi = 1 ./ var_1;
    
    % update the noise variance
    TEMPPPP=(abs(y-Z_new).^2)./abs(1+V_new/sigma2).^2+sigma2*V_new./(V_new+sigma2);
    sigma2_new=mean(TEMPPPP(:));
    sigma2=sigma2_new;
    %%   Poser mean and variance computation
    x_temp=[0,alphabet_set];
    ABSx_temp=abs(x_temp);
    X_temp=repmat(x_temp.',1,K*Nt);
    pf8=zeros(M+1,K*Nt,T);
    for j=1:T
        for i=1:M+1
            pf8(i,:,j)=exp(-abs(Ri(j,:)-X_temp(i,:)).^2./(Vi(j,:)))./Vi(j,:)/pi;
        end
    end
    pf7=zeros(M+1,K*Nt,T);
    for j=1:T
        for i=1:M+1
            if x_temp(i)==0
                pf7(i,:,j)=pf8(i,:,j).*(ones(1,K*Nt)-alpha_repmat/Nt);
            else
                pf7(i,:,j)=pf8(i,:,j).*(alpha_repmat/M/Nt);
            end
        end
    end
    PF7=sum(pf7);
    pf6=zeros(M+1,K*Nt,T);
    for j=1:T
        for i=1:M+1
            pf6(i,:,j)=pf7(i,:,j)./PF7(:,:,j);
        end
    end
    %% Calculate LLR
    [LLR_MBM]=LLR_soft_MBM1203(pf6,opt);
    [LLR_Constell]=Constell_soft_LLR(pf6,opt);
    for j=1:T
        x_hat_new(j,:)=x_temp*pf6(:,:,j);
        var_hat_new(j,:)=(ABSx_temp.^2)*pf6(:,:,j)-(abs(x_hat_new(j,:))).^2;
    end
    %% Update the activity indicator
    alpha_new_T=zeros(T,K);
    for PPP=1:T
        for i=1:K
            TEMP=pf6(:,(i-1)*Nt+1:i*Nt,PPP);
            temp=1;
            for j=1:Nt
                temp=temp*TEMP(1,j);
            end
            TEMP2=zeros(1,Nt);
            for j=1:Nt
                TEMP2(j)=temp/TEMP(1,j);
                alpha_new_T(PPP,i)=alpha_new_T(PPP,i)+sum(TEMP2(j)*TEMP(2:M+1,j));
            end
        end
    end
    alpha_new=mean(alpha_new_T,1);
    alpha_repmat_new=reshape(repmat(alpha_new,Nt,1),[1,K*Nt]);
    if t>Iter_num
        break;
    end
    x_hat=x_hat_new;
    var_hat=var_hat_new;
    alpha_repmat=alpha_repmat_new;
    
    V = V_new;
    Z = Z_new;
    t=t+1;
end






