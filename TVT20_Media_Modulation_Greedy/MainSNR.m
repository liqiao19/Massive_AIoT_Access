%% Main Function: BER Pe versus SNR
%"Compressive Sensing Based Massive Access for IoT Relying on Media Modulation Aided Machine Type Communications"
%Accepted by IEEE Transactions on vehicular technology in 2020/06/28
% Written by Li Qiao (qiaoli@bit.edu.cn) from BIT, in 2020/07/28
clear;close all;
%% Basic parameters
opt.nr=50;           % Number of antennas at the BS
opt.K=100;           % Number of UEs
opt.S=8;             % Number of active UEs
opt.mrf=2;           % Bits carriered by index modulated signals
opt.M=2^opt.mrf;     % Dimension of index modulated signals
opt.M_order=4;       % Dimensions of traditional modulation schemes    e.g., 4QAM 16QAM 64QAM; 4PSK, 8PSK,16PSK,64PSK;
opt.Normuliza=sqrt(3/2/(opt.M_order-1));
opt.sparsity=opt.S/opt.K;

opt.Na=1;            % Number of active antennas for spatial modulation
opt.Max_iter=10;     % Iteration numbers of GSP algorithm
opt.QAM_en=1;        % Using QAM
opt.PSK_en=0;        % Using PSK
opt.T=12;            % Length of a frame
opt.awgn_en=1;       % Adding AWGN noise
sim_num=1e8;         % Maximum iterations
% simuamount=5000;   % Actual simulation times
opt.snr_db_set=0:2:10;% SNR
opt.snr_num=10.^(opt.snr_db_set/10);
LLL=length(opt.snr_db_set);
Pe_BSAMP=zeros(LLL,1);
Pe_AMP=zeros(LLL,1);
Pe_TLSSCS=zeros(LLL,1);
Pe_SK_BSAMP=zeros(LLL,1);
BER_SICSSP=zeros(LLL,1);
BER_TLSSCS=zeros(LLL,1);
BER_AMP=zeros(LLL,1);
BER_OLS=zeros(LLL,1);
BER_OLSWMM=zeros(LLL,1);
BER_GSP=zeros(LLL,1);
BER_SAGETLSSCS=zeros(LLL,1);
BER_SAGEBSAMP=zeros(LLL,1);
%% Initialization
index_mr_single=zeros(opt.S,opt.T); 
index_mr_aggregate=zeros(opt.S,opt.T);
opt.tx_data=zeros(opt.S,opt.T);
tx_bin=zeros(opt.S,log2(opt.M_order),opt.T);
tx_signal=zeros(opt.S,opt.T);
%% First loop for SNR
for kk=1:length(opt.snr_db_set) 
    y=zeros(opt.nr,opt.T);
    y2=zeros(opt.nr,opt.T);
    err_BSAMP=0;
    err_AMP=0;
    err_TLSSCS=0;
    err_SK_BSAMP=0;
    ErrBit_SICSSP=0;
    ErrBit_TLSSCS=0;
    ErrBit_AMP=0;
    ErrBit_OLS=0;
    ErrBit_OLSWMM=0;
    ErrBit_GSP=0;
    ErrBit_SAGETLSSCS=0;
    ErrBit_SAGEBSAMP=0;
    if opt.snr_db_set(kk)<=2
        simuamount=3000;%break
    elseif opt.nr>2
        simuamount=6000;%break
    end
    %% Second loop for thousands of simulations
    for i2222=1:sim_num  
        % Generate the channel H
        H=(randn(opt.nr,opt.K*opt.M)+1i*randn(opt.nr,opt.K*opt.M))/sqrt(2);
        X=zeros(opt.K*opt.M,opt.T);
        % Generate the signal x
        temp_=(0:opt.K-1)*opt.M;
        opt.index_user_single=randperm(opt.K,opt.S);
        opt.index_user_single=sort(opt.index_user_single);
        index_user_single_ascend=sort(opt.index_user_single,'ascend');
        index_user_aggregate=(opt.index_user_single-1)*opt.M+1;
        
        for jj=1:opt.T
            opt.index_mr_single(:,jj)=randi([1,opt.M],1,opt.S)-1;
            index_mr_aggregate(:,jj)=index_user_aggregate'+opt.index_mr_single(:,jj);
            opt.tx_data(:,jj) = randi(opt.M_order,[opt.S,1])-1;
            tx_bin(:,:,jj)=de2bi(opt.tx_data(:,jj),log2(opt.M_order),'left-msb');
           
            if opt.QAM_en==1 % Using QAM
                tx_signal(:,jj) = qammod(opt.tx_data(:,jj),opt.M_order)*opt.Normuliza;
                opt.Set=qammod(0:opt.M_order-1,opt.M_order)*opt.Normuliza;
            elseif opt.PSK_en==1% Using PSK
                tx_signal(:,jj) = pskmod(opt.tx_data(:,jj),opt.M_order);
                opt.Set=pskmod(0:opt.M_order-1,opt.M_order);
            end
            X(index_mr_aggregate(:,jj),jj)=tx_signal(:,jj);
        end
        P_n=zeros(opt.T,1);
        if opt.awgn_en==1
            for pp=1:opt.T
                y0=H*X(:,pp);
                P_s=y0'*y0/length(y0);
                P_n(pp)=P_s/opt.snr_num(kk);
                y(:,pp)=H*X(:,pp)+sqrt(P_n(pp)/2).*(randn(opt.nr,1)+1i*randn(opt.nr,1));
            end
        else
            y=H*X;
        end
        opt.sigma2=mean(P_n);
        %%  AUD
        % StrOMP
        [support_BSAMP,NORM]=StrOMP(y,H,opt);
        % AUD lower bound: StrOMP with known number of active UEs
        [support_SK_BSAMP]=SK_StrOMP(y,H,opt);
        % TLSSCS
        sigma2=0.5*opt.S./opt.snr_num(kk);
        imax=opt.S+1;
        betaadaptive=4;
        Vthnormcsmp=betaadaptive*4*sigma2*opt.T; %Terminating threshold
        covr=y*y'/opt.T;
        [eigspace,eigvalue,v]=svd(covr);
        orthspany=eigspace(:,1:opt.T);
        [support_TLSSCS,EstSig_TLSSCS,mbm_TLSSCS] = TLSSCS(orthspany,y,H,imax,opt.M,opt.T,Vthnormcsmp);
        %%  Pe calculation
        user_num=zeros(1,opt.K);
        user_num(opt.index_user_single)=1;
        
        user_num_detected1=zeros(1,opt.K);
        user_num_detected1(support_BSAMP)=1;
        gapp1=user_num-user_num_detected1;
        err_BSAMP=err_BSAMP+nnz(gapp1);
        
        user_num_detected2=zeros(1,opt.K);
        user_num_detected2(support_TLSSCS)=1;
        gapp2=user_num-user_num_detected2;
        err_TLSSCS=err_TLSSCS+nnz(gapp2);
        
        user_num_detected3=zeros(1,opt.K);
        user_num_detected3(support_SK_BSAMP)=1;
        gapp3=user_num-user_num_detected3;
        err_SK_BSAMP=err_SK_BSAMP+nnz(gapp3);
        
        Pe_BSAMP(kk)=err_BSAMP/(i2222*opt.K);
        Pe_TLSSCS(kk)=err_TLSSCS/(i2222*opt.K);
        Pe_SK_BSAMP(kk)=err_SK_BSAMP/(i2222*opt.K);
        Pe_AMP(kk)=err_AMP/(i2222*opt.K);
        
        fprintf('SNR=%d,',opt.snr_db_set(kk));
        fprintf('sim_num=%d,',i2222);
        fprintf('nr=%d,',opt.nr);
        fprintf('T=%d,',opt.T);
        fprintf('S=%d,',opt.S);
        fprintf('Pe_BSAMP=%f,',Pe_BSAMP(kk));
        fprintf('Pe_TLSSCS=%f,',Pe_TLSSCS(kk));
        fprintf('Pe_SK_BSAMP=%f   ',Pe_SK_BSAMP(kk));
        %%  data decoding
        % SIC-SSP
        if isempty(support_BSAMP)==1
            mbm_SICSSP=0;
            EstSig_SICSSP=0;
        else
            [mbm_SICSSP,EstSig_SICSSP,OOO]=SICSSP(support_BSAMP,y,H,opt);
        end
        % TLSSCS £¨LMMSE£©
        
        % oracle-LS
        EstSig_OLS=zeros(opt.S,opt.T);
        for ppp=1:opt.T
            A=H(:,index_mr_aggregate(:,ppp));
            EstSig_OLS(:,ppp)=(A'*A)^(-1)*A'*y(:,ppp);
        end
%         % oracle-LS without media modulation (OLSWMM)

        % GSP
        if isempty(support_BSAMP)==1
            mbm_GSP=0;
            EstSig_GSP=0;
        else
            support_num=length(support_BSAMP);
            mbm_GSP=zeros(support_num,opt.T);
            EstSig_GSP=zeros(support_num,opt.T);
            support_est_M=zeros(1,support_num*opt.M);
            for qqqq=1:support_num
                support_est_M(((qqqq-1)*opt.M+1):qqqq*opt.M)=((support_BSAMP(qqqq)-1)*opt.M+1):(support_BSAMP(qqqq)*opt.M);
            end
            for lll=1:opt.T
                [aaa,support_spatial ,used_ietr]=SP_modify_std_DCS_plus_ql_v1(y(:,lll),H(:,support_est_M),support_BSAMP,opt);
                mbm_GSP(:,lll)=support_spatial;
                EstSig_GSP(:,lll)=aaa;
            end
        end
        %% BER calculation
        
        % TLSSCS
        [ErrBit_UE1,ErrBit_mbm1,ErrBit_sig1]=BERcal('TLSSCS',mbm_TLSSCS,EstSig_TLSSCS,tx_bin,support_TLSSCS,opt);
        ErrBit_TLSSCS=ErrBit_TLSSCS+ErrBit_sig1+ErrBit_mbm1+ErrBit_UE1;
        
        % Proposed
        [ErrBit_UE2,ErrBit_mbm2,ErrBit_sig2]=BERcal('BSAMP-SICSSP',mbm_SICSSP,EstSig_SICSSP,tx_bin,support_BSAMP,opt);
        ErrBit_SICSSP=ErrBit_SICSSP+ErrBit_sig2+ErrBit_mbm2+ErrBit_UE2;
        
        % BER lower bound
        mbm_OLS=0;support_OLS=0;
        [ErrBit_UE3,ErrBit_mbm3,ErrBit_sig3]=BERcal('OracleLS',mbm_OLS,EstSig_OLS,tx_bin,support_OLS,opt);
        ErrBit_OLS=ErrBit_OLS+ErrBit_sig3;
        
        % StrOMP + GSP
        [ErrBit_UE4,ErrBit_mbm4,ErrBit_sig4]=BERcal('BSAMP-GSP',mbm_GSP,EstSig_GSP,tx_bin,support_BSAMP,opt);
        ErrBit_GSP=ErrBit_GSP+ErrBit_sig4+ErrBit_mbm4+ErrBit_UE4;
        
        % benchmark 1
        
        BER_SICSSP(kk)=ErrBit_SICSSP/i2222/opt.S/opt.T/(log2(opt.M)+log2(opt.M_order));
        BER_TLSSCS(kk)=ErrBit_TLSSCS/i2222/opt.S/opt.T/(log2(opt.M)+log2(opt.M_order));
        BER_OLS(kk)=ErrBit_OLS/i2222/opt.S/opt.T/(log2(opt.M)+log2(opt.M_order));
        BER_GSP(kk)=ErrBit_GSP/i2222/opt.S/opt.T/(log2(opt.M)+log2(opt.M_order));
        
        fprintf('BER_SICSSP=%f,',BER_SICSSP(kk));
        fprintf('BER_TLSSCS=%f,',BER_TLSSCS(kk));
        fprintf('BER_GSP=%f,',BER_GSP(kk));
        fprintf('BER_OLS=%f\n',BER_OLS(kk));
%         fprintf('BER_OLSWMM=%f,\n',BER_OLSWMM(kk));
        if i2222>=simuamount
            break;
        end
    end
    save('PeCompareSNRchange.mat','Pe_BSAMP','Pe_TLSSCS','Pe_SK_BSAMP','BER_SICSSP','BER_TLSSCS','BER_GSP','BER_OLS','BER_OLSWMM');
end
%% Plot
figure;
semilogy(opt.snr_db_set, Pe_TLSSCS,'b-s');hold on;
semilogy(opt.snr_db_set, Pe_BSAMP,'r-s');hold on;
semilogy(opt.snr_db_set, Pe_SK_BSAMP,'k-s');hold on;

grid on;
legend('Pe-TLSSCS','Pe-BSAMP','Pe-SK-BSAMP');
xlabel('SNR'), ylabel('$P_e$','interpreter','latex');
title([' nr= '  ,num2str(opt.nr),' T= '  ,num2str(opt.T),' S= ' ,num2str(opt.S),  ' M-order= ',  num2str(opt.M_order),  ' QAM-en= ' ,num2str(opt.QAM_en),  ' PSK-en= ',num2str(opt.PSK_en)]);

figure;
semilogy(opt.snr_db_set, BER_OLSWMM,'b-+');hold on;
semilogy(opt.snr_db_set, BER_TLSSCS,'b-<');hold on;
semilogy(opt.snr_db_set, BER_GSP,'k-<');hold on;
semilogy(opt.snr_db_set, BER_SICSSP,'r-<');hold on;
semilogy(opt.snr_db_set, BER_OLS,'k-s');hold on;

grid on;
legend('BER-OLSWMM','BER-TLSSCS','BER-GSP','BER-SICSSP','BER-OLS');
xlabel('SNR'), ylabel('BER');
title([' nr= '  ,num2str(opt.nr),' T= '  ,num2str(opt.T),' S= ' ,num2str(opt.S),  ' M-order= ',  num2str(opt.M_order),  ' QAM-en= ' ,num2str(opt.QAM_en),  ' PSK-en= ',num2str(opt.PSK_en)]);
