%% "Massive Access in Media Modulation Based Massive Machine Type Communications"
% Written by Li Qiao (qiaoli@bit.edu.cn) in 2020/06/26
% Beijing Institute of Technology

% Cite this: L. Qiao, J. Zhang, Z. Gao, D. W. K. Ng, M. Di Renzo, and M.-S. Alouini, 
%“Massive access in media modulation based massive machine-type communications,”
%IEEE Trans. Wireless Commun., vol. PP, no. PP, Jul. 2021.

clear;close all;
%% Basic parameters
opt.nr=256;         %BS ANTENNA NUMBER
opt.K=500;          %UEs
S0=50;
LLL=length(S0);
opt.mrf=2;           %NUMBER OF RF MIRRORS
opt.M=2^opt.mrf;         
opt.M_order=4;      %M-QAM   e.g., 4QAM 16QAM 64QAM; 4PSK, 8PSK,16PSK,64PSK;
opt.Normuliza=sqrt(3/2/(opt.M_order-1));

opt.damp=0;
opt.max_iter=15;

% Conventional scheme, without media modulation
opt2.mrf=0;
opt2.M=2^opt2.mrf;
opt2.M_order=16;    
opt2.Normuliza=sqrt(3/2/(opt2.M_order-1));

opt.threshold1 = 0.500000;
opt.Na=1;      %Select one radiation pattern
opt.Max_iter=10; %Number of iterations for GSP 
opt.QAM_en=1;        %if 1, using QAM
opt.PSK_en=0;        %if 1, using PSK     
opt.T=12;            % Frame length

opt.awgn_en=1;       % if 1, AWGN
sim_num=1e8;         % Number of simul total
simuamount=10;      % Break condition

opt.snr_db_set=10;    % SNR
opt.snr_num=10.^(opt.snr_db_set/10);

%% First loop
for kk=1:LLL 
    opt.S=S0(kk);
    opt.sparsity=opt.S/opt.K;
    %% Initialization
    opt.index_mr_single=zeros(opt.S,opt.T); 
    index_mr_aggregate=zeros(opt.S,opt.T);
    opt.tx_data=zeros(opt.S,opt.T);
    tx_bin=zeros(opt.S,log2(opt.M_order),opt.T);
    tx_signal=zeros(opt.S,opt.T);
    y=zeros(opt.nr,opt.T);
 
    err_AMP6=0;
    err_MissD6=0;
    err_FalseD6=0;
    ErrBit_AMP6=0;
    ErrSym_AMP6=0;
   
    %% Second loop for thousands of simulations
    for i2222=1:sim_num 
        %Generate H
        H=(randn(opt.nr,opt.K*opt.M)+1i*randn(opt.nr,opt.K*opt.M))/sqrt(2);
        X=zeros(opt.K*opt.M,opt.T);
        %Generate X
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
            
            if opt.QAM_en==1 % 采用QAM
                tx_signal(:,jj) = qammod(opt.tx_data(:,jj),opt.M_order)*opt.Normuliza;
                opt.Set=qammod(0:opt.M_order-1,opt.M_order)*opt.Normuliza;
            elseif opt.PSK_en==1% 采用PSK
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
                P_n(pp)=P_s/opt.snr_num;
                y(:,pp)=H*X(:,pp)+sqrt(P_n(pp)/2).*(randn(opt.nr,1)+1i*randn(opt.nr,1));
            end
        else
            y=H*X;
        end
        opt.sigma2=mean(P_n);
       
        
        %% Proposed DS-AMP
        [LLR_MBM,LLR_Constell,Xe_AMP,Es_lamda_AMP]=DS_AMP(y,H,opt,X);
        %%Max-Min Normulization
        Es_lamda_AMP2=(Es_lamda_AMP-min(Es_lamda_AMP))./(max(Es_lamda_AMP)-min(Es_lamda_AMP));
        % AUD
        [ACTsetThre2]=AUD_threshold(opt,Es_lamda_AMP2,0.5);
       
        % Error Calculation
        [MissD6,FalseD6,Nuser_AMP6,Nsym6,Nbit6]=...
            ErrorCalcu(opt,Xe_AMP,ACTsetThre2,opt.index_user_single,opt.index_mr_single,opt.tx_data,tx_bin);
 
        fprintf('SNR=%d,',opt.snr_db_set);
        fprintf('sim_num=%d,',i2222);
        fprintf('nr=%d,',opt.nr);
        fprintf('T=%d,',opt.T);
        fprintf('S=%d,',opt.S);
        %%  Miss detection
        err_MissD6=err_MissD6+MissD6;
        MissD_AUD=err_MissD6/(i2222*opt.K);
        fprintf('MissD_AUD6=%f,',MissD_AUD);
        %% False detection
        err_FalseD6=err_FalseD6+FalseD6;
        FalseD_AUD=err_FalseD6/(i2222*opt.K);
        fprintf('FalseD_AUD6=%f,',FalseD_AUD);
        %%  Pe calculation
        err_AMP6=err_AMP6+Nuser_AMP6;
        Pe_AUD=err_AMP6/(i2222*opt.K);
        fprintf('Pe_AUD6=%f,',Pe_AUD);
        %% BER calculation
        ErrBit_AMP6=ErrBit_AMP6+Nbit6;
        BER_AUD=ErrBit_AMP6/i2222/opt.S/opt.T/(log2(opt.M)+log2(opt.M_order));
        fprintf('BER_AUD6=%f,',BER_AUD);
        %% SER calculation
        ErrSym_AMP6=ErrSym_AMP6+Nsym6;
        SER_AUD=ErrSym_AMP6/i2222/opt.S/opt.T;
        fprintf('SER_AUD6=%f,\n',SER_AUD);
        
        if i2222>=simuamount
            break;
        end
    end
%     save('result.mat','MissD_AUD','FalseD_AUD','Pe_AUD','BER_AUD','SER_AUD');
end
%% Plot

