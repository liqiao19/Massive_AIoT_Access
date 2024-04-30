% “Non-coherent transmission” for massive access or “Sequence modulation”
% Main function for time invariant system (TIS) considering OFDM
% Single sequence is transmitted in the time domain (TD)
% Active device error rate (ADER),BER,NMSE vs Number of observations Lseq
% Written by Li Qiao from Beijing Institute of Technology,
% Email: qiaoli@bit.edu.cn
% 2020/12/25
% 观测30，信噪比5dB
clear;close all;
tic
%% Parameters
% sequence modulation parameters
% Ob_type='Bonulli';'Gaussian';% Distribution of the elements of sequences
L_seq=30;   % time overhead of sequences
N_seq=2;      % total number of sequence
Na_seq=1;     % active number of sequence
bit=floor(log2(nchoosek(N_seq,Na_seq)));  % Embedded bits
[ICT]=FunSelct2(N_seq);
% ICT中存储了I=opt.I的所有可能的激活选择的排列组合,除了全不选和全选两种
GSMset=ICT(1:2^bit,1:Na_seq,Na_seq);  % activation patterns
T=1;      % consequent time slots (one sequence is transimitted within a time slot)
% parameters of system && channel model
ang_spread_d=10;   % degree
ang_spread_r = ang_spread_d/180*pi;  % degree to rad
K = 100;           % number of potential UDs in the cell
sparsity=0.1;
Ka=K*sparsity;   % number of active UDs
N_BSset = 16;8:10:70;16;32;         % number of BS antennas
awgn_en=1;             %0无噪声  1加性噪声
sim_num=1e8;           %每个信噪比下的仿真次数
SNR = 20;0:5:20;          % SNR in dB
snr_num=10.^(SNR/10);
LLL=length(N_BSset);   %
N =  512;1024; 64;        % number of subcarriers
Bw = 10e6;         % system bandwidth
fs = Bw;       % frequency of sampling
N_cp = 32;         % length of cyclic prefix
sigma2_alpha = 1;  % average power of path, generally, sigma2_alpha = 1;
Lp_min = 8;        % number of paths
Lp_max = 14;
f_c=6e9;           % carrier frequence Hz
v_max=0;50e3/3600;   % maximum device's velocity m/s
f_d=v_max*f_c/3e8; % maximum Doppler frequence shift Hz
% parameters of the algorithm
damp = 0.3;       % control the convergence speed (prevent algorithm from diverging)
niter = 200;      % maximum number of AMP iterations
tol = 1e-5;       % termination threshold
Nmethod=3;        % number of algorithms
Uomp=zeros(1,LLL);Somp=zeros(1,LLL);Bomp=zeros(1,LLL);NMSEomp=zeros(1,LLL);
U_ampQ=zeros(1,LLL);S_ampQ=zeros(1,LLL);B_ampQ=zeros(1,LLL);NMSE_ampQ=zeros(1,LLL);
U_ampK=zeros(1,LLL);S_ampK=zeros(1,LLL);B_ampK=zeros(1,LLL);NMSE_ampK=zeros(1,LLL);
U_ampK2=zeros(1,LLL);S_ampK2=zeros(1,LLL);B_ampK2=zeros(1,LLL);NMSE_ampK2=zeros(1,LLL);
U_ampT=zeros(1,LLL);S_ampT=zeros(1,LLL);B_ampT=zeros(1,LLL);NMSE_ampT=zeros(1,LLL);
NMSE_lmmse=zeros(1,LLL); NMSE_ls=zeros(1,LLL);
% startmatlabpool(LLL)
for i2222=1:10%80  %每个信噪比下仿真的次数
    %% generate device activity
    ActUE=sort(randperm(K,Ka));
    A_UE_index=zeros(K,1);          % K*1
    A_UE_index(ActUE)=1;            % device activity
    %% sequence modulation   bit —> index (symbol)
    %空间调制和一般化的空间调制，其调制仿真最好模拟bit流的实际传输
    %无编码下的bit传输和调制
    A_index=zeros(N_seq*K,T);         % K N_seq * T
    %transmit bits
    Data = randi([0 1],T,bit,Ka);
    % transmit symbols
    GSMindex=zeros(Ka,T);
    for Iue=1:Ka
        GSMindex(Iue,:)=bi2de(Data(:,:,Iue),'left-msb')+1;
        for t=1:T
            A_index((ActUE(Iue)-1)*N_seq+GSMset(GSMindex(Iue,t),:),t)=1;
        end
        %二进制转换到十进制，加1避免矩阵index出现0
    end
    A_seq_ind=find(A_index==1);
    %% generate sequences matrix; distribution restricted
    S0=(2*(randi([0 1],L_seq, N_seq*K)-0.5)+1i*2*(randi([0 1],L_seq, N_seq*K)-0.5))/sqrt(2*L_seq);%symetric Bonuelli
    S1=(randn(L_seq, N_seq*K)+1i*randn(L_seq, N_seq*K))/sqrt(2*L_seq);% Gaussian
    % normalize the pilot matrix
    for j = 1:K
        S0(:,j) = S0(:,j)./norm(S0(:,j));
        S1(:,j) = S1(:,j)./norm(S1(:,j));
    end
    S0_pull=S0(:,A_seq_ind);
    S0_pull=S0_pull.';
    %% generate H, H can be time-varying; UE sparsity excluded
    F=dftmtx(N)/sqrt(N);
    for kk=1:LLL  %遍历天线数
        %         L_seq=L_seqset(kk);
        N_BS=N_BSset(kk);
        Y00=zeros(N,L_seq,N_BS);
        Y=zeros(N,L_seq,N_BS);
        H=zeros(N,N_BS,Ka);
        H_a=zeros(N,N_BS,Ka);
%         for k = 1:Ka
%             Lp=randi([Lp_min,Lp_max],1,1);
%             ang_spread_r = ang_spread_d/180*pi;  % degree to rad
%             cluster_center = unifrnd(-pi/2, pi/2);
%             AOAs = unifrnd(cluster_center-ang_spread_r/2, cluster_center+ang_spread_r/2, 1, Lp);
%             AOAs(AOAs > pi/2) = AOAs(AOAs > pi/2) - pi;
%             AOAs(AOAs < -pi/2) = AOAs(AOAs < -pi/2) + pi;
%             n_BS = (0:N_BS-1).';
%             A_BS = exp(-1i*pi*n_BS*sin(AOAs));           % steering vector
%             PowerdB=[0,-1*sort(randperm(30,Lp-1))];  % multi-path power in dB for MTD k
%             Power0=10.^(PowerdB/10);
%             Power0=Power0/sum(Power0);                   % power normalization
%             [H_f2,H_angl] = channel_model_zhangshun_v3(f_d,N_BS,N,fs,Lp,F,A_BS,Power0,L_seq);
%             H(:,:,k)=H_f2;
%             H_a(:,:,k)=H_angl;
%             for i=1:N_BS
%                 Y00(:,:,i)=Y00(:,:,i)+ H_f2(:,i)*S0_pull(k,:);
%             end
%         end
        %%-----------------Begin----frequency selective----------------%%
        Lp=randi([Lp_min,Lp_max],1,Ka);
        n_BS = (0:N_BS-1).';
        n_freq=0:N-1;
        A_R=dftmtx(N_BS)/sqrt(N_BS);
        for k =1:Ka
            Dtemp=(N_cp/Bw)*rand(Lp(k),1);                 % Delay
            Delay=exp(-2*1i*pi*Dtemp*(-Bw/2+Bw*n_freq/N)); % Delay
            cluster_center = unifrnd(-pi/2, pi/2);
            AOAs = unifrnd(cluster_center-ang_spread_r/2, cluster_center+ang_spread_r/2, 1, Lp(k));
            AOAs(AOAs > pi/2) = AOAs(AOAs > pi/2) - pi;
            AOAs(AOAs < -pi/2) = AOAs(AOAs < -pi/2) + pi;
            A_BS = sqrt(1/N_BS)*exp(-1i*pi*n_BS*sin(AOAs));% steering vector
            h0=sqrt(1/2)*(randn(Lp(k),1)+1i*randn(Lp(k),1));           % channel gain
            Delay2=h0.*Delay;
            %%--------------------If Frequency-selective--------------------%%
            hk=sqrt(N_BS/Lp(k))*A_BS*Delay2;
            hk_T=hk.';
%             figure;plot(abs(hk_T))
            H(:,:,k)=hk_T;
            hk_ang = (hk.')*(A_R');
            H_a(:,:,k)=hk_ang;
            %%--------------------End Frequency-selective--------------------%%
            for i=1:N_BS
                Y00(:,:,i)=Y00(:,:,i)+ hk_T(:,i)*S0_pull(k,:);
            end
        end
        H_eff=zeros(K*N_seq,N_BS,N);
        for i=1:N
            Ht=reshape(H(i,:,:),N_BS,Ka);
            H_eff(A_seq_ind,:,i)=Ht.';
        end
        H_eff2=zeros(K*N_seq,N_BS,N);
        for i=1:N
            Ht=reshape(H_a(i,:,:),N_BS,Ka);
            H_eff2(A_seq_ind,:,i)=Ht.';
        end
        if awgn_en==0
            Y=Y00;
            Pn=1e-6;
        else
            PnV=zeros(1,N_BS);
            for i=1:N_BS
                Ps=(norm(Y00(:,:,i),'fro'))^2/ size(Y00(:,:,i),2)/ size(Y00(:,:,i),1);
                PnV(i)=Ps/snr_num;
                Y(:,:,i)=Y00(:,:,i)+sqrt(PnV(i)/2)*(randn(N,L_seq)+1i*randn(N,L_seq));
            end
            Pn=mean(PnV);
        end
        Nuser_omp=0;Nsym_omp=0;Nbit_omp=0;nmse_omp=0; %OMP
        Nuser_ampQ=0;Nsym_ampQ=0;Nbit_ampQ=0;nmse_ampQ=0; %AMP QIAO LI
        Nuser_ampK=0;Nsym_ampK=0;Nbit_ampK=0;nmse_ampK=0; %AMP KE MALONG
        Nuser_ampK2=0;Nsym_ampK2=0;Nbit_ampK2=0;nmse_ampK2=0; %AMP KE MALONG
        Nuser_ampT=0;Nsym_ampT=0;Nbit_ampT=0;nmse_ampT=0; %AMP TANG ZIHAN
        nmse_ls=0;nmse_lmmse=0;
        SS=S0_pull.';
        for j=1:N
            YY=reshape(Y(j,:,:),L_seq,N_BS);
            HH=reshape(H_eff(:,:,j),K*N_seq,N_BS);
            A_R=dftmtx(N_BS)/sqrt(N_BS);
            YY_a=YY*(A_R');
            HH_a=reshape(H_eff2(:,:,j),K*N_seq,N_BS);
            H_benchmark=reshape(H(j,:,:),N_BS,Ka);
            H_benchmark=H_benchmark.';
            varH = norm(H_benchmark(:))^2/length(H_benchmark(:));
            %% Algorithms
            %我的代码，角度域
            [He_ampQ, lambda_ampQ,pai_tmpQ,iterQ,~] = AMP_NNSPL_2D_v2(YY_a, S0, Pn,damp,niter,tol,K,N_seq);
            [EUampQ,ESampQ,EBampQ]=Perfor_AMP_angle_v3(Pn,lambda_ampQ,A_UE_index,He_ampQ,K,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data);
            Enmse_ampQ=(norm((He_ampQ-HH_a),'fro'))^2/(norm(HH_a,'fro'))^2;
            
            %柯博代码，空频域
            [He_ampK, lambda_ampK,pai_tmp,iter2,~] = GMMV_AMP(YY, S0, Pn,damp,niter,tol,K,N_seq,varH);
            [EUampK,ESampK,EBampK]=Perfor_AMP_v3(lambda_ampK,A_UE_index,He_ampK,K,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data);
            Enmse_ampK=(norm((He_ampK-HH),'fro'))^2/(norm(HH,'fro'))^2;
            
            %柯博代码，角度域
            [He_ampK2, lambda_ampK2,~, iter_ampK2] = GMMV_AMP_angle_v2(YY_a, S0, Pn,damp,niter,tol,K,N_seq,varH);
            [EUampK2,ESampK2,EBampK2]=Perfor_AMP_angle_v3(Pn,lambda_ampK2,A_UE_index,He_ampK2,K,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data);
            Enmse_ampK2=(norm((He_ampK2-HH_a),'fro'))^2/(norm(HH_a,'fro'))^2;
            
            %SOMP，空频域
            [He_omp, lambda_omp,~, iter0] = SOMP(YY, S0, Pn,K,N_seq);
            [EUomp,ESomp,EBomp]=Perfor_OMP_v3(lambda_omp,A_UE_index,He_omp,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data);
            Enmse_omp=(norm((He_omp-HH),'fro'))^2/(norm(HH,'fro'))^2;
            
            %%唐子涵代码，空频域
            [He_ampT,~,lambda_ampT] = AMP_tang(YY, S0,N_seq,sparsity,varH);
            [EUampT,ESampT,EBampT]=Perfor_AMP_tang_v3(lambda_ampT,A_UE_index,He_ampT,K,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data);
            Enmse_ampT=(norm((He_ampT-HH),'fro'))^2/(norm(HH,'fro'))^2;
            
            %             [He_lmmse] = OracleLMMSE(YY,SS,Pn);
            %             Enmse_lmmse=(norm((He_lmmse-H_benchmark),'fro'))^2/(norm(H_benchmark,'fro'))^2;
            %             [He_ls] = OracleLS(YY,SS);
            %             Enmse_ls=(norm((He_ls-H_benchmark),'fro'))^2/(norm(H_benchmark,'fro'))^2;
            %% errors accumulation
            Nuser_omp=Nuser_omp+EUomp;    Nsym_omp=Nsym_omp+ESomp;    Nbit_omp=Nbit_omp+EBomp;    nmse_omp=nmse_omp+Enmse_omp; %OMP
            Nuser_ampQ=Nuser_ampQ+EUampQ; Nsym_ampQ=Nsym_ampQ+ESampQ; Nbit_ampQ=Nbit_ampQ+EBampQ; nmse_ampQ=nmse_ampQ+Enmse_ampQ; %AMP QIAO LI
            Nuser_ampK=Nuser_ampK+EUampK; Nsym_ampK=Nsym_ampK+ESampK; Nbit_ampK=Nbit_ampK+EBampK; nmse_ampK=nmse_ampK+Enmse_ampK; %AMP KE MALONG
            Nuser_ampK2=Nuser_ampK2+EUampK2; Nsym_ampK2=Nsym_ampK2+ESampK2; Nbit_ampK2=Nbit_ampK2+EBampK2; nmse_ampK2=nmse_ampK2+Enmse_ampK2; %AMP KE MALONG
            Nuser_ampT=Nuser_ampT+EUampT; Nsym_ampT=Nsym_ampT+ESampT; Nbit_ampT=Nbit_ampT+EBampT; nmse_ampT=nmse_ampT+Enmse_ampT; %AMP TANG ZIHAN
            %             nmse_lmmse=nmse_lmmse+Enmse_lmmse; nmse_ls=nmse_ls+Enmse_ls;
        end
        Uomp(kk)=Uomp(kk)+Nuser_omp;      Somp(kk)=Somp(kk)+Nsym_omp;      Bomp(kk)=Bomp(kk)+Nbit_omp;      NMSEomp(kk)=NMSEomp(kk)+nmse_omp;
        U_ampQ(kk)=U_ampQ(kk)+Nuser_ampQ; S_ampQ(kk)=S_ampQ(kk)+Nsym_ampQ; B_ampQ(kk)=B_ampQ(kk)+Nbit_ampQ; NMSE_ampQ(kk)=NMSE_ampQ(kk)+nmse_ampQ;
        U_ampK(kk)=U_ampK(kk)+Nuser_ampK; S_ampK(kk)=S_ampK(kk)+Nsym_ampK; B_ampK(kk)=B_ampK(kk)+Nbit_ampK; NMSE_ampK(kk)=NMSE_ampK(kk)+nmse_ampK;
        U_ampK2(kk)=U_ampK2(kk)+Nuser_ampK2; S_ampK2(kk)=S_ampK2(kk)+Nsym_ampK2; B_ampK2(kk)=B_ampK2(kk)+Nbit_ampK2; NMSE_ampK2(kk)=NMSE_ampK2(kk)+nmse_ampK2;
        U_ampT(kk)=U_ampT(kk)+Nuser_ampT; S_ampT(kk)=S_ampT(kk)+Nsym_ampT; B_ampT(kk)=B_ampT(kk)+Nbit_ampT; NMSE_ampT(kk)=NMSE_ampT(kk)+nmse_ampT;
        %         NMSE_lmmse(kk)=NMSE_lmmse(kk)+nmse_lmmse;NMSE_ls(kk)=NMSE_ls(kk)+nmse_ls;
        %实时输出
        Tuser=i2222*K*N;
        Tbit=i2222*N;
        Tsym=i2222*N;
        Tnmse=i2222*N;
        fprintf('SNR=%d,',SNR);
        fprintf('sim_num=%d,',i2222);
        fprintf('N_UE=%d,',K);
        fprintf('A_UE=%d,',Ka);
        fprintf('N_BS=%d,',N_BS);
        fprintf('L_seq=%d,',L_seq);
        fprintf('N_pilot=%d,---',N_seq);
        
        fprintf('SOMP空频域:');
        fprintf('nmse=%f,',NMSEomp(kk)/Tnmse);
        fprintf('Pe=%f,',Uomp(kk)/Tuser);
        fprintf('SER=%f,',Somp(kk)/Tsym);
        fprintf('BER=%f---',Bomp(kk)/Tbit);
        fprintf('AMP_QL角度域:');
        fprintf('nmse=%f,',NMSE_ampQ(kk)/Tnmse);
        fprintf('Pe=%f,',U_ampQ(kk)/Tuser);
        fprintf('SER=%f,',S_ampQ(kk)/Tsym);
        fprintf('BER=%f---',B_ampQ(kk)/Tbit);
        fprintf('AMP_KML空频域:');
        fprintf('nmse=%f,',NMSE_ampK(kk)/Tnmse);
        fprintf('Pe=%f,',U_ampK(kk)/Tuser);
        fprintf('SER=%f,',S_ampK(kk)/Tsym);
        fprintf('BER=%f---',B_ampK(kk)/Tbit);
        fprintf('AMP_KML角度域:');
        fprintf('nmse=%f,',NMSE_ampK2(kk)/Tnmse);
        fprintf('Pe=%f,',U_ampK2(kk)/Tuser);
        fprintf('SER=%f,',S_ampK2(kk)/Tsym);
        fprintf('BER=%f---',B_ampK2(kk)/Tbit);
        fprintf('AMP_tang:');
        fprintf('nmse=%f,',NMSE_ampT(kk)/Tnmse);
        fprintf('Pe=%f,',U_ampT(kk)/Tuser);
        fprintf('SER=%f,',S_ampT(kk)/Tsym);
        fprintf('BER=%f---\n',B_ampT(kk)/Tbit);
        %         fprintf('LMMSE:');
        %         fprintf('nmse=%f,',NMSE_lmmse(kk)/Tnmse);
        %         fprintf('LS:');
        %         fprintf('nmse=%f\n',NMSE_ls(kk)/Tnmse);
    end
    fprintf('-------------------------------------\n');
end
closematlabpool
%% 计算误码率
save SeqModulation-Lseq-angle-BS-Ch.mat i2222 SNR L_seq N_seq K N Ka T N_BSset bit Uomp Somp Bomp NMSEomp U_ampQ S_ampQ B_ampQ NMSE_ampQ U_ampK S_ampK B_ampK NMSE_ampK U_ampK2 S_ampK2 B_ampK2 NMSE_ampK2 U_ampT S_ampT B_ampT NMSE_ampT

TUSER=i2222*K*N;
TBER=i2222*N;
% Tsym=i2222*Ka*T*N_block;
TNMSE=i2222*N;

X_set=N_BSset;

figure;
semilogy(X_set, NMSE_ampQ/TNMSE,'r--*');hold on;
semilogy(X_set, NMSE_ampK/TNMSE,'b--+');hold on;
semilogy(X_set, NMSE_ampK2/TNMSE,'k--+');hold on;
semilogy(X_set, NMSE_ampT/TNMSE,'k--o');hold on;
semilogy(X_set, NMSEomp/TNMSE,'g-<');hold on;
% semilogy(X_set, NMSE_lmmse/TNMSE,'k-+');hold on;
% semilogy(X_set, NMSE_ls/TNMSE,'k-<');hold on;
grid on;
legend('AMP-QL','AMP-KML空频域','AMP-KML角度域','AMP-唐子涵','SOMP')
axis([8,70,1e-3,1e2]);
xlabel('Number of antennas')
ylabel('NMSE')
title(['SNR=',num2str(SNR),'dB','  ', 'L-seq=', num2str(L_seq), ' N-seq=',num2str(N_seq)]);

figure;
semilogy(X_set, U_ampQ/TUSER,'r--*');hold on;
semilogy(X_set, U_ampK/TUSER,'b--+');hold on;
semilogy(X_set, U_ampK2/TUSER,'k--+');hold on;
semilogy(X_set, U_ampT/TUSER,'k--o');hold on;
semilogy(X_set, Uomp/TUSER,'g-<');hold on;
grid on;
legend('AMP-QL','AMP-KML空频域','AMP-KML角度域','AMP-唐子涵','SOMP')
xlabel('Number of antennas')
ylabel('$P_e$','interpreter','latex')
title(['SNR=',num2str(SNR),'dB','  ', 'L-seq=', num2str(L_seq), ' N-seq=',num2str(N_seq)]);

figure;
semilogy(X_set, B_ampQ/TBER,'r--*');hold on;
semilogy(X_set, B_ampK/TBER,'b--+');hold on;
semilogy(X_set, B_ampK2/TBER,'k--+');hold on;
semilogy(X_set, B_ampT/TBER,'k--o');hold on;
semilogy(X_set, Bomp/TBER,'g-<');hold on;
grid on;
legend('AMP-QL','AMP-KML空频域','AMP-KML角度域','AMP-唐子涵','SOMP')
grid on;
xlabel('Number of antennas')
ylabel('BER')
title(['SNR=',num2str(SNR),'dB','  ', 'L-seq=', num2str(L_seq), ' N-seq=',num2str(N_seq)]);
toc
