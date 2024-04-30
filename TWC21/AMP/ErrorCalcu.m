function [MissDetec,FalseDetec,Euser,Esym,Ebit]=ErrorCalcu(opt,Xe,ACTset,index_user_single,index_mr_single,tx_data,tx_bin)
K=opt.K; % number of users
Nt=opt.M; % number of antennas per user
T=opt.T;
% sigma2=opt.sigma2; % noise variance of channel
alphabet_set=opt.Set; % the modulation symbol set
% M = length(alphabet_set);
Normuliza=sqrt(3/2/(opt.M_order-1));
%Decode
Est=zeros(Nt,K,T);
for j=1:T
    Est(:,:,j)=reshape(Xe(j,:),[Nt,K]);
end
[symMBM0,indmbm0]=max(Est);
symMBM=reshape(symMBM0(1,:,:),[K,opt.T]);
indmbm=reshape(indmbm0(1,:,:),[K,opt.T]);
% MeanMBMsym=mean(abs(symMBM),2);
%Nuser
Es_ABS_lamda=zeros(1,K);
Es_ABS_lamda(ACTset)=1;
Act=zeros(1,opt.K);
Act(index_user_single)=1;  % true MM index
MissDetec=length(find(Es_ABS_lamda(index_user_single)==0));
Temm=setdiff(1:opt.K,index_user_single);
FalseDetec=length(find(Es_ABS_lamda(Temm)==1));
Euser=nnz(Es_ABS_lamda-Act);%error user

% Extract correctly detected UEs
XXX=find(Es_ABS_lamda==1);
[sum_act_user,sub_matrix,sub_matrix_correct]=find_equal(XXX,index_user_single,opt.S);
sub_matrix=sub_matrix(1:sum_act_user);
sub_matrix_correct=sub_matrix_correct(1:sum_act_user);
%Buser
Suser=(opt.S-sum_act_user)*opt.T;
Buser=(opt.S-sum_act_user)*opt.T*(log2(Nt)+log2(opt.M_order));
%Nmbm   based on correct user estimation
if sum_act_user==0
    Esym=Suser;
    Ebit=Buser;
    MissDetec=opt.S;
    FalseDetec=0;
else 
    MBMsym=indmbm(XXX(sub_matrix),:)-1;    %MBM symbols estimated
    CONSTlsym=symMBM(XXX(sub_matrix),:);  %constellation symbols estimated
    MM1=zeros(length(sub_matrix_correct),T);
    MM2=zeros(length(sub_matrix_correct),T);
    MM1=abs(MBMsym-index_mr_single(sub_matrix_correct,:));
    Nmbm=0;
    %Bmbm
    for j=1:opt.T
        Tx_ant_bin(:,:,j)=de2bi(index_mr_single(:,j),log2(opt.M),'left-msb');%bits of index selection
        Rx_ant_bin(:,:,j) = de2bi(MBMsym(:,j),log2(opt.M),'left-msb');% detected
    end
    Bmbm=nnz(Rx_ant_bin-Tx_ant_bin(sub_matrix_correct,:,:));
    qpskDemod = comm.QPSKDemodulator('BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio', ...
        'VarianceSource','Input port');
    [m,n]=size(CONSTlsym);
    % demodSig = qpskDemod(reshape(CONSTlsym,m*n,1),opt.sigma2);   % Demodulate
    % LLR=LLR_soft(CONSTlsym,opt);% LLR
    if opt.QAM_en==1 % adopting QAM
        rx_dec = qamdemod(CONSTlsym/Normuliza,opt.M_order);
    elseif PSK_en==1
        rx_dec = pskdemod(CONSTlsym,opt.M_order);   %decoded symbols  decimalism 
    end
    %Nconstel
    MM2=abs(rx_dec-tx_data(sub_matrix_correct,:));
    MM3=MM1+MM2;
    Nconstel=nnz(MM3);
    %Bconstel
    rx_bin =  de2bi(rx_dec,log2(opt.M_order),'left-msb');
    tx_bin_com = [];
    for pp=1:opt.T
        tx_bin_com = [tx_bin_com;tx_bin(sub_matrix_correct,:,pp)];
    end
    Bconstel=nnz(rx_bin-tx_bin_com);
    Ebit=Bconstel+Bmbm+Buser;
    Esym=Nconstel+Suser;
end
end