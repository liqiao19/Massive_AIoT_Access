function [Nuser,Nsym,Nbit]=Perfor_OMP(index_A_UE,A_UE_index,H_est,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data)
% modified in 2020.08.05
% K=opt.K;
% Ka=opt.Ka;
% N_seq=opt.N_seq;
% Na_seq=opt.Na_seq;
%% Pe calculation
Nuser=nnz(index_A_UE'-A_UE_index);%error user
%% Correct symbol calculation
XXX=find(index_A_UE~=0);
index_user_single=find(A_UE_index==1);
[sum_act_user,sub_matrix,sub_matrix_correct]=find_equal(XXX,index_user_single,Ka);
sub_matrix=sub_matrix(1:sum_act_user);
sub_matrix_correct=sub_matrix_correct(1:sum_act_user);
% 计算false alarm user
EU_FA=length(XXX)-sum_act_user;
if sum_act_user==0
%     NMSE=1;
    Nsym=1;Ka*T;
    Nbit=1;Ka*T*bit;
else
    EstedSym=zeros(sum_act_user,T);%存储正确检测用户的symbol
    ABS_H_est=abs(H_est);
    INDEX=zeros(sum_act_user*N_seq,1);%活跃用户扩展的编号
    for i=1:sum_act_user
        INDEX((i-1)*N_seq+1:i*N_seq)=(XXX(sub_matrix(i))-1)*N_seq+1:XXX(sub_matrix(i))*N_seq;
    end
    for t=1:T
        absH_eq=ABS_H_est(INDEX,(t-1)*N_BS+1:t*N_BS);
        MID=sum(absH_eq,2)/N_BS;
        Est1=reshape(MID,[N_seq,sum_act_user]);
        [mmmmm,indmbm1]=sort(Est1,'descend');
        MBMsym=indmbm1(1:Na_seq,:);%MBM符号 estimated  bit信息
        SymOneslot=MBMsym';
        EstedSym(:,t)=DecodingParfor(GSMset,SymOneslot);
    end
    Nsym=((Ka-sum_act_user+EU_FA)*T+nnz( EstedSym-GSMindex(sub_matrix_correct,:)))/(T*(Ka+EU_FA));
    Nbit=0;
    for i=1:sum_act_user
        EstedBit=de2bi(EstedSym(i,:)-1,bit,'left-msb');
        Nbit=Nbit+nnz(EstedBit-Data(:,:,sub_matrix_correct(i)));
    end
    Nbit=(Nbit+(Ka-sum_act_user+EU_FA)*T*bit)/(bit*T*(Ka+EU_FA));
    %% NMSE calculation
%     NMSE=1;
end
end