function [Nuser,Nsym,Nbit]=Perfor_AMP(lambda_new_sum,Amid2,H_est,K,Ka,N_seq,Na_seq,T,bit,N_BS,GSMset,GSMindex,Data)
% modified in 2020.08.05
% K=opt.K;
% Ka=opt.Ka;
% N_seq=opt.N_seq;
% Na_seq=opt.Na_seq;
%% Traditional AUD
Es_lamda_AMP2=lambda_new_sum;
% Es_lamda_AMP2=(lambda_new_sum-min(lambda_new_sum))./(max(lambda_new_sum)-min(lambda_new_sum));
[ACTsetThre2]=AUD_thresholdParfor(K,Es_lamda_AMP2,0.5);
Es_ABS_lamda=zeros(K,1);
Es_ABS_lamda(ACTsetThre2)=1;
% Es_ABS_lamda=zeros(opt.N_UE,1);
% Es_ABS_lamda(lambda_new_sum>0.5)=1;
% Es_ABS_lamda(lambda_new_sum<=0.5)=0;
%% Advanced AUD
% [Es_lamda_sorted,sorted_IND]=sort(lambda_new_sum,'descend');
% figure(1);subplot(121);plot(1:opt.N_UE,Es_lamda_sorted);grid on;xlabel('(a)');
% Es_lamda_sorted_diff=zeros(1,opt.N_UE-1);
% for i=1:opt.N_UE-1
%     Es_lamda_sorted_diff(i)=Es_lamda_sorted(i+1)-Es_lamda_sorted(i);
% end
% subplot(122);plot(1:(opt.N_UE-1),Es_lamda_sorted_diff);grid on;xlabel('(b)');
% close
% [ChangeNUM,ChangeIND]=min(Es_lamda_sorted_diff);
% ACT_SETTT=sorted_IND(1:ChangeIND);
% Es_ABS_lamda(ACT_SETTT)=1;
%% Pe calculation
Nuser=nnz(Es_ABS_lamda-Amid2);%error user
%% Correct symbol calculation
XXX=find(Es_ABS_lamda==1);
index_user_single=find(Amid2==1);
[sum_act_user,sub_matrix,sub_matrix_correct]=find_equal(XXX,index_user_single,Ka);
sub_matrix=sub_matrix(1:sum_act_user);
sub_matrix_correct=sub_matrix_correct(1:sum_act_user);
% 计算false alarm user
EU_FA=length(ACTsetThre2)-sum_act_user;
if sum_act_user==0
%     NMSE=1;
    Nsym=1;Ka*T;
    Nbit=1;Ka*T*bit;
else
%     Csym=zeros(opt.T,1);
    EstedSym=zeros(sum_act_user,T);%存储正确检测用户的symbol
%     EstedBit=zeros(opt.T,opt.bit,sum_act_user);%存储正确检测用户的bit
    
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
%         if N_seq~=1
%             MBMsym=sort(MBMsym,2);
%         end
        SymOneslot=MBMsym';
        EstedSym(:,t)=DecodingParfor(GSMset,SymOneslot);
        %-------------------------------------------------------------------%
%         for i=1:sum_act_user
%             if mmmmm(1,i)<1e-1%判断传的是否是全0符号
%                 EstedSym(i,t)=16;
%             end
%         end
        %-------------------------------------------------------------------%
        %     for t=1:opt.T
        %         absH_eq=ABS_H_est(:,(t-1)*opt.N_BS+1:t*opt.N_BS);
        %         MID=sum(absH_eq,2);
        %         Est1=reshape(MID,[opt.N_pilot,opt.N_UE]);
        %         [mmmmm,indmbm1]=sort(Est1,'descend');
        %         MBMsym=indmbm1(1:opt.Ia,XXX(sub_matrix));%MBM符号 estimated  bit信息
        
        %         if opt.Ia~=1
        %             MBMsym=sort(MBMsym);
        %         end
        %         Index0=opt.Index0;
        %         Index0=reshape(Index0(:,t),opt.Ia,opt.N_UE);
        %         MBMreal=Index0(:,index_user_single(sub_matrix_correct));%MBM真实
        %         if opt.Ia~=1
        %             MBMreal=sort(MBMreal);%横向不可以排序；但是只有一行时，表示横向排序
        %         end
        %         Csym(t)=sum_act_user*opt.Ia-nnz(MBMsym-MBMreal);%传输正确的符号数目
    end
    Nsym=((Ka-sum_act_user+EU_FA)*T+nnz( EstedSym-GSMindex(sub_matrix_correct,:)))/(T*(Ka+EU_FA));
%     (Ka-sum_act_user)*T+nnz( EstedSym-GSMindex(sub_matrix_correct,:));
    Nbit=0;
    for i=1:sum_act_user
        EstedBit=de2bi(EstedSym(i,:)-1,bit,'left-msb');
        Nbit=Nbit+nnz(EstedBit-Data(:,:,sub_matrix_correct(i)));
    end
    Nbit=(Nbit+(Ka-sum_act_user+EU_FA)*T*bit)/(bit*T*(Ka+EU_FA));
%     Nbit+(Ka-sum_act_user)*T*bit;
    % 没考虑多个时隙
    %% NMSE calculation
% %     Index=reshape(MBMsym,opt.Ia*sum_act_user,1)+reshape(repmat(((XXX(sub_matrix)-1).*opt.N_pilot)',opt.Ia,1),opt.Ia*sum_act_user,1);
% %     H_est_sel=H_est(sort(Index),:); %信道估计值   estimated  评价算法性能NMSE
% %     Index2=reshape(MBMreal,opt.Ia*sum_act_user,1)+reshape(repmat(((index_user_single(sub_matrix_correct)-1).*opt.N_pilot)',opt.Ia,1),opt.Ia*sum_act_user,1);
% %     Hreal=H2(sort(Index2),:);%真实信道
% %     NMSE=(norm((H_est_sel-Hreal),'fro'))^2/(norm(Hreal,'fro'))^2;
%     NMSE=1;
%     NMSE=(norm((H_est-H2),'fro'))^2/(norm(H2,'fro'))^2;
    % 没考虑多个时隙
end
end