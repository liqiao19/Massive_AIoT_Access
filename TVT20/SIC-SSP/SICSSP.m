%"Compressive Sensing Based Massive Access for IoT Relying on Media Modulation Aided Machine Type Communications"
%Accepted by IEEE Transactions on vehicular technology in 2020/06/28
% Written by Li Qiao (qiaoli@bit.edu.cn) from BIT, in 2020/07/28
function [mbm_dett,bbb,ott]=SICSSP(support_est,y,H,opt)
% global M   T Na  K  Max_iter  M_order
support_num=length(support_est);
support_est_M=zeros(1,support_num*opt.M);
for qqqq=1:support_num
    support_est_M(((qqqq-1)*opt.M+1):qqqq*opt.M)=((support_est(qqqq)-1)*opt.M+1):(support_est(qqqq)*opt.M);
end
%%  method-1 Using  SP_modify_std_DCS_plus_ql_v1 algorithm
% debug_en=0;
temppp_ = (0:support_num-1)*opt.M;
% aaa=zeros(support_num,opt.T);
bbb=zeros(support_num,opt.T);
mbm_dett=zeros(support_num,opt.T);
ott=0;  %ott, OOO  
% COMset=[];
temp_plus = [0:support_num-1]*opt.M;
for kkkk=1:opt.T      
    btmp=zeros(support_num*opt.M,1);
    
    [~,btmp,~ ,used_ietr,~,Phinew,znew,temp_plus_new] = SP_modify_std_DCS_plus_ql_v1ast(btmp,temp_plus,y(:,kkkk),H(:,support_est_M,1),opt);
    for iiijjjkkk=1:support_num-1
        [~,btmp,~ ,used_ietr,~,Phinew,znew,temp_plus_new] = SP_modify_std_DCS_plus_ql_v1ast(btmp,temp_plus_new,znew,Phinew,opt);
    end
   
    ott=ott+used_ietr;
    mbm_dett(:,kkkk)=find(btmp)'-temppp_-1;
    bbb(:,kkkk)=btmp(btmp~=0);
   
end   
end