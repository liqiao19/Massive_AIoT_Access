%"Compressive Sensing Based Massive Access for IoT Relying on Media Modulation Aided Machine Type Communications"
%Accepted by IEEE Transactions on vehicular technology in 2020/06/28
% Written by Li Qiao (qiaoli@bit.edu.cn) from BIT, in 2020/07/28
function [support,NORM] = StrOMP(z,Phi,opt)  %The output is the support
%% step 1 initialization
m = size(Phi,2); % K*M
Block_len = m;
n = size(z,2); % m x n 
a = zeros(m,n); 
NORM=zeros(30,1);
%% Initialization
active_set = [];
res = z;
iii=1;
%% step 2 iteration
while  1
    b = zeros(m, n);
    a = zeros(m, n );
    y = Phi'*res;
    [sort_index] = CoSaMP_MMV_Block(y,Block_len,opt.M); %block
    %     candidate_set = union(active_set, sort_index(1:actset_size));
    active_set_new = sort_index(1:1);
    [ T_active_set_new ] = supp_B_to_O( active_set_new,opt.M);
    for pp = 1:n
        b(T_active_set_new,pp) = Phi(:, T_active_set_new) \ z(:,pp);
    end
    active_set_old=active_set;
    active_set=union(active_set,active_set_new);
    T_active_set=supp_B_to_O(active_set,opt.M);
    for pp = 1:n
        a(T_active_set,pp) = Phi(:, T_active_set) \ z(:,pp);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  begin: select the maximum index
    %-------------------------------------------------------
    lengg=length(active_set);
    a1=zeros(m, n);
    a_abs=abs(a);
    a_abs_T_index=a_abs(T_active_set,:);
    selmax_ind_a=zeros(lengg,n);
    for pp=1:n
        [selmax_ind_a(:,pp)]=selmax_indx(a_abs_T_index(:,pp));
        inindx=selmax_ind_a(:,pp);
        a1(T_active_set(inindx),pp) = Phi(:, T_active_set(inindx)) \ z(:,pp);   
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end
    %residual
    new_res = z-Phi*a1;

    last_norm1=norm(res,'fro');
    last_norm2=norm(new_res,'fro');
    res=new_res;
    NORM(iii)=last_norm1-last_norm2;
    if NORM(iii)<2
        support=active_set_old;
        break;
    end
    iii=iii+1;
end
end