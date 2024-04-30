function [support] = SK_StrOMP(z,Phi,opt) 
%% step 1 initialization
% clear;
% load exData2.mat
% Pth=opt.threshold2;
m = size(Phi,2); % K*M
Block_len = m;
n = size(z,2); % m x n 
% NORM=zeros(30,1);
%% Initialization
active_set = [];
res = z;
%% step 2 iteration
for  iii=1:opt.S
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  begin  
    %-------------------------------------------------------
    active_set=union(active_set,active_set_new);
    opt;
    T_active_set=supp_B_to_O(active_set,opt.M);
    %-------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end
    %     %LS
    for pp = 1:n
        a(T_active_set,pp) = Phi(:, T_active_set) \ z(:,pp);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  begin  
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
    res=new_res;
    support=active_set;
end
end