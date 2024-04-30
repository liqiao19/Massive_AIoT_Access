%"Compressive Sensing Based Massive Access for IoT Relying on Media Modulation Aided Machine Type Communications"
%Accepted by IEEE Transactions on vehicular technology in 2020/06/28
% Written by Li Qiao (qiaoli@bit.edu.cn) from BIT, in 2020/07/28
function [a,btmp,support,used_iter,SSS1,Phinew,znew,temp_plus_new] = SP_modify_std_DCS_plus_ql_v1ast(btmp,temp_plus,z,Phi,opt)
d = size(Phi,2); % d-dimensinal signal
active_num=d/opt.M;
% Nf = size(Phi,3);
%Begin CoSaMP
err = 10^(-6);  % residual error
a = zeros(d,1); % Initialize the recovery signal
v = z;
it=0;
stop = 0;
T_new = [];
while ~stop
    T_old = T_new;
    y = Phi'*v;
    y_temp = abs(y).^2;
    ix_temp = reshape(y_temp,opt.M,[]);
    [tmp, ix] = sort(abs(ix_temp), 'descend');
    ix_temp2 = ix(1:opt.Na,:);
    temp_ = [0:active_num-1]*opt.M;
    Omega = temp_ + ix_temp2;
    
    T = union(Omega, T_new);
    %Signal Estimation
    b = zeros(d,1);
    b(T) = Phi(:, T) \ z;
    %Prune
    b_temp = abs(b.^2);
    ix_temp2 = reshape(b_temp,opt.M,[]);
    [tmp2, ix2] = sort(abs(ix_temp2), 'descend');
    QQTTVV=sort(tmp2(1,:),'descend');
    [AAS,SSS1,SSS2]=find_equal(QQTTVV,tmp2(1,:),active_num);%SSS1
    ix_temp3 = ix2(1:opt.Na,:);
    T_new = temp_ + ix_temp3;
    %%
    a = zeros(d, 1);
    %     btmp=zeros(d,1);
    
    a(T_new) = Phi(:, T_new) \ z;
    
    
    
    v = z - Phi(:,:)*a;
    
    %        z=znew;
    
    
    
    it = it + 1;
    
    if (it >= opt.Max_iter|| isequal(T_old,T_new) ||  norm(v)<=err*norm(z))
        stop = 1;
        CONFIrealset=temp_plus+mod(T_new-1,opt.M)+1;
        CONFIset=[CONFIrealset(SSS1==1)];
        ctmp=zeros(d,1); 
        btmp(CONFIset)=a(T_new(SSS1==1));
        if opt.QAM_en==1 % QAM
            tempp = qamdemod(btmp(CONFIset)/opt.Normuliza,opt.M_order);
            ref_sig = qammod(tempp,opt.M_order)*opt.Normuliza;
        elseif opt.PSK_en==1% PSK
            tempp = pskdemod(btmp(CONFIset),opt.M_order);
            ref_sig = pskmod(tempp,opt.M_order);
        end
        ctmp(T_new(SSS1==1))=ref_sig;%ref_sig corresponds to the true constellations
        znew=z-Phi(:,:)*ctmp;  %  update the residuals
        usersetnew=setdiff(T_new,T_new(SSS1==1));
        %Sample Update
        LLL=length(usersetnew);
        usersetneww=ceil(usersetnew/opt.M);
        %Iteration counter
        support_est_M=zeros(1,LLL*opt.M);
        for qqqq=1:LLL
            support_est_M(((qqqq-1)*opt.M+1):qqqq*opt.M)=((usersetneww(qqqq)-1)*opt.M+1):(usersetneww(qqqq)*opt.M);
        end
        Phinew=Phi(:,support_est_M);
    end
end %End CoSaMP iteration
%% SIC steps
temp_plus_new=setdiff(temp_plus,[temp_plus(SSS1==1)]);
used_iter = it;
support=[];find(btmp)';