function   [X_est,MF_out,Index]=AMP_tang(Y,pilot_Mat,N_seq,sparsity,xvar)

% This function performs Activity Detection using AMP algorithms
N_pilot=N_seq;

%% Validity Check
% check the dimensions
[rec_r,M]=size(Y); % M denotes the numebr of antennas


[P,K]=size(pilot_Mat); % D: the pilot dimension, K:  number of whole users times opt.N_pilot
if(P>K/2)
    error('pilot dimension should be sholud be smaller than potential users');
end

if(P~=rec_r)
    error('dimension of received signal and pilot_Mat are incompatible');
end

% % number of active users
% n_active=opt.A_UE;

% check that pilot_Mat has almost Frobenious norm
pilot_norm=norm(pilot_Mat,'fro')^2;
threshold=0.1;

if (pilot_norm > K*(1+threshold) || pilot_norm < K*(1-threshold))
    error('It seems that the pilot sequence is not well normalized.');
end

% the fraction of active users
lambda=sparsity;

% channel gain or varaince of signal
G=sqrt(xvar);

%% Run AMP
X=zeros(K,M);
R=zeros(P,M);

niter=5;
R=Y;
% taoSE=Tao(end);
tau_vec_est=zeros(niter,1);

% ind_all=1:K;
% ind_avg=setdiff(ind_all, active_ind);

% set the residual to Y
for i=1:niter
    % do match filtering (MF): note that the rows of X encode the activity
    % of the users
    
    % output of matched filter
    MF_out=pilot_Mat'*R+X;
    
    % estimate tau_vec from the output of the Matched Filer (genie-aided)
%           tau_vec_est(i)=mean(mean(abs(MF_out(ind_avg:end,:)).^2,2)); %saeid adopt
    
    tau_vec_est(i) =trace(R*R')/(P*M);%liu liang adopt
    
    % compute the l_2 norm of rows of matched filter output
    x_gain=sum(abs(MF_out).^2,2);
    Delta=1/tau_vec_est(i) - 1/(tau_vec_est(i)+G^2);
%     Delta=1/taoSE - 1/(taoSE+G^2);
    Pi=Delta*x_gain/M;
    Phi=log(1+G^2/tau_vec_est(i));
%     Phi=log(1+G^2/taoSE);
    Index=Pi-Phi;
    Index(Index>30)=30;
    Index(Index<-30)=-30;
    Omega=exp(M*(Index));
    Temp=reshape(Omega,[N_pilot,K/N_pilot]);
    Temp2=sum(Temp);
    Temp3=repmat(Temp2,N_pilot,1);
    Temp4=reshape(Temp3,[K,1]);
    OmegaLine=Omega./(N_pilot*(1-lambda)/lambda+Temp4);
    Theta=G^2/(G^2+ tau_vec_est(i));
%     Theta=G^2/(G^2+ taoSE);
    filter_gain=Theta*OmegaLine;
    % signal estimate
%     filter_gain=G^2/(G^2+ tau_vec_est(i)) ./ (1+ (1-lambda)/lambda * exp(M*(log(1+G^2/tau_vec_est(i)) - Delta*x_gain/M)));
    X_est=diag(filter_gain)*MF_out;
        
    % compute and update the residual
    R=Y-pilot_Mat*X_est + K/P * mean(filter_gain) * R;
       
    % signal estimate
    X=X_est;
    
    
    %stem(mean(abs(X(1:end,:)).^2,2))
    
end


end





