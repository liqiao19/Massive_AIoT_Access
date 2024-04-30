function [G2,Hk] = channel_model_zhangshun_v2_TV(f_d,N_BS,N,fs,Lp,F,A_BS,Power0,MaxSym,Delay)
% 相比于v1，多了传递的角度信息
% generate the uplink frequency selective fading channel, the BS is
% equipped with ULA while single-antenna is considered at user device.

% written by LI QIAO (qiaoli@bit.edu.cn), Beijing Institute of Technology
% version: 2021.07.04

% Inputs：
%   N_BS：number of BS antennas
%   N: number of subcarriers
%   P: pilot subcarriers
%   Bw: bandwidth
%   fs: frequency of sampling
%   fc: frequency of carrier
%   L_cp: length of cyclic prefix
%   sigma2_alpha: average power of path, generally, sigma2_alpha = 1;
%   Lp：number of path

% Outputs：
%   H_f：frequency domain channel matrix for pilots, N_BS×P
%   H_ang_f：virtual-angular domain channel matrix for pilots, N_BS×P

%% 
% f_d=opt.f_d;
% N_BS=opt.N_BS;
% N=opt.N ;
% fs=opt.fs ;
% Lp=opt.Lp;                      % number of multipath % different for diffrent devices && different for different OFDM symbols
% PowerdB=[0,-1*sort(randperm(30,Lp-1))];
% Power0=10.^(PowerdB/10);
% Power0=Power0/sum(Power0);
%% Generate the AoAs and steering vectors
% n_BS = (0:N_BS-1).';	
% index_BS = randi(N_BS,[1,1]);
% index_BS = index_BS:index_BS+Lp-1; 
% index_BS(index_BS>N_BS) = index_BS(index_BS>N_BS) - N_BS;
% theta_BS = index_BS/N_BS;
% A_BS = exp(-1i*2*pi*n_BS*theta_BS);
A_BS=reshape(A_BS,N_BS,Lp);
Power_eff=repmat(Power0,N_BS,1);
Power_eff=Power_eff.*A_BS;
%% Doppler
% H_eff3=zeros(N,N,N_BS);
G2=zeros(N*MaxSym,N,N_BS);
theta0=1:N*MaxSym;

F_d=f_d*(rand(Lp,1)-0.5)*2; 
theta=exp(1i*2*pi*F_d*theta0/fs);
Hk=zeros(N,N_BS);
%————————多径延时 生成——————————%
% N_cp=16;32;64;
% Delay = sort(randperm(N_cp-1,Lp-1)+1);
%————————多径延时 生成 END——————————%

for j=1:N_BS
    Htemp=zeros(N,MaxSym);
    for s=1:MaxSym
        Power=Power_eff(j,:);
        TEMP=zeros(N,N,Lp);
        x00=zeros(N,1);
        x00(1)=Power(1);
        thetaTemp=repmat(theta(1,(s-1)*N+1:s*N),N,1);
        TEMP(:,:,1)=thetaTemp.*toeplitz(x00,x00.');
        for i=2:Lp
            ww=Delay(i-1); %% 找到多径对应的延时域坐标
            x0=zeros(N,1);
            y0=zeros(1,N);
            x0(ww)=Power(i);
            y0(end-ww+2)=Power(i);
            thetaTemp=repmat(theta(i,(s-1)*N+1:s*N),N,1);
            TEMP(:,:,i)=thetaTemp.*toeplitz(x0,y0);
        end
        TEMPPP=sum(TEMP,3);
        Htemp(:,s)=diag(TEMPPP);
        G2((s-1)*N+1:s*N,:,j)=F*TEMPPP*F';
    end
    Hk(:,j)=mean(Htemp,2);
end
%% citural anglar domain
% H_ang_f=G2;
end
