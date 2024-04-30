%% calculate the LLR
function [LLR_MBM]=LLR_soft_MBM1203(pf6,opt)
[m,n,p]=size(pf6);
% Normuliza=sqrt(3/2/(opt.M_order-1));
% sigma2=opt.sigma2; % noise variance of channel
LLR_MBM=zeros(2*opt.K,p);
TEmp=zeros(opt.K,opt.M,p);
for i=1:n/opt.M
%     for j=1:opt.M
        TEmp(i,1,:)=sum(pf6(2:5,opt.M*(i-1)+1,:)).*pf6(1,opt.M*(i-1)+2,:).*pf6(1,opt.M*(i-1)+3,:).*pf6(1,opt.M*(i-1)+4,:);
        TEmp(i,2,:)=sum(pf6(2:5,opt.M*(i-1)+2,:)).*pf6(1,opt.M*(i-1)+1,:).*pf6(1,opt.M*(i-1)+3,:).*pf6(1,opt.M*(i-1)+4,:);
        TEmp(i,3,:)=sum(pf6(2:5,opt.M*(i-1)+3,:)).*pf6(1,opt.M*(i-1)+2,:).*pf6(1,opt.M*(i-1)+1,:).*pf6(1,opt.M*(i-1)+4,:);
        TEmp(i,4,:)=sum(pf6(2:5,opt.M*(i-1)+4,:)).*pf6(1,opt.M*(i-1)+2,:).*pf6(1,opt.M*(i-1)+3,:).*pf6(1,opt.M*(i-1)+1,:);
%     end
end
for i=1:n/opt.M
    LLR_MBM(2*(i-1)+1,:)=log((TEmp(i,1,:)+TEmp(i,2,:))./(TEmp(i,3,:)+TEmp(i,4,:)));
    LLR_MBM(2*i,:)= log((TEmp(i,1,:)+TEmp(i,3,:))./(TEmp(i,2,:)+TEmp(i,4,:)));
end
end