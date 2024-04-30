function [ACT_SETTT]=AUD_thresholdParfor(K,AUDfactors,TH)
% K=opt.K;
TempFactors=zeros(1,K);
TempFactors(AUDfactors>TH)=1;
TempFactors(AUDfactors<=TH)=0;
ACT_SETTT=find(TempFactors==1);
end