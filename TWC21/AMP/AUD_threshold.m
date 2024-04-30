function [ACT_SETTT]=AUD_threshold(opt,AUDfactors,TH)
TempFactors=zeros(1,opt.K);
TempFactors(AUDfactors>TH)=1;
TempFactors(AUDfactors<=TH)=0;
ACT_SETTT=find(TempFactors==1);
end