%X. Ma, J. Kim, D. Yuan and H. Liu, ¡°Two-level sparse structure based compressive sensing detector for 
%uplink spatial modulation with massive connectivity,¡± IEEE Commun. Lett., vol. 23, no. 9, pp. 1594-1597, Sept. 2019.
function [ pos_num_csmp,remod_ramp,rant_ramp] =TLSSCS(orthspany,y,H_wholeuser,imax, US_toant,T,Vthnorm)
[BS_ant,Tol_ant] = size(H_wholeuser);
pos_num_csmp = [];
orthprojec=eye(BS_ant); %initialization orthogonal projector
res_osmp=y;
for kk=1:imax
    % step3 orthogonal projector
Hnew=orthprojec*H_wholeuser;
% step4 User support set selection
inorm=cell2mat(arrayfun(@(x) norm(orthspany'*Hnew(:,(x-1)*US_toant+1:x*US_toant),'fro')/norm(Hnew(:,(x-1)*US_toant+1:x*US_toant),'fro'),1:size(H_wholeuser,2)/US_toant,'UniformOutput',false));
inorm(:,pos_num_csmp)=0;
[val,pos]=sort(inorm,'descend');
Js = pos(1:1);
% step5 Support Merger
pos_num_previous=pos_num_csmp;
pos_num_csmp = union(pos_num_csmp,Js);
% step6,7 Estimation
At=cell2mat(arrayfun(@(x) H_wholeuser(:,(pos_num_csmp(x)-1)*US_toant+1:pos_num_csmp(x)*US_toant),1:length(pos_num_csmp),'UniformOutput',false));
orthprojec=diag(ones(BS_ant,1))-At*((At'*At)\At'); % step7 Orthogonal
% projector
% step8 Residue update
previousres=res_osmp;
res_osmp = orthprojec*y;
% step9 10 11 12 terminating condition, adaptive sparsity level
minnorm=norm(previousres,'fro')^2-norm(res_osmp,'fro')^2;
if minnorm<Vthnorm
pos_num_csmp=pos_num_previous; %step14
break;
end
end
% step 15
if isempty(pos_num_csmp)==1
    remod_ramp=0;
    rant_ramp=0;
else
    Af=cell2mat(arrayfun(@(x) H_wholeuser(:,(pos_num_csmp(x)-1)*US_toant+1:pos_num_csmp(x)*US_toant),1:length(pos_num_csmp),'UniformOutput',false));
    theta_ls=Af'*Af\Af'*y;
    remod_ramp=zeros(length(pos_num_csmp),T);
    rant_ramp=zeros(length(pos_num_csmp),T);
    for det=1:T
        countnum=0;
        for i=1:length(pos_num_csmp)
            xeacht=abs(theta_ls(countnum+(1:US_toant),det)).^2;
            % step 16
            % if max(xeacht)<=pth
            % remod_ramp(i,det)=0;
            % rant_ramp(i,det)=0; % step 18
            % else
            % step 17 Prune support set
            [devalue,depos]=sort(xeacht,'descend');
            remod_ramp(i,det)=theta_ls(countnum+depos(1),det);
            % rant_ramp(i,det)=depos(1)+(pos_num_csmp(i)-1)*US_toant;
            rant_ramp(i,det)=depos(1)-1;
            % end
            countnum=US_toant+countnum;
        end
    end
    pos_num_csmp_ant=unique(reshape(rant_ramp,1,size(rant_ramp,1)*size(rant_ramp,2)));
    pos_num_csmp_ant(find(pos_num_csmp_ant==0))=[];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%End of Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
