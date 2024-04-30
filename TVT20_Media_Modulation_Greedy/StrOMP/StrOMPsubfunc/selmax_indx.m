function [ indx ] = selmax_indx( abs_b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
xxx=reshape(abs_b,4,[]);
% xxx=reshape([4 3 2 1 6 7 9 10 2 4 5 1]',4,[]);
[yyy,indd]=sort(xxx,'descend');
indx=indd(1,:);
for iii=1:size(indd,2)
    indx(iii)=indx(iii)+4*(iii-1);
end

end

