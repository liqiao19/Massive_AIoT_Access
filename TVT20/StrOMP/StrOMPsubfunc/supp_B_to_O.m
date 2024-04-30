function [ T_index ] = supp_B_to_O( candidate_set,M)
%SUPP_B_TO_O Summary of this function goes here
%   Detailed explanation goes here
lengg=length(candidate_set);
T_index=zeros(lengg,M);
for qqqq=1:lengg
    T_index(qqqq,:)=((candidate_set(qqqq)-1)*M+1):(candidate_set(qqqq)*M);
end
T_index=reshape( T_index,1,[]);
T_index=sort(T_index);
end

