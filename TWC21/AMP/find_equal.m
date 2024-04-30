function [equal_num,A1_submat,A2_submat]=find_equal(A1,A2,length)
equal_num=0;  
A1_submat=zeros(1,length);
A2_submat=zeros(1,length);
for llllll=1:length
    tempppp=find(A1==A2(llllll), 1);
    if isempty(tempppp)==0%
        equal_num=equal_num+1;       %number of correctly detected UEs
        A1_submat(equal_num)=tempppp;  %index of correctly detected UEs
        A2_submat(equal_num)=llllll;%
    end
end              
end