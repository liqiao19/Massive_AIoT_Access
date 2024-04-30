function [equal_num,A1_submat,A2_submat]=find_equal(A1,A2,length)
equal_num=0;  %每次都从零开始加，置零放在循环里面
A1_submat=zeros(1,length);
A2_submat=zeros(1,length);
for llllll=1:length
    tempppp=find(A1==A2(llllll), 1);
    if isempty(tempppp)==0%如果活跃用户被检测到了
        equal_num=equal_num+1;       %检正确检测出的用户个数
        A1_submat(equal_num)=tempppp;  %检测正确用户在检测用户集中的下标
        A2_submat(equal_num)=llllll;%
    end
end              
end