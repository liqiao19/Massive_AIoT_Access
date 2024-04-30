function [EstSymbol]=DecodingParfor(GSMset,ESTindex)   
% GSMset=opt.GSMset;%与6选2获得4bit有区别
    [m,n]=size(ESTindex);
    EstSymbol=zeros(m,1);
    for i=1:m
        temp=0;
        for j=1:n
            temp= temp+abs(ESTindex(i,j)-GSMset(:,j));
        end
        EstSymbol(i)=find(temp==0);
    end
end