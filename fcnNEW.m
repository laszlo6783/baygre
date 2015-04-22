function [f, IM,IMall]=fcnNEW(xx,b,LL,IM,IMall)
[~,cx]=size(xx);
xx=xx(2*LL-1:2*LL,:);
mx=reshape(xx',2*cx,1);
mm=zeros(2*cx,2);
i=1;
while i<=2*cx;
    if mx(i)==1;mm(i,1)=1;mm(i,2)=0;
    elseif mx(i)==2;mm(i,1)=0;mm(i,2)=1;
    else mm(i,1)=-1;mm(i,2)=-1;
    end
    i=i+1;
end
mm1=reshapeG(mm,2,2*cx);

		% !!! mixed logit !!!
minus=IMall(:,8*(LL-1)+1:8*LL);
plus=im(mm1,b);
IM=IM+plus-minus;

IMall(:,8*(LL-1)+1:8*LL)=plus;

s=IM;
%disp(s);
d=det(s);
if d<=0
    f=1000;
    return
else
    f=(d^(-1/(2*cx)));
    return
end
