function [kern] = kern1(rij,hk)
kern=zeros(length(rij),1);
for i=1:length(rij)
q=rij(i)/hk;
if q<=2
    fac=(7/((4*pi)*((hk)^2)));
    val=((1-(0.5*q))^4)*((2*q)+1);
    kern(i)= val*fac;
else
    kern(i)=0;
end
end

