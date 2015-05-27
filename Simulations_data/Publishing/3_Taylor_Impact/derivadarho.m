function [drho] = derivadarho(M,V1,V2,dWx,dWy,Npart,selfpart)

V1self=V1(selfpart);
V2self=V2(selfpart);
sum=0;
for i=1:length(Npart)
    numpart=Npart(i);
    dWi=[dWx(i) dWy(i)];    
    Vsumx=(V1self-V1(numpart))*dWi(1);
    Vsumy=(V2self-V2(numpart))*dWi(2);
    Vsum=Vsumx+Vsumy;
    sum= sum+(M(numpart)*Vsum);
end
drho=sum;
end
