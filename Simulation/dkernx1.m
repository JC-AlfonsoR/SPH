function [ dkernx ] = dkernx1( rij,hk,coorglobal,NP,selfp )
dkernx=zeros(length(rij),1);
for i=1:length(rij)
q=rij(i)/hk;
if q<=2
    dx=abs(coorglobal(NP(i),1)-coorglobal(selfp,1));
       
    fac=(7/((4*pi)*((hk)^2)));
    dval=-5*q*((1-(0.5*q))^3)/hk*dx;
else
    dkernx(i)=0;
end
dkernx(i)=fac*dval;
end