function [ dkerny ] = dkerny1( rij,hk,coorglobal,NP,selfp )
dkerny=zeros(length(rij),1);
for i=1:length(rij)
q=rij(i)/hk;
if q<=2
    dx=abs(coorglobal(NP(i),2)-coorglobal(selfp,2));
       
    fac=(7/((4*pi)*((hk)^2)));
    dval=-5*q*((1-(0.5*q))^3)/hk*dx;
else
    dkerny(i)=0;
end
dkerny(i)=fac*dval;
end