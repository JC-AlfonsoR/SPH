function [V1,V2] = XSPH(Npart,M,Rho,V1,V2,kern,selfpart)


V1self=V1(selfpart);
V2self=V2(selfpart);
rhoi=Rho(selfpart);

for i=1:length(Npart)
    numpart=Npart(i);
    Wij=kern(i);
    rhoj=Rho(numpart);
    mj=M(numpart);
    fact=mj/((1/2)*(rhoj+rhoi));
    Vsumx=(V1self-V1(numpart))*Wij*fact;
    Vsumy=(V2self-V2(numpart))*Wij*fact;
 
end

V1=V1self+0.5*Vsumx;
V2=V2self+0.5*Vsumy;

end

