function [Mngh] = Monaghanvisc(selfpart,NP,Rho,cs,Dist,V1,V2,dx,dy,h)

rhoa=Rho(selfpart);
rhoij=(rhoa+Rho(NP))/2;
VIJ1=V1(NP)-V1(selfpart);
VIJ2=V2(NP)-V2(selfpart);
vijdotxij=VIJ1*dx+VIJ2*dy;
Mngh=0;
if vijdotxij<0
    cij=0.5*(cs(selfpart)+cs(NP));
    muij=(h*vijdotxij)/((Dist^2)+(0.01*(h^2)));
    piij = 0.5*cij*muij + 0.5*muij*muij;
    Mngh=  piij/(rhoij^2);
end
end

