function [P,C] = EOSmie(eint,r0,C,S,rho,gamma )
a0=r0*(C^2);
b0=a0*(1+(2*(S-1)));
c0=a0*((2*(S-1))+(3*((S-1)^2)));
ratio = (rho/r0) - 1.0;
if ratio<0
    PH=a0*ratio;    
else
    PH=a0*ratio+(ratio^2)*(b0+(c0*ratio));
end

P=((1-0.5*gamma*ratio)*PH)+(rho*eint*gamma);
end

