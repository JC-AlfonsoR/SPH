function [dD] = Damageevol(rho,M,C)

num=0.4*C;
denom=((rho)/(pi()*M))^(1/2);
dD=num/denom;

end

