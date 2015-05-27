function [P,C] = EOSmie(eint,r0,C,S,rho,gamma )
% EOSmie calcula la presion hidrostatica para una particula, el calculo se
% hace con la ecuacion de estado de Mie-Gruniesen
% 
%
% [P,C] = EOSmie(eint, r0, C, S, rho, gamma)
%
%   Inputs
%   eint        double      Energia interna
%   r0          double      Densiadad inicial
%   C           double      Huggoniot
%   S           double      Huggoniot
%   rho         double      Densiad de la particula i
%   gamma       double      XSPH
%
%   Outputs
%   P           double      Presion hidrostatica
%   C           double      ~~~

a0 = r0*(C^2);
b0 = a0*(1+(2*(S-1)));
c0 = a0*((2*(S-1))+(3*((S-1)^2)));
ratio = (rho/r0) - 1.0;   %eta en ref [1]
if ratio<0
    PH = a0*ratio;    
else
    PH = a0*ratio+(ratio^2)*(b0+(c0*ratio));
end

P = ((1-0.5*gamma*ratio)*PH)+(rho*eint*gamma);
%
c2 = (1/r0)*((a0*(1-gamma*ratio)) +...
    b0*(2*ratio-1.5*(gamma*ratio^2))+...
    c0*(3*ratio^2-2*gamma*ratio^3)+...
    gamma*eint*(1+ratio) + gamma*P*eint*(1+ratio));
C = c2^0.5;

end

%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentacion de
materiales fragiles con el metodo de particulas suavizadas (SPH), Uniandes,
2015, pg 9 (Eucacion de estado);

[2] G.R. Liu & M.B. Liu, Smoothed Particle Hydrodynamics - a meshfree particle
ethod, World Scientifics Publishing Co., 2003. (Eq 8.9)
pg 313
%}
%}