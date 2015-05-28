function [dV1,dV2] = Momentumeq2d(dev11,dev12,dev21,dev22,P,rho,dkernx,dkerny,m,Npart,selfpart,cs,Dist,coorglobal,h,V1,V2)
%Momentumeq2d
%
% Resuelve la ecuaci�n de Momentum para particula selfpart. Calcula la 
% derivadad total de las velocidades. Emplea aproximaci�n por part�culas.
%
%
%   Inputs
%   dev11       [npart x 1]     Esfuerzos 11
%   dev12       [npart x 1]     Esfuerzos 12
%   dev21       [npart x 1]     Esfuerzos 21
%   dev22       [npart x 1]     Esfuerzos 22
%   P           [npart x 1]     Presion Hidrostatica    
%   rho         [npart x 1]     Densidad 
%   dkernx      [1 x n]         Derivada x de kernel evaluada para las n
%                               particulas vecinas
%   dkerny      [1 x n]         Derivada y de kernel evaluada para las n
%                               particulas vecinas
%   m           [npart x 1]     Masa de las particulas
%   Npart       [1 x n]         Identidad de las n particulas vecinas
%   selpart     int             Identidad de la particula actual
%   cs          [npart x 1]     ~~~
%   Dist        [1 x n]         Distancia de las n particulas vecinas
%   coorglobal  [npart x 2]     Coordenadas globales XY de las particulas
%   h           double          radio para dominio de soporte
%   V1          [npart x 1]     Velocidad en 1 de las particulas
%   V2          [npart x 1]     Velocidad en 2 de las particulas
%
%   Outputs
%   dV1         double          Derivada temporal de V1
%   dV2         double          Derivada temporal de V2
%
rhoa = rho(selfpart);
rhoa21 = 1/(rhoa*rhoa);
sigma11a = P(selfpart) + dev11(selfpart);
sigma12a = dev12(selfpart);
sigma21a = dev21(selfpart);
sigma22a = P(selfpart) + dev22(selfpart);
dV1 = 0;
dV2 = 0;
for i=1:length(Npart)
    NP = Npart(i);
    rhob = rho(NP);
    mb = m(NP);
    sigma11b = P(NP) + dev11(NP);
    sigma12b = dev12(NP);
    sigma21b = dev21(NP);
    sigma22b = P(NP) + dev22(NP);
    rhob21 = 1/(rhob*rhob);
    dW = [dkernx dkerny];
    dx = abs(coorglobal(Npart(i),1) - coorglobal(selfpart,1));
    dy = abs(coorglobal(Npart(i),2) - coorglobal(selfpart,2));
    Mgh = Monaghanvisc(selfpart, NP, rho, cs, Dist(i), V1, V2, dx, dy, h);
    dV1 = (dV1 + mb*(sigma11a*rhoa21+sigma11b*rhob21+Mgh)*dW(1) +...
        mb*(sigma12a*rhoa21+sigma12b*rhob21)*dW(2));
    dV2 = (dV2 + mb*(sigma21a*rhoa21+sigma21b*rhob21)*dW(1) +...
        mb*(sigma22a*rhoa21+sigma22b*rhob21+Mgh)*dW(2));
end
%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentaci�n de
materiales fr�giles con el m�todo de part�culas suavizadas (SPH), Uniandes,
2015. (Ecuacion 23 - Conservacion de momentum dsicretizada)

[2] G.R. Liu & M.B. Liu, Smoothed Particle Hydrodynamics - a meshfree particle
ethod, World Scientifics Publishing Co., 2003. (Ecuacion 4.41 - 
Conservacion de Momentum)
%}