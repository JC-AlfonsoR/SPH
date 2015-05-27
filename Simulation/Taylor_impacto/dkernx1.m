function [ dkernx ] = dkernx1( rij,hk,coorglobal,NP,selfp )
% Calcula la derivada en x del kernel para las n particulas vecinas que se
% encuentran a distnacias rij. El calculo se hace con dominio soporte de
% radio hk.
%
% [dkernx] = dkernx1(rij, hk, coorglobal, NP, selfp)
%
%   Inputs
%   rij     [1 x n]     Vector con las n distancias de las n partículas
%                       vecinas a la partícula i
%   hk      double      Radio del dominio soporte para partícula i
%   coorglobal      [npart x 2]
%                       Coordenadas x,y globales de todas las particulas
%   Np      [1 x n]     Identificación de las n partículas vecinas a la
%                       particula i
%   slefp   int         Identificación de la partícula i
%
%   Outputs
%   dkernx  [n x 1]     Vector con la magnitud de la derivada en x evaluada
%                       para las n partículas vecinas de la partícula i

dkernx = zeros(length(rij),1); %[n x 1]     vector de 0s

for i=1:length(rij)
    q = rij(i)/hk;  %q  double
    
    if q<=2
        dx = abs(coorglobal(NP(i),1)-coorglobal(selfp,1));

        fac = (7/((4*pi)*((hk)^2)));
        dval = -5*q*((1-(0.5*q))^3)/hk*dx;
    else
        dkernx(i) = 0;
    end
    
    dkernx(i) = fac*dval;
end
%% Referencias
%{
[1] G.R. Liu & M.B. Liu, Smoothed Particle Hydrodynamics - a meshfree particle
ethod, World Scientifics Publishing Co., 2003. (Tabla 3.1 - Resumen de kernels)
pg 102 - quintic smoothing function
%}