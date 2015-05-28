function [ dkerny ] = dkerny1( rij,hk,coorglobal,NP,selfp )
%Calcula la derivada en y del kernel para las n particulas vecinas que se
%encuentran a distnacias rij. El calculo se hace con dominio soporte de
%radio hk.
%
% [dkerny] = dkerny1(rij, hk, coorglobal, NP, selfp)
%
%   Inputs
%   rij     [1 x n]     Vector con las n distancias de las n part�culas
%                       vecinas a la part�cula i
%   hk      double      Radio del dominio soporte para part�cula i
%   coorglobal      [npart x 2]
%                       Coordenadas x,y globales de todas las particulas
%   Np      [1 x n]     Identificaci�n de las n part�culas vecinas a la
%                       particula i
%   slefp   int         Identificaci�n de la part�cula i
%
%   Outputs
%   dkerny  [n x 1]     Vector con la magnitud de la derivada en y evaluada
%                       para las n part�culas vecinas de la part�cula i

dkerny=zeros(length(rij),1);

for i=1:length(rij)
    q = rij(i)/hk;      %q  double
    
    if q<=2
        dy = abs(coorglobal(NP(i),2)-coorglobal(selfp,2));

        fac = (7/((4*pi)*((hk)^2)));
        dval = -5*q*((1-(0.5*q))^3)/hk*dy;
    else
        dkerny(i) = 0;
    end
    
    dkerny(i) = fac*dval;
end