function [dD] = Damageevol(rho,M,C)
% Damagevol
%
%
%
% [dD] = Damageevol(rho, M, C)
%
%   Inputs
%   rho         double      Densidad de la partícula
%   M           double      Masa de la partícula
%   C           double      Velocidad de crecimiento de la grieta
%
%   Outputs
%   dD          double      Derivada del damage
%
num = 0.4*C; 
denom =((rho)/(pi()*M))^(1/2);
dD = num/denom;

% Calcula el crecimiento de una grieta en el delta de tiempo
end

%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentación de
materiales frágiles con el método de partículas suavizadas (SPH), Uniandes,
2015. (pg 10, Ecuación 21)
%}