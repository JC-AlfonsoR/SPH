function [dD] = Damageevol(rho,M,C)
%Damagevol
%
%
%
%[dD] = Damageevol(rho, M, C)
%
%   Inputs
%   rho         double      Densidad de la part�cula
%   M           double      Masa de la part�cula
%   C           double      ---
%
%   Outputs
%   dD          double      Derivada del da�o
%
num = 0.4*C;
denom =((rho)/(pi()*M))^(1/2);
dD = num/denom;
end

%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentaci�n de
materiales fr�giles con el m�todo de part�culas suavizadas (SPH), Uniandes,
2015. (pg 10, Ecuaci�n 21)
%}

%% Comentarios Finales
%{
La funci�n parece diferente de la que se reporta en el documneto.
Preguntarle a Daniel.
%}