function [drho] = derivadarho(M,V1,V2,dWx,dWy,Npart,selfpart)
% derivadarho
%
% Aproxima la derivada de la densidad en la particula selfpart a traves del
% método de SPH
%
%[drho] = derivadarho(M, V1, V2, dWx, dWy, Npart, selfpart)
%
%   Inputs
%   M           [npart x 1]     Masa de las particulas
%   V1          [npart x 1]     Velocidad en 1 de las particulas
%   V2          [npart x 1]     Velocidad en 2 de las particulas
%   dWx         [1 x n]         Derivada de kernel en x, evaluada para las
%                             n particulas vecinas de la particula selfpart
%   dWy         [1 x n]         Derivada de kernel en y, evaluada para las
%                             n particulas vecinas de la aprticula selfpart
%   Npart       [1 x n]         Identidad de las n particulas vecinas de la
%                               particula selpart
%   selfpart    int             Identidad de la partícula actual
%
%   Output
%   drho        double          Derivada de la densidad en la particula
%   actual
%
%

V1self = V1(selfpart); %    double      velocidad 1 de la particula actual
V2self = V2(selfpart); %    double      velocidad 2 de la particula actual
sum = 0;
for i=1:length(Npart)
    numpart = Npart(i);
    dWi = [dWx(i) dWy(i)];    
    Vsumx = (V1self-V1(numpart))*dWi(1);
    Vsumy = (V2self-V2(numpart))*dWi(2);
    Vsum = Vsumx + Vsumy;
    sum = sum +(M(numpart)*Vsum);
end
drho = sum;
end
%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentación de
materiales frágiles con el método de partículas suavizadas (SPH), Uniandes,
2015. (Ecuacion 22 - Ecuación de continuidad discretizada)
%}
%% Comentarios
%{
---
%}