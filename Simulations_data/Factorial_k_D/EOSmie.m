function [P,C] = EOSmie(eint,r0,C,S,rho,gamma )
% EOSmie calcula la presi�n hidrost�tica para una part�cula, el calculo se
% hace relacionando las dem�s variables de estado
% de estado.
% 
%
% [P,C] = EOSmie(eint, r0, C, S, rho, gamma)
%
%   Inputs
%   eint        double      Energ�a interna
%   r0          double      Densiadad inicial
%   C           double      Huggoniot
%   S           double      Huggoniot
%   rho         double      Densiad de la part�cula i
%   gamma       double      XSPH
%
%   Outputs
%   P           double      Presi�n hidrostatica
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
% Estas lineas existen en el archivo EOSmie[conflicto].m, pero no existen 
% en el EOSmie.m - Preungtarle a Daniel por estas lineas.
%}
end

%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentaci�n de
materiales fr�giles con el m�todo de part�culas suavizadas (SPH), Uniandes,
2015, pg 9 (Eucaci�n de estado);
%}

%% Comentarios Finales
%{
No estoy seguro de que la funci�n trabaje bien.Cuando se llama
desde el main se espera una sola salida, pero la funci�n tiene 2 salidas.
%}