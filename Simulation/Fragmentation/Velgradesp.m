function [dv1dx1,dv1dx2,dv2dx1,dv2dx2 ] = Velgradesp(rho,m,V1,V2,dwx,dwy,Npart,Spart)
% Usando aproximaci�n por las particulas vecinas, se calculan las derivadas
% espaciales para las componentes de velocidad V1 y V2 en la part�cula
% Spart
%
% [dv1dx1, dv1dx2, dv2dx1, dv2dx2] = 
%                   Velgradesp(rho, m, V1, V2, dwx, dwy, Npart, Spart)    
%
%   Inputs
%   rho     [npart x 1]     Densidad de todas las part�culas
%   m       [npart x 1]     Masa de todas las part�culas
%   V1      [npart x 1]     Velocidad en direcci�n 1
%   V2      [npart x 1]     Velocidad en direcci�n 2
%   dwx     [n x 1]         Derivada de kernel en x para las n part�culas
%                               vecinas
%   dwy     [n x 1]         Derivada de kernel en y para las n part�culas
%                               vecinas
%   Npart   [1 x n]         Identidad de las n part�culas vecinas
%   Spart   int             Identidad de la part�cula actual
%
%   Outputs
%   dv1dx1      double      Derivada de V1 en dir 1 para part�cula Spart
%   dv1dx2      double      Derivada de V1 en dir 2 para part�cula Spart
%   dv2dx1      double      Derivada de V2 en dir 1 para part�cula Spart
%   dv2dx2      double      Derivada de V2 en dir 2 para part�cula Spart

%% Calculos

% Incializar las salidas
dv1dx1 = 0;
dv1dx2 = 0;
dv2dx1 = 0;
dv2dx2 = 0;

for i=1:length(Npart);
    fact = m(Npart(i))/rho(Npart(i));
    dv1 = V1(Npart(i));%V1(Spart)-V1(Npart(i));
    dv2 = V1(Npart(i));%V2(Spart)-V2(Npart(i));
    
    dx1 = dwx(i);
    dx2 = dwy(i);
    
    dv1dx1 = dv1dx1 + (dv1*dx1*fact);
    dv1dx2 = dv1dx2 + (dv1*dx2*fact);
    dv2dx1 = dv2dx1 + (dv2*dx1*fact);
    dv2dx2 = dv2dx2 + (dv2*dx2*fact);

end

%% Referencias
%{
G.R. Liu & M.B. Liu, Smoothed Particle Hydrodynamics - a meshfree particle
ethod, World Scientifics Publishing Co., 2003. (pg 44, eq 2.25)
%}
%% Comentarios Finales
%{
Preguntarle a Daniel por esta funci�n, al parecer esta diferente a la
ecuaci�n del libro.
En los terminos dv1 y dv2 se calcula una diferencia entre las velocidades
de la aprt�cula i y cada una de las part�culas vecinas j. El dv1 y el dv2
se multiplican por el kernel y por el factor de masa/densidad.
En la ecuaci�n del libro, la �nica diferencia es que no se calculan las
diferencias dv1 y dv2 sino que se aproximan con la velocidad evaluada en
cada una de las part�culas.
%}