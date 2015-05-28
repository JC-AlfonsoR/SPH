function [dv1dx1,dv1dx2,dv2dx1,dv2dx2 ] = Velgradesp(rho,m,V1,V2,dwx,dwy,Npart,Spart)
% Usando aproximación por las particulas vecinas, se calculan las derivadas
% espaciales para las componentes de velocidad V1 y V2 en la partícula
% Spart
%
% [dv1dx1, dv1dx2, dv2dx1, dv2dx2] = 
%                   Velgradesp(rho, m, V1, V2, dwx, dwy, Npart, Spart)    
%
%   Inputs
%   rho     [npart x 1]     Densidad de todas las partículas
%   m       [npart x 1]     Masa de todas las partículas
%   V1      [npart x 1]     Velocidad en dirección 1
%   V2      [npart x 1]     Velocidad en dirección 2
%   dwx     [n x 1]         Derivada de kernel en x para las n partículas
%                               vecinas
%   dwy     [n x 1]         Derivada de kernel en y para las n partículas
%                               vecinas
%   Npart   [1 x n]         Identidad de las n partículas vecinas
%   Spart   int             Identidad de la partícula actual
%
%   Outputs
%   dv1dx1      double      Derivada de V1 en dir 1 para partícula Spart
%   dv1dx2      double      Derivada de V1 en dir 2 para partícula Spart
%   dv2dx1      double      Derivada de V2 en dir 1 para partícula Spart
%   dv2dx2      double      Derivada de V2 en dir 2 para partícula Spart

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
Preguntarle a Daniel por esta función, al parecer esta diferente a la
ecuación del libro.
En los terminos dv1 y dv2 se calcula una diferencia entre las velocidades
de la aprtícula i y cada una de las partículas vecinas j. El dv1 y el dv2
se multiplican por el kernel y por el factor de masa/densidad.
En la ecuación del libro, la única diferencia es que no se calculan las
diferencias dv1 y dv2 sino que se aproximan con la velocidad evaluada en
cada una de las partículas.
%}