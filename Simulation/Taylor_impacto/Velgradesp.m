function [dv1dx1,dv1dx2,dv2dx1,dv2dx2 ] = Velgradesp(rho,m,V1,V2,dwx,dwy,Npart,Spart)
% Usando aproximacion por las particulas vecinas, se calculan las derivadas
% espaciales para las componentes de velocidad V1 y V2 en la particula
% Spart
%
% [dv1dx1, dv1dx2, dv2dx1, dv2dx2] = 
%                   Velgradesp(rho, m, V1, V2, dwx, dwy, Npart, Spart)    
%
%   Inputs
%   rho     [npart x 1]     Densidad de todas las particulas
%   m       [npart x 1]     Masa de todas las particulas
%   V1      [npart x 1]     Velocidad en direccion 1
%   V2      [npart x 1]     Velocidad en direccion 2
%   dwx     [n x 1]         Derivada de kernel en x para las n particulas
%                               vecinas
%   dwy     [n x 1]         Derivada de kernel en y para las n particulas
%                               vecinas
%   Npart   [1 x n]         Identidad de las n particulas vecinas
%   Spart   int             Identidad de la particula actual
%
%   Outputs
%   dv1dx1      double      Derivada de V1 en dir 1 para particula Spart
%   dv1dx2      double      Derivada de V1 en dir 2 para particula Spart
%   dv2dx1      double      Derivada de V2 en dir 1 para particula Spart
%   dv2dx2      double      Derivada de V2 en dir 2 para particula Spart

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
