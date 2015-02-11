%% BulletImpact_2_data simmulation
% Recrea los datos simulados de BulletImpact_2
% 
%
% La simulación BulletImpact_data fue corrida con el script BulletImpact_2.m
% con las sigueintes especificaciones:
%
%
%
%
%
%
%% Importar datos
clear;close all; clc;
% Abrir el BulletImpact_data.mat y extraer la información en el arreglo
% correspondiente.
%       Todos los datos son arreglos de orden 3, con 1407 datos que
%       corresponden a las 1407 particulas de la simulación. Todos los
%       arreglos tienen tamaño 200 en la dirección 3, que corresponden a
%       los 200 steps de tiempo en la simulación.
   
BulletImpact = open('BulletImpact_2_data.mat');
Coord = BulletImpact.Coordenadas;
V1 = BulletImpact.Velocidad1;
V2 = BulletImpact.Velocidad2;
Dev11 = BulletImpact.Esfuerzos11;
Dev12 = BulletImpact.Esfuerzos12;
Dev21 = BulletImpact.Esfuerzos21;
Dev22 = BulletImpact.Esfuerzos22;
Rho = BulletImpact.Densidad;
P = BulletImpact.Presion;

%% Calculos
V = sqrt(V1.^2+V2.^2);
%% Graficar Datos
steps = 174;
% The step #175 aborted the simulation


for i=1:2
   figure(1)
   %plot(Coord(:,1,i), Coord(:,2,i),'.b')
   scatter(Coord(:,1,i), Coord(:,2,i),10,V(:,1,i),'filled')
   axis([-0.015 0.015 -0.01 0.03])
   title('Simulacion BulletImpact.mat - V [m/s]')
   xlabel('x  [m]')
   ylabel('y [m]')
   colorbar
   drawnow
end
%% Comentarios finales
%{
% Febrero 10
Existen problemas con el colorbar porque se cambia en cada paso de tiempo.
Definir un colorbar constante para que todos los pasos de tiempo se
grafiquen con el mismo colorbar
%}