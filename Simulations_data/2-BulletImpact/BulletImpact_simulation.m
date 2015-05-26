%% BulletImpact_data simmulation
% Recrea los datos simulados de BulletImpact
% Se requiere que ela archivo tr_data.mat este en la misma carpeta.
% trdata disponible en: https://github.com/JC-AlfonsoR/SPH/raw/master/...
%                               ...Simulations_data/Bulletimpact_data.mat
%
% La simulaci�n BulletImpact_data fue corrida con el script BulletImpact.m
% con las sigueintes especificaciones:
%
%
%
%
%
%
%% Importar datos
clear;close all; clc;
% Abrir el BulletImpact_data.mat y extraer la informaci�n en el arreglo
% correspondiente.
%       Todos los datos son arreglos de orden 3, con 1407 datos que
%       corresponden a las 1407 particulas de la simulaci�n. Todos los
%       arreglos tienen tama�o 200 en la direcci�n 3, que corresponden a
%       los 200 steps de tiempo en la simulaci�n.
   
BulletImpact = open('Bulletimpact_data.mat');
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



   figure(1)
   subplot(1,2,1)
   scatter(Coord(:,1,1), Coord(:,2,1),10,V(:,1,1),'filled')
   axis([-0.01 0.01 0 0.025])
   title('Inicial')
   xlabel('x  [m]')
   ylabel('y [m]')
   colorbar
   
   subplot(1,2,2)
   scatter(Coord(:,1,steps), Coord(:,2,steps),10,V(:,1,steps),'filled')
   axis([-0.01 0.01 0 0.025])
   title('Final')
   xlabel('x  [m]')
   ylabel('y [m]')
   colorbar


%% Comentarios finales
%{
% Febrero 10
Existen problemas con el colorbar porque se cambia en cada paso de tiempo.
Definir un colorbar constante para que todos los pasos de tiempo se
grafiquen con el mismo colorbar
%}