%% Bullet_Impact
%
% Simulacion de impacto
%
% Version 1
%
% Marzo 18 - 2015
%
% Author: *J. Camilo Alfonso R.*
%
% *Problema Especial IMEC*
%
% Simulacion de impacto entre proyectil ductil y objetivo fragil
%
% _Se omiten tildes para evitar problemas de compatibilidad en ecoding_
clear all; clc;
%% Definir Geometria del Objetivo
% Todas las unidades son dadas en el sistema internacional de unidades
%
% Geometria del objetivo
dy = 3e-4; %    Separacion entre particulas
dx = dy;   %     
k = 2.0;   %    Constante para expandir radio de soporte
h = k*dx;  %    Radio de soporte

T_width = 0.0076;               % Ancho del objetivo
T_height = 0.0076;              % Alto del objetivo
T_x = -T_width : dx : T_width;
T_y = -T_height : dy : T_height;
[X,Y] = meshgrid(T_x, T_y);     % Matriz con la malla de las posiciones x,y para las particulas
Target = [X(:),Y(:)];           %   | x_1 , y_1 |   En esta matriz organiza    
                                %   | .   , .   |   todas las posiciones de las
                                %   | x_n , y_n |   particulas.
                                %   [T_np x 2]

T_np = size(Target,1);          % Numero de particulas en el objetivo

%% Asignacion de defectos puntuales en Objetivo
% Constantes para asginacion de fallas en basalto
m = 3;
k = 7;
V = dx*dy*1;    % Volumen infinitesimal

Nflaws = T_np*log(T_np);        % Numero de defectos puntuales a asignar
Nflaws = round(Nflaws);         
assign_flaws = randi(T_np,Nflaws,1,'uint32'); % [Nflaws x 1]
                                %  Nflaws numeros aleatorios entre 1 y T_np
[Flaws{1:T_np,1}] = deal([]);     %  Flaws [cell array] {T_np x 1}
                                %  Cell array con ''T_np matrices vacias
                                %  Cell de pre-alocacion para fallas
                                
for i = 1:Nflaws
    Flaws{assign_flaws(i),1}(size(Flaws{assign_flaws(i),1})+1) = ...
        (i/(k*V))^(1/m);
end
%%
% *assign_flaws*   [Nflaws x 1]    vector que contiene _Nflaws_ posiciones de
% las fallas a generar. Las posiciones estan dadas como numeros enteros en
% relacion al numero de cada particula.
% 
% *Flaws*       cell{T_np x 1}    Contiene _T_np_ matrices vacias que
% identifican las fallas para cada particula. *Flaws* se va llenando 
% de forma aleatoria con las posiciones que indica _assing_flaws_. La
% primera vez que se pasa por una amtriz de Flaws, le asigna el numero
% (i/(k*V))^(1/m) formando una matriz 1x1. La segunda vez que se pasa por
% la misma matriz, aumenta la dimension de la matriz en una sola direccion
% para asignar otro numero. Asi, si semi-aleatoriamente, el numero _k_
% aparecio _n_ veces de _assign_flaws_, la celda _Flaws_ en su posicion _k_
% debe contener una matriz nx1 con numeros asignados. Los numeros asignados
% corresponden a las deformaciones de activacion para los defectos
% puntuales de cada particula
