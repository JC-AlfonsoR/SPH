%% tr_data simmulation
% Recrea los datos simulados de tr_data.mat
% Se requiere que ela archivo tr_data.mat este en la misma carpeta.
% trdata disponible en: https://github.com/JC-AlfonsoR/SPH/raw/master...
%                                       .../Simulations_data/tr_data.mat
%
% La simulación tr_data.mat fue corrida con el script tr.m con las
% sigueintes especificaciones:
%
%
%
%
%
%
%% Importar datos
% Abrir el archivo tr_data.mat y extraer la información en el arreglo
% correspondiente.
%       Todos los datos son arreglos de orden 3, con 1407 datos que
%       corresponden a las 1407 particulas de la simulación. Todos los
%       arreglos tienen tamaño 200 en la dirección 3, que corresponden a
%       los 200 steps de tiempo en la simulación.
   
tr_data = open('tr_data.mat');
Coord = tr_data.Coordenadas;
V1 = tr_data.Velocidad1;
V2 = tr_data.Velocidad2;
Dev11 = tr_data.Esfuerzos11;
Dev12 = tr_data.Esfuerzos12;
Dev21 = tr_data.Esfuerzos21;
Dev22 = tr_data.Esfuerzos22;
Rho = tr_data.Densidad;
P = tr_data.Presion;

%% Graficar Datos
steps = 174;
% The step #175 aborted the simulation


for i=1:steps
   figure(1)
   %plot(Coord(:,1,i), Coord(:,2,i),'.b')
   scatter(Coord(:,1,i), Coord(:,2,i),10,P(:,1,i),'filled')
   axis([-0.02 0.02 -0.02 0.04])
   title('Simulacion tr-data.mat - Presion Hidrostatica')
   xlabel('x  [m]')
   ylabel('y [m]')
   colorbar
   drawnow
end