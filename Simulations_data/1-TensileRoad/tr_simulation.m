%% tr_data simmulation
% Recrea los datos simulados de tr_data.mat
% Se requiere que ela archivo tr_data.mat este en la misma carpeta.
% trdata disponible en: https://github.com/JC-AlfonsoR/SPH/raw/master/Simulations_data/tr_data.mat
%% Importar datos
% Abrir el archivo tr_data.mat y extraer la información en el arreglo
% correspondiente.
%       Todos los datos son arreglos de dimension 3, con 1407 datos que
%       corresponden a las 1407 particulas de la simulación. Todos los
%       arreglos tienen tamano 200 en la direccion 3, que corresponden a
%       los 200 steps de tiempo en la simulacion.
   
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



   figure(1)
   subplot(1,2,1)
   scatter(Coord(:,1,1), Coord(:,2,1),10,V2(:,1,1),'filled')
   axis([-0.02 0.02 -0.02 0.04])
   title('Velocidad Inicial')
   xlabel('x  [m]')
   ylabel('y [m]')
   
   axis([-0.01,0.01,-0.01,0.04])
   colorbar
    
   subplot(1,2,2)
   scatter(Coord(:,1,steps), Coord(:,2,steps),10,V2(:,1,steps),'filled')
   axis([-0.02 0.02 -0.02 0.04])
   title('Velocidad Final')
   xlabel('x  [m]')
   ylabel('y [m]')
   caxis([-1500,1500])
   axis([-0.01,0.01,-0.01,0.04])
   colorbar
   