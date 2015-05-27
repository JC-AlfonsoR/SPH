%% Taylor_impact
% Con este script se busca simular el experimento de la barra de taylor.
% El experimento consiste en impactar una barra ductil contra una
% superficie infinitamente rigida
clear all; close all; clc
% Todas las unidades en el sistema internacional de unidades
%% Geometria

dy = 0.000384848;   %   [1x1] separacion entre particulas
dx = dy;
hdx = 2.0;          %   [#] cte para expandir radio
h = hdx*dx;         %   [#] Radio para dominio soporte

bar_width = 0.0076; %   [#] Ancho de la barra
xbar = -bar_width/2:dx:bar_width/2+dx;  % [#:dx:#] posiciones en x
ybar = 0:dy:0.0254;         %   [#:dy:#]    posiciones en y
[X,Y] = meshgrid(xbar,ybar);%   [Matriz]    Matrices que con la malla 
                            %       de las posiciones para las particulas

coorbar = [X(:),Y(:)];      %   | x_1 , y_1 |   En esta matriz organiza    
                            %   | .   , .   |   todas las posiciones de las
                            %   | x_n , y_n |   particulas.
                            %   [npart x 2]

num_bar = length(coorbar(:,1));
% Particulas virtuales de frontera
frontera_x = -bar_width:dx:bar_width+dx;
frontera_y = -2*dx;
[X,Y] = meshgrid(frontera_x,frontera_y);
xx = [coorbar(:,1);X(:)];
yy = [coorbar(:,2);Y(:)];
coorbar = [xx(:),yy(:)];

%% Propiedades Fisicas

r0 = 7850;          %   [#] densidad volumetrica del acero
m0 = dx*dy*r0;      %   [#] masa de una particula. (dz=1)
v_s = 300;          %   [#] velocidad incial
    
    ss = 4699;      %   [#]  ---
    C = 3630;       %   [#] Huggoniot
    S = 1800;       %   [#]  ---
    
    gamma = 1.81;   %   [#]  ---
    alpha = 0.5;    %   [#]
    beta = 0.5;     %   [#] XSPH
    eta = 0.01;     %   [#]
    eps = 0.5;      %   [#]  ---

G = 8*10^10;        %   [#] Modulo Cortante
Yo = 6*10^8;        %   [#] Esfuerzo de fluencia
E = ss^2*r0;        %   [#] Modulo de Young - calculado con parametros 
                    %        de Huggoniot
ro2 = 2750;         %   [#] Densidad aluminio ~ creo que es la densidad 
                    %        del proyectil
m = 3;          %   [#] parametros de Weibull para generacion de
    k = 7;          %   [#] fractuas al interior del material
    
V = dx*dy;          %   [#] Volumen infinitesimal                    

%% Asginacion de fallas
%

numpart = length(coorbar(:,1));    %   [#] Numero de particulas. ~size(1,coorbar)
Nflaws = numpart*log(numpart);%[#] Numero de fallas en el material
Nflaws = round(Nflaws);         % ~Modelo de fallas
Flawpart = randi(numpart,Nflaws,1,'uint32'); %[Nflaws x 1]
                            %   vector que contiene # aleatorio entre 1 y
                            %   numpart para cada una de las Nflaws
[Flaws{1:numpart,1}] = deal(zeros(0)); % Flaws [cell array] {numrpart x 1}
                            %   Cell array con 'numrpart' matrices vacias
                            %   Cell de prealocaci�n para fallas


for i = 1:length(Flawpart)                            
    Flaws{Flawpart(i),1}(size(Flaws{Flawpart(i),1})+1) = (i/(k*V))^(1/m);
end
% Flawpart  [Nflaws x 1]     vector que contiene Nflaws posiciones de las 
% fallas a generar. Las posiciones estan dadas como numeros enteros en 
% relacion al numero de cada particula.
% Flaws     cell{numrpart x 1}  contiene 'numrpart' matrices vacias que 
% identifican las fallas para cada particula.
% Flaws se va llenando con las posciones que indica Flawpart. La primera
% vez que pasa por una matriz de Flaws, le imprime el numero 
% entero (i/(k*v))^(1/m) formando una matriz 1x1. La segunda vez que pasa 
% por la misma matriz, aumenta la dimension de la matriz en una sola
% direccion para imprimir otro numero entero. Asi, si semialeatoriamente,
% el numero 'k' aparecio 'n' veces en Flawpart, la celda Flaws en su
% posicion 'k' debe contener una matriz nx1 de numero enteros.

%% Matrices de propiedades fisicas
V1 = zeros(numpart,1);          %[numpart x 1] Velocidad 1 -> 0's
V2 = ones(numpart,1);           %[numpart x 1] Velocidad 2 -> 1�s

M = (ones(numpart,1))*dx*dy*r0; %[numpart x 1] Masa de las particulas
Rho = (ones(numpart,1))*r0;     %[numpart x 1] Densidad de las particulas
drho = zeros(numpart,1);        %[numpart x 1] Derivada de la densidad
D = zeros(numpart,1);           %[numpart x 1] Da�o de cada part�cula [0 1]
cs = ss*ones(numpart,1);        %[numpart x 1] ~Velocidad del sonido

% Aceleracion
dV1 = zeros(numpart,1);     %[numpart x 1] Derivada total de V1
dV2 = zeros(numpart,1);     %[numpart x 1] Derivada total de V2
    % Derivadas Espaciales
    dv1dx1 = zeros(numpart,1);  %[numpart x 1] Derivada V1 en dir 1 (x)
    dv1dx2 = zeros(numpart,1);  %[numpart x 1] Derivada V1 en dir 2 (y)
    dv2dx1 = zeros(numpart,1);  %[numpart x 1] Derivada V2 en dir 1 (x)
    dv2dx2 = zeros(numpart,1);  %[numpart x 1] Derivada V2 en dir 2 (y)

P = zeros(numpart,1);       %[numpart x 1] Presion Hidrostatica

% Esfuerzos
dev11 = zeros(numpart,1);   % -----------
dev12 = zeros(numpart,1);   % Esfuerzos Cortantes
dev21 = zeros(numpart,1);   % [numpart x1]
dev22 = zeros(numpart,1);   % ------------
    % Derivadas Espaciales
    ddev11 = zeros(numpart,1);  %------------
    ddev12 = zeros(numpart,1);  % Derivadas Espaciales de los esfuerzos
    ddev21 = zeros(numpart,1);  % [numpart x 1]
    ddev22 = zeros(numpart,1);  %------------

% Deformaciones
eps11 = zeros(numpart,1);   %------------
eps12 = zeros(numpart,1);   % Deformacion Unitaria
eps21 = zeros(numpart,1);   %   [numpart x 1]
eps22 = zeros(numpart,1);   %------------

eint = zeros(numpart,1);    % Energ�a interna y derivada
deint = zeros(numpart,1);   %   [numpart x 1]

%% Distribucion de Velocidad
%Vdistr = linspace(-v_s,v_s,67); %[1 x 67] Distribucion de velocidad
                                % vector fila, espacio discretizado en 67
                                % puntos desde -v_s hasta v_s
%for i=1:(numpart/67) 
%    V2((i-1)*67+1:i*67) = Vdistr';
%end
V2(1:num_bar) = -10;


%% Inicio de Simulacion
dt = 1*10^-7;       %   [#] Paso de tiempo
tf = 1e-6*20;       %   [#] Tiempo final

t = 0;              %   int     Tiempo inicial
steps = tf/dt;      %   int     Numero de pasos

for ti=1:steps
    
    [Nearpart,Dist] = rangesearch(coorbar,coorbar,h);
    kern = Dist;
    kern = cellfun(@(x) x*0, kern, 'un',0);
    dkernx = kern;  %   cell {npart x 1} 
    dkerny = kern;  %   Derivadas del kern iniciadas en 0s
    
    for i = 1:num_bar
        kern{i} = kern1(Dist{i},h);
        dkernx{i} = dkernx1(Dist{i}, h, coorbar, Nearpart{i}, i);
        dkerny{i} = dkerny1(Dist{i}, h, coorbar, Nearpart{i}, i);        
        
        P(i) = EOSmie(eint(i), r0, C, S, Rho(i), gamma);
        
        [dv1dx1(i), dv1dx2(i), dv2dx1(i), dv2dx2(i)] = ...
            Velgradesp(Rho, M, V1, V2, dkernx{i,1}, dkerny{i,1},...
            Nearpart{i,1}, i);
        
        [eps11(i), eps12(i), eps21(i), eps22(i)] = ...
            deform(dv1dx1(i), dv1dx2(i), dv2dx1(i), dv2dx2(i));
        
        [ddev11(i), ddev12(i), ddev21(i), ddev22(i)] = ...
            devstresshooke(dv1dx2(i), dv2dx1(i),...
                dev11(i), dev12(i), dev21(i), dev22(i),...
                eps11(i), eps21(i), eps22(i), G);
        
        dev11(i) = dev11(i) + ddev11(i)*dt;
        dev12(i) = dev12(i) + ddev12(i)*dt;
        dev21(i) = dev21(i) + ddev21(i)*dt;
        dev22(i) = dev22(i) + ddev22(i)*dt;
        
        J = (dev11(i)^2) + 2*dev12(i)*dev21(i) + (dev22(i)^2);
       fact = sqrt(2*Yo/3);
       
       if J>Yo*3/2
           scalar = fact/sqrt(J);
           dev11(i) = dev11(i)*scalar;
           dev12(i) = dev12(i)*scalar;
           dev21(i) = dev21(i)*scalar;
           dev22(i) = dev22(i)*scalar;
       end
       
       drho(i) = derivadarho(M, V1, V2,...
            dkernx{i,1}, dkerny{i,1}, Nearpart{i,1}, i);
    end
    Rho(1:num_bar) = Rho(1:num_bar) + dt*drho(1:num_bar);
    Rho = real(Rho);
    
    for i=1:num_bar
        % Ecuacion de Momentum
        [dV1(i),dV2(i)] = Momentumeq2d(dev11, dev12, dev21, dev22,...
            P, Rho, dkernx{i,1}, dkerny{i,1}, M, Nearpart{i,1},...
            i, cs, Dist{i,1}, coorbar, h, V1, V2);
        
        % Conservacion de Energia
        deint(i) = Deint(dev11(i), dev12(i), dev21(i), dev22(i),...
            P, Rho, dkernx{i,1}, dkerny{i,1}, M, Nearpart{i,1},...
            i, cs, Dist{i,1}, coorbar, h, V1, V2,...
            eps11(i), eps12(i), eps21(i), eps22(i));
    end
    
    
    % Repulsion
    n1 = 12;
    n2 = 4;
    D_repulsion  = mean(sqrt(V1.^2+V2.^2));
    r0_repulsion = h*2;
      for i = num_bar+1:numpart
      [ddV1,ddV2] = Repulsion_frontera(coorbar, i, Nearpart, Dist,...
           D_repulsion, r0_repulsion, n1, n2, M);
       dV1 = dV1 + ddV1;
       dV2 = dV2 + ddV2;
      end
    
    V1(1:num_bar) = V1(1:num_bar) + dt*dV1(1:num_bar);
    V2(1:num_bar) = V2(1:num_bar) + dt*dV2(1:num_bar);
    eint(1:num_bar) = eint(1:num_bar) + dt*deint(1:num_bar);
    V1 = real(V1);
    V2 = real(V2);
    eint = real(eint);
    
    for i=1:num_bar
       %Velocidad del sonido
       cs(i) = Miespeedofsound(eint(i), r0, C, S, Rho(i), gamma);
       %XSPH
       [V1(i), V2(i)] = XSPH(Nearpart{i,1}, M, Rho, V1, V2, kern{i,1}, i);
       % Actualizar posiciones
       coorbar(i,:) = coorbar(i,:) + [V1(i), V2(i)]*dt;
    end
    
    
    %Graficas
    t = t + dt;
    figure(2)
    scatter(coorbar(:,1), coorbar(:,2),10,P,'filled')
    colorbar()
    axis([-0.01 0.01 -0.01 0.03])
    axis('equal')
    
    drawnow
end
%% Mostrar condicion inicial
scatter(coorbar(:,1),coorbar(:,2),10,V2,'filled')
colorbar()
axis('equal')