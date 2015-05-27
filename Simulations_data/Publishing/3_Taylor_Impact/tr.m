%% Jc_tensile
%{
%{
Este script consiste en una copia de el escript tensileroadsph.m autoria de
ing. Daniel Luna. 
Desarrollo este script para entender la estructura general del c�digo de
referencia.
%}
clear all;clc;
%% Constantes
%{
Todas las unidades en sistema internacional de unidades
%}
dy = 0.000384848;   %   [1x1] separacion entre particulas
dx = dy;
hdx = 5.0;          %   [#] cte para expandir radio
h = hdx*dx;         %   [#] Radio para dominio soporte
r0 = 7850;          %   [#] densidad volumetrica del acero
m0 = dx*dy*r0;      %   [#] masa de una particula. (dz=1)
v_s = 1.0;          %   [#] velocidad incial
    
    ss = 4699;      %   [#]  ---
    C = 3630;       %   [#] Huggoniot
    S = 1800;       %   [#]  ---
    
    gamma = 1.81;   %   [#]  ---
    alpha = 0.5;    %   [#]
    beta = 0.5;     %   [#] XSPH
    eta = 0.01;     %   [#]
    eps = 0.5;      %   [#]  ---
    
bar_width = 0.0076; %   [#] Ancho de la barra
G = 10*10^10;        %   [#] Modulo Cortante
Yo = 6*10^8;        %   [#] Esfuerzo de fluencia
E = ss^2*r0;        %   [#] Modulo de Young - calculado con parametros 
                    %        de Huggoniot
ro2 = 2750;         %   [#] Densidad aluminio ~ creo que es la densidad 
                    %        del proyectil
platestart = -2*bar_width; %[#] Dominio 2D debe empezar en esta posici�n 
plateend = 2*bar_width;    %[#] Dominio 2D debe terminar en esta posici�n
    
    m = 3;          %   [#] parametros de Weibull para generacion de
    k = 7;          %   [#] fractuas al interior del material
    
V = dx*dy;          %   [#] Volumen infinitesimal
dt = 1*10^-7;       %   [#] Paso de tiempo
tf = 3e-6*20;       %   [#] Tiempo final

%% Definir Geoemtr�a
%{
En esta secci�n se crean las posiciones de las particulas y se asginan las
fallas.
%}
xbar = -bar_width/2:dx:bar_width/2+dx;  % [#:dx:#] posiciones en x
ybar = 0:dy:0.0254;         %   [#:dy:#]    posiciones en y
[X,Y] = meshgrid(xbar,ybar);%   [Matriz]    Matrices que con la malla 
                            %       de las posiciones para las particulas

coorbar = [X(:),Y(:)];      %   | x_1 , y_1 |   En esta matriz organiza    
                            %   | .   , .   |   todas las posiciones de las
                            %   | x_n , y_n |   particulas.
                            %   [npart x 2]

numrpart = length(X(:));    %   [#] Numero de particulas. ~size(1,coorbar)
Nflaws = numrpart*log(numrpart);%[#] Numero de fallas en el material
Nflaws = round(Nflaws);         % ~Modelo de fallas
Flawpart = randi(numrpart,Nflaws,1,'uint32'); %[Nflaws x 1]
                            %   vector que contiene # aleatorio entre 1 y
                            %   numpart para cada una de las Nflaws
coorglobal = coorbar;       %   [Matriz] copia de la matriz coorbar
numpart = numrpart;         %   [#] copia de numero de particulas

[Flaws{1:numrpart,1}] = deal(zeros(0)); % Flaws [cell array] {numrpart x 1}
                            %   Cell array con 'numrpart' matrices vacias
                            %   Cell de prealocaci�n para fallas


for i = 1:length(Flawpart)                            
    Flaws{Flawpart(i),1}(size(Flaws{Flawpart(i),1})+1) = (i/(k*V))^(1/m);
end
% Flawpart  [Nflaws x 1]     vector que contiene Nflaws posiciones de las 
% fallas a generar. Las posiciones estan dadas como numeros enteros en 
% relaci�n al n�mero de cada part�cula.
% Flaws     cell{numrpart x 1}  contiene 'numrpart' matrices vacias que 
% identifican las fallas para cada particula.
% Flaws se va llenando con las posciones que indica Flawpart. La primera
% vez que pasa por una matriz de Flaws, le imprime el n�mero 
% entero (i/(k*v))^(1/m) formando una matriz 1x1. La segunda vez que pasa 
% por la misma matriz, aumenta la dimension de la matriz en una sola
% direccion para imprimir otro n�mero entero. As�, si semialeatoriamente,
% el numero 'k' aparecio 'n' veces en Flawpart, la celda Flaws en su
% posicion 'k' debe contener una matriz nx1 de n�mero enteros.

%% Otras Matrices
% ~Luna no espececifica comentarios para estas matrices
M = (ones(numpart,1))*dx*dy*r0; %[numpart x 1] Masa de las particulas
Rho = (ones(numpart,1))*r0;     %[numpart x 1] Densidad de las particulas
drho = zeros(numpart,1);        %[numpart x 1] Derivada de la densidad
D = zeros(numpart,1);           %[numpart x 1] Da�o de cada part�cula [0 1]
cs = ss*ones(numpart,1);        %[numpart x 1] ~Velocidad del sonido

%% Velocidades
V1 = zeros(numpart,1);          %[numpart x 1] Velocidad 1 -> 0's
V2 = ones(numpart,1);           %[numpart x 1] Velocidad 2 -> 1�s

Vdistr = linspace(-v_s,v_s,67); %[1 x 67] Distribucion de velocidad
                                % vector fila, espacio discretizado en 67
                                % puntos desde -v_s hasta v_s
for i=1:(numrpart/67) 
    V2((i-1)*67+1:i*67) = Vdistr';
end
% El recorrido se realiza numpart/67 = 1407/67 = 21 veces. En cada uno de
% los recorridos, Vdistr se convierte a columna y se guarda en V2. As�,
% Vdistr se repite 21 veces dentro de V2. Creo que el #67 se debe a la
% forma en la que esta creada la malla, indicando que la distribuci�n de
% velocidad se realiza sobre 67 particulas en fila.

Velplot = zeros(size(X));       % [Matriz] Matriz para graficar. Tiene el 
                                %   mismo tama�o que la matriz de malla X.
for i=1:numpart
   Velplot(i) = V2(i);
end
% La matriz Velplot se llena con los datos del vector V2. De esta forma,
% Velplot queda con 21 columnas de 67 filas. Cada Columna conteniendo el
% vector Vdistr

%image([-100 100], [-100 100], Velplot, 'CDataMapping', 'scaled');
%xlabel('X (normalizado a 100)'); ylabel('Y (normalizado a 100)');colorbar;
% Grafica Velplot como una im�gen. En x se presentan 21 puntos y en y se
% presentan 67 puntos. La velocidad no cambia con x.

%% Derivadas
%{
Aca se inicializan las matrices de todas las derivadas
%}

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
%% Incio de la simulaci�n
%{
Aca se corre la simulaci�nn
%}
display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
display('   Inicio de la Simulaci�n   ');
display('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');

t = 0;              %   int     Tiempo inicial
steps = tf/dt;      %   int     Numero de pasos

% Matrices para guardar informaci�n de la simulaci�n.
Coordenadas = zeros([size(coorbar), steps]);
Velocidad1 = zeros([size(V1),steps]);
Velocidad2 = zeros([size(V2),steps]);
Presion = zeros([size(P),steps]);
%Esfuerzos11 = zeros([size(dev11), steps]);
%Esfuerzos12 = zeros([size(dev12), steps]);
%Esfuerzos21 = zeros([size(dev21), steps]);
%Esfuerzos22 = zeros([size(dev22), steps]);
Densidad = zeros([size(Rho), steps]);
%}
for ti=1:steps
%%
    fprintf('P %d/%d %t E_int %d %t E_mec %d\n',ti,steps,sum(eint), ...
        0.5*sum(M.*(V1^2+V2^2)));
    
    %~~~~~ Inicializaci�n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    for i = 1:numpart 
    
       
        %%
        % Deformacion unitaria
        
%   Para la part�cula i, calcula lasc componentes del tensor dedeformaci�n

        %%
        % Derivadas de Esfuerzos Cortantes
        
%  Para la part�cula i, calcula las componentes del tensor de esfuerzos              

        
%{
    Una vez calculada la derivada de los esfuerzos para la particula i, 
    los esfuerzos para la part�cula i se calculan con un incremento
    infinitesimal de tiempo.
%}
       %%
       % Von Mises
       
       % * preguntar po est parte ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*****
       %Calculo del da�o
        %ei = J/((1-D(i))*E);
        %{
        for j=1:length(Flaws{i});
            count = 0;
            if ei>length(Flaws{i}(j)) %Posible error, Flaws{i}(j) es un 
                %double, por lo que length(flaws{i}(j)) siempre es 1.
                count = count+1;
            end
        end
        
        dD = (Damageevol(Rho(i), M(i), cs(i)))*count;
        D(i) = (dD^2)*dt;
        %}  
        %dev11(i) = dev11(i);%*(1-D(i));
        %dev12(i) = dev12(i);%*(1-D(i));
        %dev22(i) = dev22(i);%*(1-D(i));
        %dev22(i) = dev22(i);%*(1-D(i));
        % Alivio de Esfuerzos generado por da�o
        %%
        % Ecuaci�n de Continuidad
        
    % Calcula la derivada de la densidad para la particula i
    end
    
   
   
    
    %% Correcciones
    
    %%
    t = t + dt;
    figure(2)
    scatter(coorbar(:,1), coorbar(:,2),10,V2,'filled')
    colorbar()
    axis([-0.01 0.01 -0.01 0.03])
    axis('equal')
    
    drawnow
    
    %Guardar Info de la simulaci�n.
    Coordenadas(:,:,ti) = coorbar;
    Velocidad1(:,:,ti) = V1;
    Velocidad2(:,:,ti) = V2;
    Presion(:,:,ti) = P;
    Esfuerzos11(:,:,ti) = dev11;
    Esfuerzos12(:,:,ti) = dev12;
    Esfuerzos21(:,:,ti) = dev21;
    Esfuerzos22(:,:,ti) = dev22;
    Densidad(:,:,ti) = Rho;
    %}  
    % Con estas matrices grabo las posiciones y propiedades de todas las
    % particulas para cada incremento de tiempo
end
%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentaci�n de
materiales fragiles con el m�todo de particulas suavizadas (SPH), Uniandes,
2015.

[2] G.R. Liu & M.B. Liu, Smoothed Particle Hydrodynamics - a meshfree particle
ethod, World Scientifics Publishing Co., 2003.
%}

%% Comentarios JC
%{
Hasta aca entiendo todo
Hay algunas modificaciones que quisiera hacer. Todas consisten en
modificaicones de estilo, por ejemplo: crear una matriz de 0's en lugar de
crear una matriz de 1's y multiplicarla por 0
%}
