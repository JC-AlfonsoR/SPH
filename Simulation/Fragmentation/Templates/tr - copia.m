%% Jc_tensile
%{
Este script consiste en una copia de el escript tensileroadsph.m autoria de
ing. Daniel Luna. 
Desarrollo este script para entender la estructura general del c�digo de
referencia.
%}
%% Constantes
%{
Todas las unidades en sistema internacional de unidades
%}


v_s = 300;          %   [#] velocidad incial
    

ro2 = 2750;         %   [#] Densidad aluminio ~ creo que es la densidad 
                    %        del proyectil


%% Velocidades

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


%% Simulacion

t = 0;              %   int     Tiempo inicial
steps = tf/dt;      %   int     Numero de pasos


for ti=1:steps
    
    for i = 1:numpart 
        %%
        % Ecuaci�n de estado Mie Gruniensen
        P(i) = EOSmie(eint(i), r0, C, S, Rho(i), gamma);
%{
    Para la part�cula i, se calcula la presi�n hidrost�tica en funci�n de
    las demas variables de Estado. El calculo se hace usando la Ecuaci�n de
    Mie-Grniensen. Las variables de estado consisten en enrg�a interna,
    densidad y Parametros de la curva de Huggoniot.
%}
        %%
        % Derivada Espacial de Velocidad
        [dv1dx1(i), dv1dx2(i), dv2dx1(i), dv2dx2(i)] = ...
            Velgradesp(Rho, M, V1, V2, dkernx{i,1}, dkerny{i,1},...
            Nearpart{i,1}, i);
%{
    Para la part�cula i, se calculan las derivadas espaciales en
    direcciones 1 y 2 para las componentes de velocidad V1 y V2. El c�lculo
    se hace con aproximaci�n por part�culas vecinas.
%}
        %%
        % Deformacion unitaria
        [eps11(i), eps12(i), eps21(i), eps22(i)] = ...
            deform(dv1dx1(i), dv1dx2(i), dv2dx1(i), dv2dx2(i));
%   Para la part�cula i, calcula lasc componentes del tensor dedeformaci�n

        %%
        % Derivadas de Esfuerzos Cortantes
        [ddev11(i), ddev12(i), ddev21(i), ddev22(i)] = ...
            devstresshooke(dv1dx2(i), dv2dx1(i),...
                dev11(i), dev12(i), dev21(i), dev22(i),...
                eps11(i), eps21(i), eps22(i), G);
%  Para la part�cula i, calcula las componentes del tensor de esfuerzos              

        dev11(i) = dev11(i) + ddev11(i)*dt;
        dev12(i) = dev12(i) + ddev12(i)*dt;
        dev21(i) = dev21(i) + ddev21(i)*dt;
        dev22(i) = dev22(i) + ddev22(i)*dt;
%{
    Una vez calculada la derivada de los esfuerzos para la particula i, 
    los esfuerzos para la part�cula i se calculan con un incremento
    infinitesimal de tiempo.
%}
       %%
       % Von Mises
       J = (dev11(i)^2) + 2*dev12(i)*dev21(i) + (dev22(i)^2);
       fact = sqrt(2*Yo/3);
       
       if J>Yo*3/2
           scalar = fact/sqrt(J);
           dev11(i) = dev11(i)*scalar;
           dev12(i) = dev12(i)*scalar;
           dev21(i) = dev21(i)*scalar;
           dev22(i) = dev22(i)*scalar;
       end
       % * preguntar po est parte ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*****
       %Calculo del da�o
        ei = J/((1-D(i))*E);
        
        for j=1:length(Flaws{i});
            count = 0;
            if ei>length(Flaws{i}(j)) %Posible error, Flaws{i}(j) es un 
                %double, por lo que length(flaws{i}(j)) siempre es 1.
                count = count+1;
            end
        end
        
        dD = (Damageevol(Rho(i), M(i), cs(i)))*count;
        D(i) = (dD^2)*dt;
            
        dev11(i) = dev11(i)*(1-D(i));
        dev12(i) = dev12(i)*(1-D(i));
        dev22(i) = dev22(i)*(1-D(i));
        dev22(i) = dev21(i)*(1-D(i));
        % Alivio de Esfuerzos generado por da�o
        %%
        % Ecuaci�n de Continuidad
        drho(i) = derivadarho(M, V1, V2,...
            dkernx{i,1}, dkerny{i,1}, Nearpart{i,1}, i);
    % Calcula la derivada de la densidad para la particula i
    end
    
    % Paso de tiempo
    for i=1:numpart;
        Rho(i) = Rho(i) + dt*drho(i); % Cambia la densidad en el tiempo
    end
    
    for i=1:numpart
        % Ecuacion de Momentum
        [dV1(i),dV2(i)] = Momentumeq2d(dev11, dev12, dev21, dev22,...
            P, Rho, dkernx{i,1}, dkerny{i,1}, M, Nearpart{i,1},...
            i, cs, Dist{i,1}, coorglobal, h, V1, V2);
        
        % Conservacion de Energia
        deint(i) = Deint(dev11(i), dev12(i), dev21(i), dev22(i),...
            P, Rho, dkernx{i,1}, dkerny{i,1}, M, Nearpart{i,1},...
            i, cs, Dist{i,1}, coorglobal, h, V1, V2,...
            eps11(i), eps12(i), eps21(i), eps22(i));
    end
    %%
    V1(1:length(dV1),1) = V1(1:length(dV1),1) + dt*dV1;
    V2(1:length(dV2),1) = V2(1:length(dV2),1) + dt*dV2;
    eint(1:length(deint),1)=eint(1:length(deint),1) + deint*dt; 
    % Cambia la energia interna y las velocidades acorde al paso de tiempo
    % y las derivadas calculadas
    
    %% Correcciones
    for i=1:numpart
       %Velocidad del sonido
       cs(i) = Miespeedofsound(eint(i), r0, C, S, Rho(i), gamma);
       %XSPH
       [V1(i), V2(i)] = XSPH(Nearpart{i,1}, M, Rho, V1, V2, kern{i,1}, i);
       % Actualizar posiciones
       coorbar(i,:) = coorbar(i,:) + [V1(i), V2(i)]*dt;
    end
    %%
    t = t + dt;
    figure(2)
    plot(coorbar(:,1), coorbar(:,2),'.b')
    axis([-0.0254 0.0254 -0.0254*2 0.0254*2])
    drawnow
    
    %%Guardar Info de la simulaci�n.
    Coordenadas(:,:,ti) = coorbar;
    Velocidad1(:,:,ti) = V1;
    Velocidad2(:,:,ti) = V2;
    Presion(:,:,ti) = P;
    Esfuerzos11(:,:,ti) = dev11;
    Esfuerzos12(:,:,ti) = dev12;
    Esfuerzos21(:,:,ti) = dev21;
    Esfuerzos22(:,:,ti) = dev22;
    Densidad(:,:,ti) = Rho;
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
%}