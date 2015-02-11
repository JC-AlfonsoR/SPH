clc
clear all; close all
% SPH programa principal Taylor bar 
%%
%Separacion entre particulas
dy = 0.000384848/2;
dx=dy;
% dyg=0.00384848;
% dxg=dyg;
% radio del dominio de soporte
hdx = 2.0;
h = hdx * dx;
%Densidad inicial
r0 = 2900;
%Masa planar
m0 = dx * dy * r0;
%Velocidad inicial
v_s = 300;
%Constantes de Hugoniot
ss = 4699;
C = 3630;
S = 1800;
%Constantes XSPH
gamma = 1.81;
alpha = 0.5;
beta = 0.5;
eta = 0.01;
eps = 0.5;
%propiedades del material
bar_width=0.0076;
%Modulo cortante
G = 8*10^10;
%Esfuerzo de fluencia
Yo = 6*10^8;
%Modulo de young
E=ss*ss*r0;
ro2= 2750;
platestart = -2*bar_width;
plateend = 2*bar_width;
%weibull constant
m=3;
k=7;
V=dx*dy;
%dt
dt = 1*10^-7;
tf = 0.000001*20;

%%
%Barra
xbar = -bar_width/2:dx:bar_width/2+dx;
ybar = 0:dx:0.0254;
[X, Y]=meshgrid(xbar,ybar);
coorbar = [X(:) Y(:)];
numrpart=length(X(:));
% numpart= numipart+numrpart;
Nflaws= numrpart*log(numrpart);
Nflaws=round(Nflaws);
% Vector de asignacion de fallas cada falla (posicion del vector)
% le corresponde un numero de particula (valor en el vector) 
Flawpart= randi(numrpart,Nflaws,1,'uint32');
coorglobal=coorbar; 
numpart=numrpart;
%Prealocacion de la celula de fallas
[Flaws{1:numrpart,1}]=deal(zeros(0));
%Asignacion del valor de falla
for i=1:length(Flawpart)
Flaws{Flawpart(i),1}(size(Flaws{Flawpart(i),1})+1)=(i/(k*V))^(1/m);
end 

%%
M=(ones(numpart,1))*dx*dy*r0;
Rho=(ones(numpart,1))*r0;
drho=(ones(numpart,1))*0;
D=(zeros(numpart,1));
%Velocidad del sonido
cs=(ones(numpart,1))*ss;
%%
%Velocidades
V1=(ones(numpart,1))*0;
V2=(ones(numpart,1));
Vdistr=linspace(v_s,0,20);
%{
for i=1:length(Vdistr)
    Bstart=27;
    Bfinish=40;
    V1(((i-1)*133) + (47+(i-1)) : ((i-1)*133) + (87-(i-1))) = Vdistr(i);
    
end
%}
for i = 1:length(Vdistr)
    l = 40;
    i0 = 18;
    V1(i0+133*(i-1) +i: i0+133*(i-1)+l -i) = Vdistr(i);
    V1(4*i0+133*(i-1) +i: 4*i0+133*(i-1)+l -i ) = Vdistr(i);
end
%
Velplot=zeros(size(X));
for i=1:numpart
Velplot(i)=V1(i);
end
image([-100 100], [-100 100], Velplot, 'CDataMapping', 'scaled');
colorbar
%
%%
%Accel
dV1=(ones(numpart,1))*0;
dV2=(ones(numpart,1))*0;
%Derivadas espaciales
dv1dx1=(ones(numpart,1))*0;
dv1dx2=(ones(numpart,1))*0;
dv2dx1=(ones(numpart,1))*0;
dv2dx2=(ones(numpart,1))*0;
%Presion hidroestatica
P=(ones(numpart,1))*0;
%Esfuerzos cortantes
dev11=(ones(numpart,1))*0;
dev12=(ones(numpart,1))*0;
dev21=(ones(numpart,1))*0;
dev22=(ones(numpart,1))*0;
%derivada esfuerzos cortantes
ddev11=(ones(numpart,1))*0;
ddev12=(ones(numpart,1))*0;
ddev21=(ones(numpart,1))*0;
ddev22=(ones(numpart,1))*0;

%deformaciones
eps11=(ones(numpart,1))*0;
eps12=(ones(numpart,1))*0;
eps21=(ones(numpart,1))*0;
eps22=(ones(numpart,1))*0;

eint=(ones(numpart,1))*0;        
deint=(ones(numpart,1))*0;

%---------------Inicio de simulacion --------------------
%%
t=0;
steps=tf/dt;

% Matrices para guardar informaci?n de la simulaci?n.
Coordenadas = zeros([size(coorbar), steps]);
Velocidad1 = zeros([size(V1),steps]);
Velocidad2 = zeros([size(V2),steps]);
Presion = zeros([size(P),steps]);
Esfuerzos11 = zeros([size(dev11), steps]);
Esfuerzos12 = zeros([size(dev12), steps]);
Esfuerzos21 = zeros([size(dev21), steps]);
Esfuerzos22 = zeros([size(dev22), steps]);
Densidad = zeros([size(Rho), steps]);

for ti=1:steps
%%
    fprintf('P %d / %d\n',ti,steps)
     [Nearpart,Dist]=rangesearch(coorbar,coorbar,h);
     kern=Dist;
     kern=cellfun(@(x) x*0,kern,'un',0);
     dkernx=kern;
     dkerny=kern;
%%    
     for i=1:numrpart
         kern{i}=kern1(Dist{i},h);
         dkernx{i}=dkernx1(Dist{i},h,coorbar,Nearpart{i},i);
         dkerny{i}=dkerny1(Dist{i},h,coorbar,Nearpart{i},i);
%%   
         % Ecuacion de estado Mie Gruniensen
        P(i)=EOSmie(eint(i),r0,C,S,Rho(i),gamma);
%%       %Derivadas espaciales de la velocidad
        [dv1dx1(i),dv1dx2(i),dv2dx1(i),dv2dx2(i)]=Velgradesp(Rho,M,V1,V2,dkernx{i,1},dkerny{i,1},Nearpart{i,1},i);
%%       %Calcvulo de deformaciones unitarias   
        [eps11(i),eps12(i),eps21(i),eps22(i)]=deform(dv1dx1(i),dv1dx2(i),dv2dx1(i),dv2dx2(i));
%%          %Derivadas de los esfuerzos cortantes
        [ddev11(i),ddev12(i),ddev21(i),ddev22(i)]=devstresshooke(dv1dx2(i),dv2dx1(i),dev11(i),dev12(i),dev21(i),dev22(i)...
            ,eps11(i),eps21(i),eps22(i),G);
        
        dev11(i)=dev11(i)+ddev11(i)*dt;
        dev12(i)=dev12(i)+ddev12(i)*dt;
        dev21(i)=dev21(i)+ddev21(i)*dt;
        dev22(i)=dev22(i)+ddev22(i)*dt;
        
%%        %Criterio de falla: Von mises
        J=(dev11(i)^2)+(2*dev12(i)*dev21(i))+(dev22(i)^2);
        fact=sqrt(2*Yo/3);
        if J>Yo*3/2
            scalar=fact/sqrt(J);
            dev11(i)=dev11(i)*scalar;
            dev12(i)=dev12(i)*scalar;
            dev22(i)=dev22(i)*scalar;
        end
%%      Ecuacion de continuidad
     drho(i)=derivadarho(M,V1,V2,dkernx{i,1},dkerny{i,1},Nearpart{i,1},i);
     
     end
     
     %% Paso de tiempo
     for i=1:numpart;
     Rho(i)=Rho(i)+(dt*drho(i));
     end
%%
     for i=1:numpart
        % Ecuacion de momentum
        [dV1(i),dV2(i)]=Momentumeq2d(dev11,dev12,dev21,dev22,P,Rho,dkernx{i,1},dkerny{i,1},M,Nearpart{i,1},i,cs,Dist{i,1},coorbar,h,V1,V2);
        % Ecuacion de conservacion de energia  
        deint(i)=Deint(dev11(i),dev12(i),dev21(i),dev22(i),P,Rho,dkernx{i,1},dkerny{i,1},M,Nearpart{i,1},i,cs,Dist{i,1},coorbar,h,V1,V2,eps11(i),eps12(i),eps21(i),eps22(i));
     end
     %%
     V1(1:length(dV1),1)=V1(1:length(dV1),1)+dt*dV1;
     V2(1:length(dV2),1)=V2(1:length(dV2),1)+dt*dV2;
     eint(1:length(deint),1)=eint(1:length(deint),1)+deint*dt; 
     %% Correciones
     for i=1:numpart
         %Velocidad del sonido
         cs(i)= Miespeedofsound(eint(i),r0,C,S,Rho(i),gamma);
         %XSPH
         [V1(i),V2(i)]=XSPH(Nearpart{i,1},M,Rho,V1,V2,kern{i,1},i);
         %Actualizar posiciones
         coorbar(i,:)=coorbar(i,:)+[V1(i) V2(i)]*dt;
         
     end
     %%
     
    t=t+dt;
    
    plot(coorbar(:,1),coorbar(:,2),'.b');
    axis([-0.0254 0.0254 -0.0254*1 0.0254*1])
    t=t+dt;
     
    drawnow
    
    %%Guardar Info de la simulaci?n.
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
%}