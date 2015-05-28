%imapct_data_graphics
clear all; close all; clc

load('Impact_data.mat')

n = 30;
% Generar gifs
%{
Velocidad2 = abs(real(Velocidad2));
f_name = 'Velociadad_2.gif';
for i = 1:n
   
   scatter(Coordenadas_x(:,i),Coordenadas_y(:,i),10,Velocidad2(:,i),'filled'); 
   colorbar()
   caxis([0 1e3]);
   title(f_name)
   drawnow
   
   frame = getframe(1);
   im = frame2im(frame);
   [imind,cm] = rgb2ind(im,256);
   
   if i == 1
       imwrite(imind,cm,f_name,'gif', 'Loopcount',inf);
   else
       imwrite(imind,cm,f_name,'gif','WriteMode','append');
   end
   
end
%}

% Subplots
Velocidad2 = abs(real(Velocidad2));
c = 1;
figure(1)
for i = [15,19,21,23,25,30]
    %if mod(i-1,5) == 0
   subplot(2,3,c); c=c+1;
   scatter(Coordenadas_x(:,i),Coordenadas_y(:,i),5,Presion(:,i),'filled'); 
   caxis([0 1e7]);
   xlim([-4e-3,0.5e-3])
        ylim([-3e-3,3e-3])
   axis('equal')
   axis off
   if c==3
       title('Presion [Pa]')
   end
   
   %subplot(3,3,7:9)
   %scatter(Coordenadas_x(:,n),Coordenadas_y(:,n),5,Velocidad2(:,n),'filled'); 
   %caxis([0 1e3]);
   %colorbar('south')
   %caxis([0 1e3]);
   %axis off
   %subplot(3,3,9)
   %colorbar()
   %axis off
end
figure(2)
   caxis([0 1e7]);
   colorbar('south')
   axis off