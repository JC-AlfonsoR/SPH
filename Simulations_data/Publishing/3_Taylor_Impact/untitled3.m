for i = 1:600
   scatter(Coordenadas(:,1,i), Coordenadas(:,2,i),10,Velocidad2(:,1,i),'filled')
    colorbar()
    axis([-0.01 0.01 -0.01 0.03])
    axis('equal')
    
    drawnow
   %axis('equal')
end