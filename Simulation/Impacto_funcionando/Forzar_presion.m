function [ P_nuevo ] = Forzar_presion( Presion, Vecinos, P_0, kernel, f )
% fuerzo la presion de la particula a ser como la presion de las particulas
% vecinas
n = length(Vecinos);
N = length(Presion);
p = f;
P_nuevo = P_0;
for i = 1:n
    indice = Vecinos(i);
    P_nuevo = P_nuevo + abs(p*kernel(i)*Presion(indice))/N;
end

end

