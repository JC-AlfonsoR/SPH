function [kern] = kern1(rij,hk)
% Calcula el valor de Kernel par las partículas que se encuentran a 
%distancias rij. El calculo se hace con un dominio de radio hk
%
%[kern] = kern1(rij, hk)
%    
%    Inputs
%    rij     [1 x n]     Vector con las n distancias de las n particulas
%                        vecinas a la particula i
%    hk      double      Radio del dominio soporte para partícula i
%
%    Outputs
%    kern    [n x 1]     Vector con la magnitud del kernel evaluado para las
%                        n partículas vecinas
%

kern = zeros(length(rij),1);        % [n x 1]   vector de 0's

for i=1:length(rij)
    q = rij(i)/hk;  %q  double     R en ecuaciones de ref
    if q<=2
        fac = (7/((4*pi)*((hk)^2)));    
        val = ((1-(0.5*q))^4)*((2*q)+1);
        % Kernel de orden 5
        kern(i) = val*fac;
    else
        kern(i)=0;
    end
end

%% Comentarios Finales
%{
El kernel que de la función no es el Kernel de Johnson et al.,1996b
El de Johnson es de orden 2, mientras que este kernel es de orden 5.

%La forma del kernel fue verificada con el sigueinte codigo:
a = (0:.1:5); h = 2;
k = kern1(a,h);
plot(a,k); xlabel('Distancia R'); ylabel('Valor Kernel')

Se obtiene la forma general del kernel, pero no las propiedades
especificadas en referencias. A medida que h tiende a cero, se obtiene una
gráfica que en R=0 es mayor a 1. Para esta condición, el kernel no se
parece al delta de dirac a medida que h tiende a 0
~~ Verificar que la integral del kernel =~ 1

~~~> preguntarle a Daniel que kernel empleó
%}