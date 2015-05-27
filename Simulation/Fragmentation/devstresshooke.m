function [ddev11,ddev12,ddev21,ddev22 ] = devstresshooke(dv1dx2,dv2dx1,dev11,dev12,dev21,dev22,eps11,eps21,eps22,mu)
%devstreesshooke
%
%
%
%[ddev11, ddev12, ddev21, ddev22] =
%           devstresshooke(dv1dx2, dv2dx1, dev11, dev12, dev21, dev22,...
%                           eps11, eps21, eps22, mu)
%
%   Inputs
%   dv1dx2      double      Derivada de V1 en dir 2
%   dv2dx1      double      Derivada de V2 en dir 1
%   dev11       double      Componente Esfuerzos en 11
%   dev12       double      Componente Esfuerzos en 12
%   dev21       double      Componente Esfuerzos en 21
%   dev22       double      Componente Esfuerzos en 22
%   eps11       double      Deformación unitaria en 11
%   eps21       double      Deformación unitaria en 21
%   eps22       double      Deformación unitaria en 22
%   mu          double      Modulo de Cortante
%
%   Outputs
%   ddev11      double      Derivada de Esfuerzos en 11
%   ddev12      double      Derivada de Esfuerzos en 12
%   ddev21      double      Derivada de Esfuerzos en 21
%   ddev22      double      Derivada de Esfuerzos en 22
%

 omega12 = 0.5 * (dv1dx2-dv2dx1);
 omega21 = -omega12;
 fact = 2*mu;
 traza = (1/3)*(eps11+eps22);

 ddev11 = fact*(eps11-traza) + (dev12*omega12) + (dev21*omega12);

 ddev22 = fact*(eps22-traza) + (dev21*omega21) + (dev12*omega12);
       
 ddev12 = fact*(eps21) + (dev11*omega21) + (dev22*omega21);
 ddev21 = 12;
        
end
%% Referencias
%{
[1] D.Luna & A. Gonzalez, Estudio computacional de la fragmentación de
materiales frágiles con el método de partículas suavizadas (SPH), Uniandes,
2015. (pg 9, Ecuaciones 12 - 14 );
%}

%% Comentarios Finales
%{
Preungtarle a Daniel como trabaja esta función
%}