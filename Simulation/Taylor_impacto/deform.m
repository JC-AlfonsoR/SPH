function [eps11,eps12,eps21,eps22] = deform(dv1dx1,dv1dx2,dv2dx1,dv2dx2)
% deform
%
% Calcula el tensor de deformación unitaria para un partícula. 
%
% [eps11, eps12, eps21, eps22] = deform(dv1dx1, dv1dx2, dv2dx1, dv2dx2)
%
%   Inputs
%   dv1dx1      double      Derivada de V1 en dir 1
%   dv1dx2      double      Derivada de V1 en dir 2
%   dv2dx1      double      Derivada de V2 en dir 1
%   dv2dx2      double      Derivada de V2 en dir 2
%
%   Outputs
%   eps11       double      Deformación unitaria en dir 11
%   eps12       double      Deformación unitaria en dir 12
%   eps21       double      Deformación unitaria en dir 21
%   eps22       double      Deformación unitaria en dir 22
%
        eps11 = dv1dx1;
        eps12 =0.5*(dv1dx2 + dv2dx1);
        eps21 = eps12;
        eps22 = dv2dx2;
end
%% Referencias
%{
http://es.wikipedia.org/wiki/Tensor_deformaci%C3%B3n
%}