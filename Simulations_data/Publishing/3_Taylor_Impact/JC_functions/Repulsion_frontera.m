function [ dV1,dV2 ] = Repulsion_frontera( Particles, selfpart, Near_Part, Dist, D,...
    r0, n1, n2, M)
%% Repulsion
% [ dV1,dV2 ] = Repulsion_frontera( Particles, selfpart, Near_Part, Dist, D,...
%   r0, n1, n2, M)
%
% Inputs
% Particles [N_part*2] Posiciones x,y de las particulas
% Selfpart              Identidad de la particula
% Near_part {n*1}       Identidad de las n particulas vecinas
% Dist      {n*1}       Distancia de las n particulas vecinas

n = length(Near_Part{selfpart});
N_part = length(Particles);
F1 = zeros(N_part,1);
F2 = zeros(N_part,1);

for i = 1:n
    particula = Near_Part{selfpart}(i);
    rij = Dist{selfpart}(i);
    xij = Particles(particula,1) - Particles(selfpart,1);
    yij = Particles(particula,2) - Particles(selfpart,2);
    
    if r0/rij <=1
        F1(particula) = D*((r0/rij)^n1 - (r0/rij)^n2) * xij/rij^2;
        F2(particula) = D*((r0/rij)^n1 - (r0/rij)^n2) * yij/rij^2;
    end
end

dV1 = F1.*M;
dV2 = F2.*M;

end

