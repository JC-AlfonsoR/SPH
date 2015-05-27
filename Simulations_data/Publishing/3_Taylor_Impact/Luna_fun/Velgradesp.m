function [dv1dx1,dv1dx2,dv2dx1,dv2dx2 ] = Velgradesp(rho,m,V1,V2,dwx,dwy,Npart,Spart)

dv1dx1=0;
dv1dx2=0;
dv2dx1=0;
dv2dx2=0;

for i=1:length(Npart);
    fact=m(Npart(i))/rho(Npart(i));
    dv1=V1(Spart)-V1(Npart(i));
    dv2=V2(Spart)-V2(Npart(i));
    
    dx1=dwx(i);
    dx2=dwy(i);
    
    dv1dx1=dv1dx1+(dv1*dx1*fact);
    dv1dx2=dv1dx2+(dv1*dx2*fact);
    dv2dx1=dv2dx1+(dv2*dx1*fact);
    dv2dx2=dv2dx2+(dv2*dx2*fact);

end

