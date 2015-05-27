function [deint] = Deint(dev11,dev12,dev21,dev22,P,rho,dkernx,dkerny,m,Npart,selfpart,cs,Dist,coorglobal,h,V1,V2,eps11,eps12,eps21,eps22)
rhoa=rho(selfpart);
sumdeint=0;
V1s=V1(selfpart);
V2s=V2(selfpart);
for i=1:length(Npart)
    NP=Npart(i);
    rhob=rho(NP);
    mb=m(NP);
    dW=[dkernx dkerny];
    dx=abs(coorglobal(Npart(i),1)-coorglobal(selfpart,1));
    dy=abs(coorglobal(Npart(i),2)-coorglobal(selfpart,2));
    Mgh=Monaghanvisc(selfpart,NP,rho,cs,Dist(i),V1,V2,dx,dy,h);
    VIJ1=V1(NP)-V1s;
    VIJ2=V2(NP)-V2s;
    dvijdotdwij=VIJ1*dW(1)+VIJ2*dW(2);
    sumdeint=sumdeint+(mb*((P(selfpart)/(rhoa^2))+(P(NP)/(rhob^2))+Mgh)*dvijdotdwij);
    
end    

 devdoteps=(dev11*eps11)+(dev12*eps12)+(dev21*eps21)+(dev22*eps22);
 deint=1/2*sumdeint+(devdoteps/rhoa);

end
