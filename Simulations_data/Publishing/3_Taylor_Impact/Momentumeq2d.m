function [dV1,dV2] = Momentumeq2d(dev11,dev12,dev21,dev22,P,rho,dkernx,dkerny,m,Npart,selfpart,cs,Dist,coorglobal,h,V1,V2)
rhoa=rho(selfpart);
rhoa21=1/(rhoa*rhoa);
sigma11a=P(selfpart)+dev11(selfpart);
sigma12a=dev12(selfpart);
sigma21a=dev21(selfpart);
sigma22a=P(selfpart)+dev22(selfpart);
dV1=0;
dV2=0;
for i=1:length(Npart)
    NP=Npart(i);
    rhob=rho(NP);
    mb=m(NP);
    sigma11b=P(NP)+dev11(NP);
    sigma12b=dev12(NP);
    sigma21b=dev21(NP);
    sigma22b=P(NP)+dev22(NP);
    rhob21=1/(rhob*rhob);
    dW=[dkernx dkerny];
    dx=abs(coorglobal(Npart(i),1)-coorglobal(selfpart,1));
    dy=abs(coorglobal(Npart(i),2)-coorglobal(selfpart,2));
    Mgh=Monaghanvisc(selfpart,NP,rho,cs,Dist(i),V1,V2,dx,dy,h);
    dV1=dV1+mb*(sigma11a*rhoa21+sigma11b*rhob21+Mgh)*dW(1)+...
        mb*(sigma12a*rhoa21+sigma12b*rhob21)*dW(2);
    dV2=dV2+mb*(sigma21a*rhoa21+sigma21b*rhob21)*dW(1)+...
        mb*(sigma22a*rhoa21+sigma22b*rhob21+Mgh)*dW(2);
end

