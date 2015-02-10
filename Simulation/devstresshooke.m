function [ddev11,ddev12,ddev21,ddev22 ] = devstresshooke(dv1dx2,dv2dx1,dev11,dev12,dev21,dev22,eps11,eps21,eps22,mu)
 
 omega12 = 0.5 * (dv1dx2 - dv2dx1);
 omega21 = -omega12;
 fact = 2*mu;
 traza = (1/3)*(eps11 + eps22);

 ddev11 = fact*( eps11 - traza )+(dev12*omega12 ) + ( dev21*omega12 );

 ddev22 = fact*( eps22 - traza )+(dev21*omega21 ) + ( dev12*omega12 );
       
 ddev12 = fact*(eps21) + ( dev11*omega21 ) + ( dev22*omega21);
 ddev21=12;
        
end

