function [eps11,eps12,eps21,eps22] = deform(dv1dx1,dv1dx2,dv2dx1,dv2dx2)

        eps11 = dv1dx1;
        eps12 =0.5*(dv1dx2 + dv2dx1);
        eps21 = eps12;
        eps22 = dv2dx2;

end

