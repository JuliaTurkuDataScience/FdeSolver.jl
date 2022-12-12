%jacobian of Three-Species Lotka-Volterra Model 
%https://www.mdpi.com/2073-8994/13/3/368
function dJ=JfunLV(t,x,param)

a1=param.a1;
a2=param.a2;
a3=param.a3;
a4=param.a4;
a5=param.a5;
a6=param.a6;
a7=param.a7;

dJ=zeros(3,3);

    dJ(1,1) = a1-2*a2*x(1)-x(2)-x(3);
    dJ(1,2) = -x(1);
    dJ(1,3) = -x(1);
    dJ(2,1) =  a4*x(2);
    dJ(2,2) =  1-a3+a4*x(1);
    dJ(2,3) =  0;
    dJ(3,1) =  a6*x(3);
    dJ(3,2) =  a7*x(3);
    dJ(3,3) =  a6*x(1)-a5+a7*x(2)+1;

end
