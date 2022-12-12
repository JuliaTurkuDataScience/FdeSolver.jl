%Equation of Three-Species Lotka-Volterra Model 
%https://www.mdpi.com/2073-8994/13/3/368

function dx=funLV(t,x,param)

a1=param.a1;
a2=param.a2;
a3=param.a3;
a4=param.a4;
a5=param.a5;
a6=param.a6;
a7=param.a7;

dx=zeros(3,1);

    dx(1) = x(1)*(a1-a2*x(2)-x(3));
    dx(2) = x(2)*(1-a3+a4*x(1));
    dx(3) = x(3)*(1-a5+a6*x(1)+a7*x(2));
end
