function dx=funSIR(t,x,param)

b=param.b;
g=param.g;
S=x(1);I=x(2);

dx=zeros(3,1);

    dx(1) = - b .* S .* I;
    dx(2) = b .* S .* I - g .* I;
    dx(3) = g .* I;
end
