function dJ=Jfun(t,x,param)

b=param.b;
g=param.g;

dJ=zeros(3,3);

    S=x(1);
    I=x(2);

    dJ(1,1) = - b * I;
    dJ(1,2) = - b * S;
    dJ(1,3) =  0;
    dJ(2,1) =  b * I;
    dJ(2,2) =  b * S - g;
    dJ(2,3) =  0;
    dJ(3,1) =  0;
    dJ(3,2) =  g;
    dJ(3,3) =  0;

end
