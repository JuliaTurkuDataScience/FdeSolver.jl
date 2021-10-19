function dx=fun(t,x)

global b N Ki

dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(b(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end

function fi=fi_Xk(i, x)
global n N Kij
fi=1;
K=1:N;K(i)=[];
for j=1:N-1
    k=K(j);
fi=fi*(Kij(k).^n/(Kij(k).^n+x(k).^n));
end

end
