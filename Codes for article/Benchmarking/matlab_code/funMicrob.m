function dx=funMicrob(t,x,param)

b=param.b;
N=param.N;
Ki=param.Ki;

dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(b(i).*fi_Xk1(i, x,param)-Ki(i).*x(i));
end
end

function fi=fi_Xk1(i, x, param)

n=param.n;
N=param.N;
Kij=param.Kij;

fi=1;
K=1:N;K(i)=[];
for j=1:N-1
    k=K(j);
fi=fi*(Kij(k).^n/(Kij(k).^n+x(k).^n));
end

end
