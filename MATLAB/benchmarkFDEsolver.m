clear
clc

H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];
Bench1=zeros(length(H),2,2);

alpha = 0.6; lambda = -10 ;
f_fun = @(t,y,lam) lam * y;
J_fun = @(t,y,lam) lam ;
param = lambda ;
t0 = 0 ; T = 5 ; y0 = 1 ;

for i=1:length(H)
    h=H(i);
Bench1(i,1,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param));
Bench1(i,2,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param));

if i>=4
[t,y]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param);
[t1,y1]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param);

Exact=mlf(alpha,1,lambda*t.^alpha);

Bench1(i,1,2)=norm(y-Exact);
Bench1(i,2,2)=norm(y1-Exact);
end
end
%%
figure; sgtitle('Benchmark for example 1')
subplot(1,2,1)
loglog(H,Bench1(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Predictor-corrector PI rules","Implicit PI trapezoidal rule")
subplot (1,2,2)
loglog(H(4:end),Bench1(4:end,:,2),'LineWidth',3);
csvwrite('Bench1.csv',Bench1)
ylabel("Square norm of errors")
xlabel("step size")
%%
clear
clc

H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];
Bench2=zeros(length(H),2,2);

f_fun = @(t,y,al) 40320/gamma(9-al)*t.^(8-al) -3*gamma(5+al/2)/gamma(5-al/2)*t.^(4-al/2)+9/4*gamma(al+1) +(3/2*t.^(al/2)-t.^4).^3 - y.^(3/2) ;
J_fun = @(t,y,al) -3/2.*y.^(1/2) ;
alpha = 0.5 ;
param = alpha ;
t0 = 0 ; T = 1 ;
y0 = 0 ;

for i=1:length(H)
    h=H(i);
Bench2(i,1,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param));
Bench2(i,2,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param));

% if i>=4
[t,y]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param);
[t1,y1]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param);

Exact=t.^8 -3.*t.^(4+alpha/2)+ 9/4*t.^alpha;

Bench2(i,1,2)=norm(y-Exact);
Bench2(i,2,2)=norm(y1-Exact);
% end
end
%%
figure; sgtitle('Benchmark for example 2')
subplot(1,2,1)
loglog(H,Bench2(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Predictor-corrector PI rules","Implicit PI trapezoidal rule")
subplot (1,2,2)
loglog(H,Bench2(:,:,2),'LineWidth',3);
csvwrite('Bench2.csv',Bench2)
ylabel("Square norm of errors")
xlabel("step size")

%%

clear
clc

H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];
Bench3=zeros(length(H),2,2);
alpha = 0.7 ;
t0 = 0 ; T = 2 ;
y0 = 0 ;
param = alpha ;
J_p_fun=@(t,y,param) 0;
P_fun=@p_fun;

for i=1:length(H)
    h=H(i);
Bench3(i,1,1) = timeit(@() fde_pi12_pc(alpha,P_fun,t0,T,y0,h,param));
Bench3(i,2,1) = timeit(@() fde_pi2_im(alpha,P_fun,J_p_fun,t0,T,y0,h,param));

% if i>=4
[t,y]=fde_pi12_pc(alpha,P_fun,t0,T,y0,h,param);
[t1,y1]=fde_pi2_im(alpha,P_fun,J_p_fun,t0,T,y0,h,param);

syms x
 ExcatPiecW = piecewise(x<= 1, x ,x > 1, x-(x-1)^2);
Exact=subs( ExcatPiecW,x,t);

Bench3(i,1,2)=norm(y-Exact);
Bench3(i,2,2)=norm(y1-Exact);
% end
end
%%
figure; sgtitle('Benchmark for example 3')
subplot(1,2,1)
loglog(H,Bench3(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Predictor-corrector PI rules","Implicit PI trapezoidal rule")
subplot (1,2,2)
P=loglog(H,Bench3(:,1,2), H,Bench3(:,2,2),'--');
set(P,'LineWidth',3)
csvwrite('Bench3.csv',Bench3)
ylabel("Square norm of errors")
xlabel("step size")
   
%%
clear 
clc

H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];
Bench4=zeros(length(H),2,2);

alpha = [0.5, 0.2, 0.6] ;
f_fun = @(t,y) [(((y(2)-0.5).*(y(3)-0.3)).^(1/6) + sqrt(t))/sqrt(pi) ; ...
gamma(2.2)*(y(1)-1) ; ...
gamma(2.8)/gamma(2.2)*(y(2)-0.5) ] ;
J_fun = @(t,y) [0 , (y(2)-0.5).^(-5/6).*(y(3)-0.3).^(1/6)/6/sqrt(pi) , ...
(y(2)-0.5).^(1/6).*(y(3)-0.3).^(-5/6)/6/sqrt(pi) ; ...
gamma(2.2) , 0 , 0 ; 0 , gamma(2.8)/gamma(2.2) , 0 ] ;
t0 = 0 ; T = 5 ;
y0 = [ 1 ; 0.500000001 ; 0.300000001 ] ;


for i=1:length(H)
    h=H(i);
Bench4(i,1,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h));
Bench4(i,2,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h));

% if i>=4
[t,y]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h);
[t1,y1]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h);

Exact=[t+1; t.^(1.2)+0.5; t.^(1.8)+0.3];

Bench4(i,1,2)=norm(y-Exact);
Bench4(i,2,2)=norm(y1-Exact);
% end
end
%%
figure; sgtitle('Benchmark for example 4')
subplot(1,2,1)
loglog(H,Bench4(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Predictor-corrector PI rules","Implicit PI trapezoidal rule")
subplot (1,2,2)
loglog(H,Bench4(:,:,2),'LineWidth',3);
csvwrite('Bench4.csv',Bench4)
ylabel("Square norm of errors")
xlabel("step size")


%%
clear
clc

H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];
Bench5=zeros(length(H),2,2);

alpha = 1.5;param=alpha;
f_fun = @(t,y,al) t.^(al)*y+ 4*sqrt(t/pi)-t.^(2+al);
J_fun = @(t,y,al) t.^(al);
t0 = 0 ; T = 1; y0 = [0 0];

for i=1:length(H)
    h=H(i);
Bench5(i,1,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param));
Bench5(i,2,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param));


[t,y]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param);
[t1,y1]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param);

Exact=t.^2;

Bench5(i,1,2)=norm(y-Exact);
Bench5(i,2,2)=norm(y1-Exact);

end
%%
figure; sgtitle('Benchmark for example 5')
subplot(1,2,1)
loglog(H,Bench5(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Predictor-corrector PI rules","Implicit PI trapezoidal rule")
subplot (1,2,2)
loglog(H,Bench5(:,:,2),'LineWidth',3);
csvwrite('Bench5.csv',Bench5)
ylabel("Square norm of errors")
xlabel("step size")

%%
clear
clc

H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];
Bench6=zeros(length(H),2,2);

k=16;m=4;
alpha = 2; param=[k,m];
f_fun = @(t,y,param) -param(1)/param(2) *y;
J_fun = @(t,y,param) -param(1)/param(2);
t0 = 0 ; T = 10; y0=[1,1];

for i=1:length(H)
    h=H(i);
Bench6(i,1,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param));
Bench6(i,2,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param));


[t,y]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param);
[t1,y1]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param);

Exact=y0(1).*cos(sqrt(k/m).*t)+y0(2)./((sqrt(k/m))).*sin(sqrt(k/m).*t);
Bench6(i,1,2)=norm(y-Exact);
Bench6(i,2,2)=norm(y1-Exact);

end
%%
figure; sgtitle('Benchmark for example 6')
subplot(1,2,1)
loglog(H,Bench6(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Predictor-corrector PI rules","Implicit PI trapezoidal rule")
subplot (1,2,2)
loglog(H,Bench6(:,:,2),'LineWidth',3);
csvwrite('Bench6.csv',Bench6)
ylabel("Square norm of errors")
xlabel("step size")
