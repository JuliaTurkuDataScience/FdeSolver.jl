clear
clc
alpha = 0.6; lambda = -10 ;
f_fun = @(t,y,lam) lam * y;
J_fun = @(t,y,lam) lam ;
param = lambda ;
t0 = 0 ; T = 5 ; y0 = 1 ; h = 2^(-5) ;


tic
[t,y] = fde_pi12_pcIM(alpha,f_fun,J_fun,t0,T,y0,h,param,2);
toc
tic
[t1,y1] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param);
toc
tic
[t2,y2] = fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param,2);
toc

%%
clear
clc
h = 2^(-5) ;

f_fun = @(t,y,al) 40320/gamma(9-al)*t.^(8-al) - 3*gamma(5+al/2)/gamma(5-al/2)*t.^(4-al/2)+9/4*gamma(al+1) + (3/2*t.^(al/2)-t.^4).^3 - y.^(3/2) ;
J_fun = @(t,y,al) -3/2.*y.^(1/2) ;
alpha = 0.5 ;
param = alpha ;
t0 = 0 ; T = 1 ;
y0 = 0 ;


tic
[t,y] = fde_pi12_pcIM(alpha,f_fun,J_fun,t0,T,y0,h,param);
toc
tic
[~,y1] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param);
toc
tic
[~,y2] = fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param);
toc


yExact = @(t) t.^8 -3 .* t .^(4 +alpha/2) +9/4 * t.^(alpha);
Y_exact = yExact(t);

norm((Y_exact)-(y),inf)
norm((Y_exact)-(y1),inf)
norm((Y_exact)-(y2),inf)

%%
% 
alpha = [0.8,0.7] ;
A = 1 ; B = 3 ;
param = [ A , B ] ;
f_fun = @(t,y,par) [ ...
par(1) - (par(2)+1)*y(1) + y(1)^2*y(2) ; ...
par(2)*y(1) - y(1)^2*y(2) ] ;
% J_fun = @(t,y,par) [ ...
% -(par(2)+1) + 2*y(1)*y(2) , y(1)^2 ; ...
% par(2) - 2*y(1)*y(2) , -y(1)^2 ] ;
t0 = 0 ; T = 700 ; h = 2^(-5) ;
y0 = [ 1.2 ; 2.8] ;
tic
[t1,y1] = fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param);
toc

%%

clear
clc

h = 2^(-8) ;
alpha = [0.5, 0.2, 0.6] ;
f_fun = @(t,y) [(((y(2)-0.5).*(y(3)-0.3)).^(1/6) + sqrt(t))/sqrt(pi) ; ...
gamma(2.2)*(y(1)-1) ; ...
gamma(2.8)/gamma(2.2)*(y(2)-0.5) ] ;
J_fun = @(t,y) [0 , (y(2)-0.5).^(-5/6).*(y(3)-0.3).^(1/6)/6/sqrt(pi) , ...
(y(2)-0.5).^(1/6).*(y(3)-0.3).^(-5/6)/6/sqrt(pi) ; ...
gamma(2.2) , 0 , 0 ; ...
0 , gamma(2.8)/gamma(2.2) , 0 ] ;
t0 = 0 ; T = 5;
y0 = [ 1 ; 0.500000001 ; 0.300000001 ] ;

tic
[t,y] = fde_pi12_pcIM(alpha,f_fun,J_fun,t0,T,y0,h,[]);
toc
tic
[t1,y1] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,[],[],100);
toc
tic
[t2,y2] = fde_pi12_pc(alpha,f_fun,t0,T,y0,h,[],2);
toc


yExact = @(t) [t + 1;
                  t.^(1.2) + .5;
                  t.^(1.8)+.3];
              
Y_exact = yExact(t);

norm(real(Y_exact)-real(y),'fro')
norm(real(Y_exact)-real(y1),'fro')
norm(real(Y_exact)-real(y2),'fro')

% plot(t, Y_exact, t1,real(y1),':', t,real(y),'--', t,real(y2), '-.')


%%

clear 
clc


alpha = [0.8,0.7] ;
A = 1 ; B = 3 ;
param = [ A , B ] ;
f_fun = @(t,y,par) [par(1) - (par(2)+1)*y(1) + y(1)^2*y(2) ; ...
par(2)*y(1) - y(1)^2*y(2) ] ;
J_fun = @(t,y,par) [-(par(2)+1) + 2*y(1)*y(2) , y(1)^2 ; ...
par(2) - 2*y(1)*y(2) , -y(1)^2 ] ;
t0 = 0 ; T = 100 ; h = 2^(-5) ;
y0 = [ 1.2 ; 2.8] ;


tic
[t,y] = fde_pi12_pcIM(alpha,f_fun,J_fun,t0,T,y0,h,param);
toc
tic
[t1,y1] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param,[],2);
toc
tic
[t2,y2] = fde_pi12_pc(alpha,f_fun,t0,T,y0,h,param,2);
toc

% T = 100 ; h = 2^(-14) ;
% [tex,Y_exact] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h,param,10^-14);

% plot(t3, y3, t1,real(y1),':', t,real(y),'--', t,real(y2), '-.')

norm(real(y1)-real(y),inf)
norm(real(y1)-real(y2),inf)
% 
% 
% norm(real(Y_exact(:,indx))-real(y),'fro')
% norm(real(Y_exact(:,indx))-real(y1),'fro')
% norm(real(Y_exact(:,indx))-real(y2),'fro')
% norm(real(Y_exact(:,indx))-real(y3),'fro')
% 
