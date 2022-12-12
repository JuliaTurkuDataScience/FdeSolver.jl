clear
clc

H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];% setting step-size

%equation
F = @(t,y,al) 40320/gamma(9-al)*t.^(8-al) -3*gamma(5+al/2)/gamma(5-al/2)*t.^(4-al/2)+9/4*gamma(al+1) +(3/2*t.^(al/2)-t.^4).^3 - y.^(3/2) ;
% jacobian of the equation
JF = @(t,y,al) -3/2.*y.^(1/2) ;
%inputs
alpha = 0.5 ; %order of derivative
param = alpha ; %parameter
t0 = 0 ; T = 1 ; % initial and final time
y0 = 0 ; % initial value

% benchmarking
Bench=zeros(length(H),4,2);
for i=1:length(H)
    h=H(i);
%computting the time
Bench(i,1,1) = timeit(@() fde_pi1_ex(alpha,F,t0,T,y0,h,param));
Bench(i,2,1) = timeit(@() fde_pi12_pc(alpha,F,t0,T,y0,h,param));
Bench(i,3,1) = timeit(@() fde_pi1_im(alpha,F,JF,t0,T,y0,h,param));
Bench(i,4,1) = timeit(@() fde_pi2_im(alpha,F,JF,t0,T,y0,h,param));
%computting the error
[t1,y1]=fde_pi1_ex(alpha,F,t0,T,y0,h,param);
[t2,y2]=fde_pi12_pc(alpha,F,t0,T,y0,h,param);
[t3,y3]=fde_pi1_im(alpha,F,JF,t0,T,y0,h,param);
[t4,y4]=fde_pi2_im(alpha,F,JF,t0,T,y0,h,param);

Exact=t1.^8 -3.*t1.^(4+alpha/2)+ 9/4*t1.^alpha;Ex=Exact(2:end);

% Bench(i,1,2)=norm((y1(2:end)-Ex)./Ex);% Let's try relative error
% Bench(i,2,2)=norm((y2(2:end)-Ex)./Ex);
% Bench(i,3,2)=norm((y3(2:end)-Ex)./Ex);
% Bench(i,4,2)=norm((y4(2:end)-Ex)./Ex);
Bench(i,1,2)=norm(y1-Exact);
Bench(i,2,2)=norm(y2-Exact);
Bench(i,3,2)=norm(y3-Exact);
Bench(i,4,2)=norm(y4-Exact);

end
%% plotting
figure; 
sgtitle('Benchmark for nonStiff')
subplot(1,2,1)
loglog(H,Bench(:,:,1), 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
legend("Explicit PI of rectanguar","Predictor-corrector PI rules","Implicit PI of rectanguar","Implicit PI trapezoidal rule")
subplot (1,2,2)
loglog(H,Bench(:,:,2),'LineWidth',3);
ylabel("Square norm of errors")
xlabel("step size")

figure
loglog(Bench(:,:,1),Bench(:,:,2), 'LineWidth',3);
xlabel("Execution time (Sc)")
ylabel("Square norm of errors")
legend("Explicit PI of rectanguar","Predictor-corrector PI rules","Implicit PI of rectanguar","Implicit PI trapezoidal rule")
%%
writematrix(Bench,'BenchNonStiff.csv')