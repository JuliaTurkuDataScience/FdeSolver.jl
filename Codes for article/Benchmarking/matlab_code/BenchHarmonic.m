clear
clc

% setting step-size
H=[2^(-2) 2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];

%Equation
F= @(t,y,param) -param(1)/param(2) *y;
%Jacobian of the equation
JF = @(t,y,param) -param(1)/param(2);
%inputs
k=16;m=4;
alpha = 2; param=[k,m];
t0 = 0 ; T = 10; y0=[1,1];

%benchmarking
Bench=zeros(length(H),2,2);
for i=1:length(H)
    h=H(i);
%computting the runtime
Bench(i,1,1) = timeit(@() fde_pi1_ex(alpha,F,t0,T,y0,h,param));
Bench(i,2,1) = timeit(@() fde_pi12_pc(alpha,F,t0,T,y0,h,param));
Bench(i,3,1) = timeit(@() fde_pi1_im(alpha,F,JF,t0,T,y0,h,param));
Bench(i,4,1) = timeit(@() fde_pi2_im(alpha,F,JF,t0,T,y0,h,param));
%computing the error
[t1,y1]=fde_pi1_ex(alpha,F,t0,T,y0,h,param);
[t2,y2]=fde_pi12_pc(alpha,F,t0,T,y0,h,param);
[t3,y3]=fde_pi1_im(alpha,F,JF,t0,T,y0,h,param);
[t4,y4]=fde_pi2_im(alpha,F,JF,t0,T,y0,h,param);
%exact solution
Exact=y0(1).*cos(sqrt(k/m).*t1)+y0(2)./((sqrt(k/m))).*sin(sqrt(k/m).*t1);

Bench(i,1,2)=norm(y1-Exact);
Bench(i,2,2)=norm(y2-Exact);
Bench(i,3,2)=norm(y3-Exact);
Bench(i,4,2)=norm(y4-Exact);

end
%% plot
figure; 
sgtitle('Benchmark for Harmonic')
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
writematrix(Bench,'BenchHarmonic.csv')