%% bench 
    
clear
clc
%setting step-size
H=[2^(-4) 2^(-5) 2^(-6) 2^(-7) 2^(-8)];
%inputs
t0 = 0 ; T = 60; y0=[1;1;1];
alpha=[1,.9,.7];
param.a1=3;
param.a2=3;
param.a3=3;
param.a4=5;
param.a5=3;
param.a6=3;
param.a7=3;
%equation
F=@funLV;
%jacobian of Equation
JF=@JfunLV;
%fine solution
[tex,Exact]=fde_pi2_im(alpha,F,JF,t0,T,y0,2^(-10),param, 1e-12);


Bench=zeros(length(H),4,2);
for i=1:length(H)
    h=H(i);
%computing the runtime
Bench(i,1,1) = timeit(@() fde_pi1_ex(alpha,F,t0,T,y0,h,param));
Bench(i,2,1) = timeit(@() fde_pi12_pc(alpha,F,t0,T,y0,h,param,4,1e-8));
Bench(i,3,1) = timeit(@() fde_pi1_im(alpha,F,JF,t0,T,y0,h,param,1e-8));
Bench(i,4,1) = timeit(@() fde_pi2_im(alpha,F,JF,t0,T,y0,h,param,1e-8));
%computing the error
[t1,y1]=fde_pi1_ex(alpha,F,t0,T,y0,h,param);
[t2,y2]=fde_pi12_pc(alpha,F,t0,T,y0,h,param,4,1e-8);
[t3,y3]=fde_pi1_im(alpha,F,JF,t0,T,y0,h,param,1e-8);
[t4,y4]=fde_pi2_im(alpha,F,JF,t0,T,y0,h,param,1e-8);

Bench(i,1,2)=norm(y1-Exact(:,1:2^(7-i):end));
Bench(i,2,2)=norm(y2-Exact(:,1:2^(7-i):end));
Bench(i,3,2)=norm(y3-Exact(:,1:2^(7-i):end));
Bench(i,4,2)=norm(y4-Exact(:,1:2^(7-i):end));

end

%% plot
figure; 
sgtitle('Benchmark for SIR')
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
ylabel("Execution time (Sc)")
ylabel("Square norm of errors")
legend("Explicit PI of rectanguar","Predictor-corrector PI rules","Implicit PI of rectanguar","Implicit PI trapezoidal rule")
%% save the data
writematrix(Bench,'BenchLV.csv')
writematrix(Exact,'M_Exact_LV3.csv')