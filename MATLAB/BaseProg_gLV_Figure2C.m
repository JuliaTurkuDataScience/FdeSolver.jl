clear
clc
global n N b Ki Kij

n=2; % Hill coefficient
N=3; % number of species
Kij=0.1*ones(N); % interaction matrix
Ki=1*ones(N,1); % death rate
t0=0; T=100; % t0 is start time, and T is the end time
b=[1, .95, 1.05]; % growth rate
x0=[.3,.1,.2]; % initial conditions
% alpha=0.5*ones(N,1); % order of derivatives
alpha=[1,1,1]; % order of derivatives

h=0.01; % step size for computing
tic
[t, x] = fde_pi12_pc(alpha,@fun,t0,T,x0',h);
toc

%%
% semilogx(t,x)
figure
p2e=plot(t,x(1,:),'b',t,x(2,:),'r',t,x(3,:),'g');set(p2e,'LineWidth',2)
