clear 
clc

% number of species
N=20; 

Kij=4*rand(N,N);

% Order of derivatives,  0<mu(i)=<1    
mu=ones(1,N);   

% Hill coefficient
n=2; 

% death rate
Ki=1*rand(N,1); 

%  final time
T=5000;


%Growth rates
b=2*rand(N,1); 

% Initial conditions

X0=2*rand(N,1); 

H=[2^(-5), 2^(-6),  2^(-7), 2^(-8)];
t0=0; % initial time



param.n=n;
param.N=N;
param.Ki=Ki;
param.Kij=Kij;
param.b=b;
Fun=@fun1;

Bench=zeros(length(H));

for i=1:length(H)
h=H(i);

Bench(i) = timeit(@() fde_pi12_pc(mu,Fun,t0,T,X0,h,param));

end

figure
loglog(H,Bench, 'LineWidth',3);
ylabel("Execution time (Sc)")
xlabel("step size")
csvwrite('BenchLongTerm.csv',Bench)