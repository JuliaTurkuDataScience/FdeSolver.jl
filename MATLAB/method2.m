function [t,x,B]=method2(Interaction,Perturbation,InitialCondition)
global T mu X0 b bb N n Ki

switch Interaction
    case 'Predefined'
        load('Kij.mat');
        %------------------------------------------------------------------
    case 'Random'
        load('KijRand.mat');
end

%--------------------------------------------------------------------------
    
switch Perturbation
    case 'False'
        Fun=@fun;
        %------------------------------------------------------------------
    case    'False&b*'
        load('bMean.mat');
        bb=b;
        Fun=@fun;
        %------------------------------------------------------------------
    case 'Pulse'
        load('bRand.mat');
        Fun=@fun2;
        % Check if the final time passes the pulse
        if T<330
           warning('MATLAB:FinalTimeNotProper',...
              'The final time T=%d should be more than 155 to pass the pulse',T);
        end
        %------------------------------------------------------------------
    case 'OUP'
        Fun=@fun3;
        load('b15OUP.mat');
        bb=@(t,N)b(t,N);
        % Check if the final time is suitable with produced OUP
        if T>700
            error('MATLAB:FinalTimeNotCompatible', ...
                'The final time T=%d should be less than 700', T);
        end        
end
B=b;

%--------------------------------------------------------------------------

switch InitialCondition
    case 'Uniform(0,0.1)'
        load('x0.mat');
        X0=x0;
        %------------------------------------------------------------------
    case 'Equilibrium'
        load('X0Equilibrium.mat');
        X0=x0;
        % Check the conditions for fixed points
         % Check the conditions fixed points
        if n==4 && all(Ki==ones(N,1)) && ...
            strcmp(Perturbation,'Pulse') && strcmp(Interaction,'Random')
        else
            warning('MATLAB:FixedPoints', ...
                ['Check the conditions for using InitialCondition=Equilibrium: '...
                '---> Interaction=Random, n=4, ki=1*ones(N,1), and Perturbation=Pulse']);
        end           
end

%----------------------------------------------------------------------

t0=0; % initial time
h=0.01; % step size for computing


% solver for fractional differential equation
[t, x] = fde_pi12_pc(mu,Fun,t0,T,X0,h);
end
% =========================================================================
% =========================================================================
function dx=fun(~,x)

global b N Ki

dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(b(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function dx=fun2(t,x)
global b N Ki 

%%% Pulse pertubation
if t>100 && t<155
    B=b;
B(1:7)=.5*b(1:7);
else
    B=b;
end
%%%

dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function dx=fun3(t,x)
global N Ki b

dx=zeros(N,1);
%%% OUP pertubation
    if t==0        
    B=b(1,1:N);        
    else
    B=b(ceil(t),1:N);        
    end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function fi=fi_Xk(i, x)
global n N Kij
fi=1;
K=1:N;K(i)=[];
for j=1:N-1
    k=K(j);
fi=fi*(Kij(i,k).^n/(Kij(i,k).^n+x(k).^n));
end
end
