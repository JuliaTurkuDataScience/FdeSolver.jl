function [t, y] = flmm2(alpha,fdefun,Jfdefun,t0,tfinal,y0,h,param,method,tol,itmax)
%FLMM2  Solves an initial value problem for a non-linear fractional
%       differential equation (FDEs) by means of some implicit fractional
%       linear multistep methods (FLMMs) of the second order.
%
%   FLMMs are a generalization to FDEs of classical linar multistep methods
%   and the were introduced in [1]. This code implements the 3 different 
%   implicit FLMMs of the second order described in [2]. The specific
%   method is selected by means of the parameter METHOD: 1 for the
%   generalization of the classical Trapezoidal rule (also known as the
%   Tustin method); 2 for the generalization of the Newton-Gregory formula;
%   3 for the generalization of the Backward Differentiation Formula (BDF).
%   By default the BDF is selected when no value to the parameter METHOD is
%   assigned. 
%
%   The discrete convolution quadrature rule are evaluated by means of the
%   FFT algorithm described in [3] allowing to keep the computational cost
%   proportional to N*log(N)^2 instead of N^2 as in the classical
%   implementation; N is the number of time-points in which the solution is
%   evaluated, i.e. N = (TFINAL-T)/H.  
%
%   Please, cite this code as [1,2]if you need.
% 
%   [T,Y] = FLMM2(ALPHA,FDEFUN,JFDEFUN,T0,TFINAL,h) integrates the initial
%   value problem for the FDE, or the system of FDEs, of order ALPHA > 0
%      D^ALPHA Y(t) = FDEFUN(T,Y(T))
%      Y^(k)(T0) = Y0(:,k+1), k=0,...,m-1
%   where m is the smallest integer grater than ALPHA and D^ALPHA is the
%   fractional derivative according to the Caputo's definition. FDEFUN is a
%   function handle corresponding to the vector field of the FDE for a the
%   scalar time variable T and the vector state variable Y (note that
%   FDEFUN(T,Y) must return a column vector). JFDEFUN is a function handle
%   corresponding to the Jacobian of the vector field FDEFUN (note that
%   JFDEFUN(T,Y) must return a matrix). The set of initial conditions Y0 is
%   a matrix with a number of rows equal to the size of the problem (hence
%   equal to the number of rows of the output of FDEFUN) and a number of
%   columns depending on ALPHA and given by m. The step-size H>0 is assumed 
%   constant throughout the integration.  
%
%   [T,Y] = FLMM2(ALPHA,FDEFUN,JFDEFUN,T0,TFINAL,H,PARAM) solves as above
%   with the additional set of parameters for FDEFUN and JFDEFUN as
%   FDEFUN(T,Y,PARAM) and JFDEFUN(T,Y,PARAM).
%
%   [T,Y] = FLMM2(ALPHA,FDEFUN,JFDEFUN,T0,TFINAL,H,PARAM,METHOD) allows to
%   select a specific method to be used for the integration (refer to [2]
%   for a description of each method). According to the value of METHOD one
%   of the following method is selected (the default value is 3): 
%      METHOD = 1 : Trapezoidal method
%      METHOD = 2 : Newton-Gregory formula
%      METHOD = 3 : BDF-2 
%
%   [T,Y] = FLMM2(ALPHA,FDEFUN,T0,TFINAL,H,PARAM,METHOD,TOL) use the
%   tolerance TOL when solving the internal nonlinear system of equations.
%   The defalut value for TOL is 1.0e-6.
%
%   [T,Y] = FLMM2(ALPHA,FDEFUN,T0,TFINAL,H,PARAM,METHOD TOL,ITMAX) allows 
%   to specify the maximum number of iterations for solving the internal
%   nonlinear system of equations. The defalut value for ITMAX is 100
%
%
% REFERENCES
%
%   [1] C. Lubich, Discretized fractional calculus, SIAM J. Numer. Anal.
%   17(3) (1986), 704–719 
%
%   [2] R. Garrappa, Trapezoidal methods for fractional differential
%   equations: Theoretical and computational aspects. Mathematics and
%   Computers in Simulation, in press,
%   http://dx.doi.org/10.1016/j.matcom.2013.09.012 
%
%   [3] E. Hairer, C. Lubich, M. Schlichte, Fast numerical solution of
%   nonlinear Volterra convolution equations, SIAM J. Sci. Statist. Comput.
%   6 (3) (1985) 532-541.
%
%
%
%   Copyright (c) 2014, Roberto Garrappa, University of Bari, Italy
%   roberto dot garrappa at uniba dot it
%   Homepage: http://www.dm.uniba.it/Members/garrappa
%   Revision: 1.0 - Date: June 27 2014

% Check inputs
if nargin < 11
    itmax = 100 ;
    if nargin < 10
        tol = 1.0e-6 ;
        if nargin < 9
            method = 3 ;
            if nargin < 8
                param = [] ;
            end
        end
    end
end

% Check order of the FDE
if alpha <= 0
    error('MATLAB:flmm2:NegativeOrder',...
        ['The order ALPHA of the FDE must be positive. The value ' ...
        'ALPHA = %f can not be accepted. See FLMM2.'], alpha);
end

% Check the step--size of the method
if h <= 0
    error('MATLAB:flmm2:NegativeStepSize',...
        ['The step-size H for the method must be positive. The value ' ...
        'H = %e can not be accepted. See FLMM2.'], h);
end

% Structure for storing initial conditions
ic.t0 = t0 ;
ic.y0 = y0 ;
ic.m_alpha = ceil(alpha) ;
ic.m_alpha_factorial = factorial(0:ic.m_alpha-1) ;

% Structure for storing information on the problem
Probl.ic = ic ;
Probl.fdefun = fdefun ;
Probl.Jfdefun = Jfdefun ;
Probl.problem_size = size(y0,1) ;
Probl.param = param ;

% Check number of initial conditions
if size(y0,2) < ic.m_alpha
    error('MATLAB:flmm2:NotEnoughInputs', ...
        ['A not sufficient number of assigned initial conditions. ' ...
        'Order ALPHA = %f requires %d initial conditions. See FLMM2.'], ...
        alpha,ic.m_alpha);
end

% Check compatibility size of the problem with size of the vector field
f_temp = f_vectorfield(t0,y0(:,1),Probl) ;
if Probl.problem_size ~= size(f_temp,1)
    error('MATLAB:flmm2:SizeNotCompatible', ...
        ['Size %d of the problem as obtained from initial conditions ' ...
        '(i.e. the number of rows of Y0) not compatible with the ' ...
        'size %d of the output of the vector field FDEFUN. ' ...
        'See FLMM2.'], Probl.problem_size,size(f_temp,1));
end

% Number of points in which to evaluate the solution or the weights
r = 16 ;
N = ceil((tfinal-t0)/h) ;
Nr = ceil((N+1)/r)*r ;
Q = ceil(log2((Nr)/r)) - 1 ;
NNr = 2^(Q+1)*r ;

% Preallocation of some variables
y = zeros(Probl.problem_size,N+1) ;
fy = zeros(Probl.problem_size,N+1) ;
zn = zeros(Probl.problem_size,NNr+1) ;

% Evaluation of convolution and starting weights of the FLMM
[flmm.omega,flmm.w,flmm.s] = Weights(alpha,NNr+1,method) ;
flmm.halpha = h^alpha ;
flmm.tol = tol ; flmm.itmax = itmax ;

% Initializing solution and proces of computation
t = (0 : N)*h ;
y(:,1) = y0(:,1) ;
fy(:,1) = f_vectorfield(t0,y0(:,1),Probl) ;
[y, fy] = FirstApproximations(t, y, fy, flmm, Probl) ; 
[y, fy] = Triangolo(flmm.s+1, r-1, 0, t, y, fy, zn, N, flmm, Probl) ;

% Main process of computation by means of the FFT algorithm
nx0 = 0 ; ny0 = 0 ; 
ff = zeros(1,2^(Q+2),1) ; ff(1:2) = [0,2] ;
for q = 0 : Q
    L = 2^q ;
    [y, fy] = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, ...
        zn, N, flmm, Probl ) ;
    ff(1:4*L) = [ff(1:2*L) , ff(1:2*L-1) , 4*L] ;
end


% Evaluation solution in TFINAL when TFINAL is not in the mesh
if tfinal < t(N+1)
    c = (tfinal - t(N))/h ;
    t(N+1) = tfinal ;
    y(:,N+1) = (1-c)*y(:,N) + c*y(:,N+1) ;
end

t = t(1:N+1) ; y = y(:,1:N+1) ;

end


% =========================================================================
% =========================================================================
% r : dimension of the basic square or triangle
% L : factor of resizing of the squares
function [y, fy] = DisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, ...
    zn, N , flmm, Probl)

nxi = nx0 ; nxf = nx0 + L*r - 1 ;
nyi = ny0 ; nyf = ny0 + L*r - 1 ;
is = 1 ;
s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;

i_triangolo = 0 ; stop = 0 ;
while ~stop
    
    stop = nxi+r-1 == nx0+L*r-1 | (nxi+r-1>=Nr-1) ; 
    
    zn = Quadrato(nxi, nxf, nyi, nyf, fy, zn, flmm, Probl) ;
    
    [y, fy] = Triangolo(nxi, nxi+r-1, nxi, t, y, fy, zn, N, flmm, Probl) ;
    i_triangolo = i_triangolo + 1 ;
    
    if ~stop
        if nxi+r-1 == nxf   % Il triangolo finisce dove finisce il quadrato -> si scende di livello
            i_Delta = ff(i_triangolo) ;
            Delta = i_Delta*r ;
            nxi = s_nxf(is)+1 ; nxf = s_nxf(is)  + Delta ;
            nyi = s_nxf(is) - Delta +1; nyf = s_nxf(is)  ;
            s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;
        else % Il triangolo finisce prima del quadrato -> si fa un quadrato accanto
            nxi = nxi + r ; nxf = nxi + r - 1 ; nyi = nyf + 1 ; nyf = nyf + r  ;
            is = is + 1 ;
            s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;
        end
    end
    
end

end

% =========================================================================
% =========================================================================
function zn = Quadrato(nxi, nxf, nyi, nyf, fy, zn, flmm, Probl)

coef_beg = nxi-nyf ; coef_end = nxf-nyi+1 ;
funz_beg = nyi+1 ; funz_end = nyf+1 ;

vett_coef = flmm.omega(coef_beg+1:coef_end+1) ;
vett_funz = [ fy(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
zzn = real(FastConv(vett_coef,vett_funz)) ;
zn(:,nxi+1:nxf+1) = zn(:,nxi+1:nxf+1) + zzn(:,nxf-nyf:end-1) ;

end



% =========================================================================
% =========================================================================
function [y, fy] = Triangolo(nxi, nxf, j0, t, y, fy, zn, N, flmm, Probl)

%if nxi < 20 , y(:,1:8) , fy(:,1:8) , pause , end

for n = nxi : min(N,nxf)
    
    n1 = n + 1 ; 
    St = StartingTerm(t(n1),Probl.ic) ;
    
    Phi = zeros(Probl.problem_size,1) ;
    for j = 0 : flmm.s
        Phi = Phi + flmm.w(j+1,n1)*fy(:,j+1) ;
    end
    for j = j0 : n-1
        Phi = Phi + flmm.omega(n-j+1)*fy(:,j+1) ;
    end
    Phi_n = St + flmm.halpha*(zn(:,n1) + Phi) ;
    
    yn0 = y(:,n) ; fn0 = f_vectorfield(t(n1),yn0,Probl) ;
    Jfn0 = Jf_vectorfield(t(n1),yn0,Probl) ;
    Gn0 = yn0 - flmm.halpha*flmm.omega(1)*fn0 - Phi_n ;
    stop = 0 ; it = 0 ;
    while ~stop
        
        JGn0 = eye(Probl.problem_size,Probl.problem_size) - ...
            flmm.halpha*flmm.omega(1)*Jfn0 ;
        yn1 = yn0 - JGn0\Gn0 ;
        fn1 = f_vectorfield(t(n1),yn1,Probl) ;
        Gn1 = yn1 - flmm.halpha*flmm.omega(1)*fn1 - Phi_n ;
        it = it + 1 ;
        
        stop = norm(yn1-yn0,inf) < flmm.tol | norm(Gn1,inf) <  flmm.tol;
        if it > flmm.itmax && ~stop
            warning('MATLAB:flmm2:NonConvergence',...
                strcat('Internal iterations do not convergence to the ', ...
                ' tolerance %9.2e in %d iterations. Try with a less ', ...
                ' restrictive tolarence or a bigger numeber of max ', ... 
                ' iterations. Reducing the step-size might also help.'), ...
                flmm.tol,flmm.itmax) ;
            stop = 1 ;
        end
        
        yn0 = yn1 ; Gn0 = Gn1 ;
        if ~stop
            Jfn0 = Jf_vectorfield(t(n1),yn0,Probl) ;
        end
        
    end
    y(:,n1) = yn1 ;
    fy(:,n1) = fn1 ;
end

end

% =========================================================================
% =========================================================================
function [y, fy] = FirstApproximations(t, y, fy, flmm, Probl)

m = Probl.problem_size ; s = flmm.s ;
Im = eye(m,m) ; Ims = eye(m*s,m*s) ;

Y0 = zeros(s*m,1) ; F0 = Y0 ; B0 = Y0 ;
for j = 1 : s
    Y0((j-1)*m+1:j*m,1) = y(:,1) ;
    F0((j-1)*m+1:j*m,1) = f_vectorfield(t(j+1),y(:,1),Probl) ;
    St = StartingTerm(t(j+1),Probl.ic) ;
    B0((j-1)*m+1:j*m,1) = ...
        St + flmm.halpha*(flmm.omega(j+1)+flmm.w(1,j+1))*fy(:,1) ;
end

W = zeros(s,s) ;
for i = 1 : s
    for j = 1 : s
        if i >= j
            W(i,j) = flmm.omega(i-j+1) + flmm.w(j+1,i+1) ;
        else
            W(i,j) = flmm.w(j+1,i+1) ;
        end
    end
end
W = flmm.halpha*kron(W,Im) ; 

G0 = Y0 - B0 - W*F0 ;
JF = zeros(s*m,s*m) ;
for j = 1 : s
    JF((j-1)*m+1:j*m,(j-1)*m+1:j*m) = Jf_vectorfield(t(j+1),y(:,1),Probl) ;
end

stop = 0 ; it = 0 ;
while ~stop
    
    JG = Ims - W*JF ;
    Y1 = Y0 - JG\G0 ;
    
    for j = 1 : s
        F1((j-1)*m+1:j*m,1) = f_vectorfield(t(j+1),Y1((j-1)*m+1:j*m,1),Probl) ;
    end
    G1 = Y1 - B0 - W*F1 ;
    
    it = it + 1 ;
    
    stop = norm(Y1-Y0,inf) < flmm.tol | norm(G1,inf) <  flmm.tol ;
    if it > flmm.itmax && ~stop
        warning('MATLAB:flmm2:NonConvegence',...
            strcat('Internal iterations do not convergence to the tolerance ', ...
            '%e in %d iterations. Try with a larger tolarence of a bigger',...
            'numeber of max iterations'),flmm.tol,flmm.itmax) ;
        stop = 1 ;
    end
    
    Y0 = Y1 ; G0 = G1 ;
    if ~stop
        for j = 1 : s
            JF((j-1)*m+1:j*m,(j-1)*m+1:j*m) = Jf_vectorfield(t(j+1),Y1((j-1)*m+1:j*m,1),Probl) ;
        end
    end
    
end

for j = 1 : s
    y(:,j+1) = Y1((j-1)*m+1:j*m,1) ;
    fy(:,j+1) = F1((j-1)*m+1:j*m,1) ;
end

end



% =========================================================================
% =========================================================================
function z = FastConv(x, y)

Lx = length(x) ; Ly = size(y,2) ; problem_size = size(y,1) ;
if Lx ~= Ly
    disp('Warning: dimension of 1x and y in FastConv does not agree') ;
    %z = 0 ; return
end
r = Lx ;


z = zeros(problem_size,r) ;
X = fft(x,r) ;
for i = 1 : problem_size
    Y = fft(y(i,:),r) ;
    Z = X.*Y ;
    z(i,:) = ifft(Z,r) ;
end

end


function [omega,w,s] = Weights(alpha,N,method)
% =========================================================================
% Evaluate convolution and starting weights for fractional linear multistep
% methods of order 2 for approximating a Riemann-Liouville integral of
% order alpha. Convolution weights are evaluated as the coefficients in the
% formal power series of the generating function. Starting weights are
% evaluated by solving the linear systems resulting from imposing that the
% formula is exact for some fractional powers.
% =========================================================================
% alpha  : fractional order of the integral
% N      : number of coefficients to be evaluated
% method : convolution quadrature to be implemented
%          1 : Trapezoidal method
%          2 : Newton-Gregory formula
%          3 : BDF-2

% =========================================================================
% Evaluation of convolution weights
switch method
    
    % Trapezoidal method with generating function ((1+x)/2/(1-x))^alpha
    case 1 
        omega1 = zeros(1,N+1) ; omega2 = omega1 ;
        omega1(1) = 1 ; omega2(1) = 1 ;
        alpha_minus_1 = alpha - 1 ; alpha_plus_1 = alpha + 1 ;
        for n = 1 : N
            omega1(n+1) = (alpha_plus_1/n - 1)*omega1(n) ;
            omega2(n+1) = (1 + alpha_minus_1/n)*omega2(n) ;
        end
        x = fft([omega1,zeros(size(omega1))]) ;
        y = fft([omega2,zeros(size(omega2))]) ;
        omega = ifft(x.*y) ; 
        omega = omega(1:N+1)/2^alpha ;

    % Newton-Gregory formula with generating function (1-x)^(-alpha)*(1-alpha/2*(1-x))
    case 2 
        omega1 = zeros(1,N+1) ; omega = omega1 ;
        alphameno1 = alpha - 1 ;
        omega1(1) = 1 ;
        for n = 1 : N
            omega1(n+1) = (1 + alphameno1/n)*omega1(n) ;
        end
        omega(1) = 1-alpha/2 ;
        omega(2:N+1) = (1-alpha/2)*omega1(2:N+1) + alpha/2*omega1(1:N) ;
     
	% BDF-2 with generating function (2/3/(1-4x/3+x^2/3))^alpha
    case 3 
        omega = zeros(1,N+1) ;
        onethird = 1/3 ; fourthird = 4/3 ;
        twothird_oneminusalpha = 2/3*(1-alpha) ;
        fourthird_oneminusalpha = 4/3*(1-alpha) ;
        omega(1) = 1 ; omega(2) = fourthird*alpha*omega(1) ;
        for n = 2 : N
            omega(n+1) = (fourthird - fourthird_oneminusalpha/n)*omega(n) + ...
                (twothird_oneminusalpha/n - onethird)*omega(n-1) ;
        end
        omega = omega*((2/3)^(alpha)) ;      
end

% =========================================================================
% Evaluation of starting weights

% Set of the real powers
k = floor(1/abs(alpha)) ;
if abs(k - 1/alpha) < 1.0e-12
    A = (0:k)*abs(alpha) ;
else
    A = [(0:k)*abs(alpha), 1] ;
end
s = length(A) - 1 ;

% Generation of the matrix and the right hand--side vectors of the system
nn = 0 : N ;
V = zeros(s+1,s+1) ; jj_nu = zeros(s+1,N+1) ; nn_nu_alpha = jj_nu ;
for i = 0 : s
    nu = A(i+1) ;
    V(i+1,:) = (0:s).^nu ;
    jj_nu(i+1,:) = nn.^nu ;
    if alpha > 0
        nn_nu_alpha(i+1,:) = gamma(nu+1)/gamma(nu+1+alpha)*nn.^(nu+alpha) ;
    else
        if i == 0
            nn_nu_alpha(i+1,:) = zeros(1,N+1) ;
        else
            nn_nu_alpha(i+1,:) = gamma(nu+1)/gamma(nu+1+alpha)*nn.^(nu+alpha) ;
        end
    end
end
temp = FastConv([omega,zeros(size(omega))],[jj_nu,zeros(size(jj_nu))]) ;
b = nn_nu_alpha - temp(:,1:N+1) ;

% Solution of the linear system with multiple right-hand side
w = V\b ;

end


% =========================================================================
% =========================================================================
function f = f_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.fdefun,t,y) ;
else
    f = feval(Probl.fdefun,t,y,Probl.param) ;
end

end

% =========================================================================
% =========================================================================
function f = Jf_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.Jfdefun,t,y) ;
else
    f = feval(Probl.Jfdefun,t,y,Probl.param) ;
end

end

% =========================================================================
% =========================================================================
function ys = StartingTerm(t,ic)

ys = zeros(size(ic.y0,1),1) ;
for k = 1 : ic.m_alpha
    ys = ys + (t-ic.t0)^(k-1)/ic.m_alpha_factorial(k)*ic.y0(:,k) ;
end

end