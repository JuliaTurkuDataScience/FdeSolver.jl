function [t, y] = MT_FDE_PI1_Im(al,lam,f_fun,J_fun,t0,T,y0,h,param,tol,itmax)
% -------------------------------------------------------------------------
% MT_FDE_PI1_Im  Solves on the interval [t0,T] an initial value problem
%                for the multiterm fractional differential equation
%
%    lam_Q D^(al_Q) y(t) + ... + lam_1 D^(al_1) y(t) = f(t,y(t))
%    y(0) = y0(1), y'(0) = y0(2), ... y^m(0) = y0(m)
%
%                with m the smallest integer greater than max(al_1,...,al_Q). 
%                The problem is solved by means of the implicit  
%                product-integration rule of rectangular type having order 
%                1 of convergence. 
%
% Further information on this code are available in [1]; please, cite this
% code as [1] (see later on). 
%
% -------------------------------------------------------------------------
% Description of input parameters
% -------------------------------------------------------------------------
% al        : vector of the fractional orders alpha of each term of the
%             multiterm FDE. The orders in al = (al_1, ..., al_Q) must not
%             be necessarly sorted in descending or ascending order.
% lam       : vector of the coefficients lam of each term of the multiterm 
%             FDE which correspondis to the respective al-order drivatives
%             lam = (lam_1, ..., lam_Q). It is requested that the element
%             of lam corresponding to max(al_1,...,al_Q) is not zero.
% f_fun     : function defining the source or nonlinear term f(t,y(t)).
%             f_fun must be a function handle corresponding to the vector
%             field of the multiterm FDE for the scalar variable t and the
%             (possibly vector) state variable Y as independent variable.
%             If necessary, f_fun(t,y,param) can use a vector of parameters
%             (see later on). Note that f_fun(t,y) or f_fun(t,y,param) must
%             return a scalar when the f_fun is a scalar function or a
%             column vector when f_fun is a vector function.
% J_fun     : Jacobian of the the source or nonlinear term f(t,y(t)). J_fun
%             must be a function handle corresponding to the Jacobian of
%             the vector field f_fun. Note that J_fun(t,y) must return a
%             matrix. If necessary, J_fun(t,y,param) can use a vector of
%             parameters (see later on).  
% t0, T     : starting and final time of the integration interval.
% y0        : initial conditions to be assigned to uniquely define the
%             solution. The set of initial conditions y0 must be a matrix
%             with a number of rows equal to the size of the problem (hence
%             equal to the number of rows of the output of f_fun) and a
%             number of columns equal to the smallest integer greater than
%             max(al_1,...,al_Q).  
% h         : step-size for integration and must be real and positive
% param     : vector of possible parameters for the evaluation of the
%             source or nonlinear term f(t,y(t)) and its Jacobian. This an
%             optional parameter: it can be omitted or it is possible to
%             use the empty vector [] if the function f_fun(t,y) and its
%             Jacobian do not have any parameter. 
% tol       : tolerance for the solution of the internal nonlinear system
%             of equations. This is an optional parameter and when not
%             specified the defalut value 1.0e-6 is assumed. 
% itmax     : number of maximum iterations for the solution of the internal
%             nonlinear system of equations. This is an optional parameter
%             and when not specifiad the defalut value 100 is assumed. 
%
% -------------------------------------------------------------------------
% Description of output parameters
% -------------------------------------------------------------------------
% t          : vector of nodes on the interval [t0,T] in which the
%              numerical solution is evaluated
% y          : matrix with the values of the solution evaluated in t
%
% -------------------------------------------------------------------------
% Possible usages
% -------------------------------------------------------------------------
% [t, y] = MT_FDE_PI1_Im(al,lam,f_fun,J_fun,t0,T,y0,h)
% [t, y] = MT_FDE_PI1_Im(al,lam,f_fun,J_fun,t0,T,y0,h,param)
% [t, y] = MT_FDE_PI1_Im(al,lam,f_fun,J_fun,t0,T,y0,h,param,tol)
% [t, y] = MT_FDE_PI1_Im(al,lam,f_fun,J_fun,t0,T,y0,h,param,tol,itmax)
%
% -------------------------------------------------------------------------
% References and other information
% -------------------------------------------------------------------------
%
%  [1] Garrappa R.: Numerical Solution of Fractional Differential
%  Equations: a Survey and a Software Tutorial, Mathematics 2018, 6(2), 16
%  doi: https://doi.org/10.3390/math6020016
%  downloadable pdf: http://www.mdpi.com/2227-7390/6/2/16/pdf   
%
%  Revision: 1.0 - Date: November, 23 2017
%
%  Please, report any problem or comment to :
%          roberto dot garrappa at uniba dot it
%
%  Copyright (c) 2017
%
%  Author:
%   Roberto Garrappa (University of Bari, Italy)
%   roberto dot garrappa at uniba dot it
%   Homepage: http://www.dm.uniba.it/Members/garrappa
%
% -------------------------------------------------------------------------

% Check inputs
if nargin < 11
    itmax = 100 ;
    if nargin < 10
        tol = 1.0e-6 ;
        if nargin < 9
            param = [] ;
            if nargin < 8
                error('MATLAB:NumberOfParameters',...
                    'The number of input parameters is not sufficient.') ;
            end
        end
    end
end

% Check order of the multiterm FDE 
if any(al < 0)
    error('MATLAB:NegativeOrder',...
        ['The orders ALPHA of the muliterm FDE must be positive. On of the' ...
         'values can not be accepted.']);
end 

% Check the step--size of the method
if h <= 0
    error('MATLAB:NegativeStepSize',...
        ['The step-size H for the method must be positive. The value ' ...
         'H = %e can not be accepted.'], h);
end 

% Ascending ordering (w.r.t. to the fractional order) of the terms of the
% equation and transformation in normal form
Q = length(al) ;
[al,i_al] = sort(al,'ascend') ;
lam = lam(i_al) ;
al_Q = al(end) ; al_i = al(1:end-1) ;
lam_Q = lam(end) ; lam_rat_i = lam(1:end-1)/lam_Q ;
m_Q = ceil(al_Q) ; m_i = ceil(al(1:end-1)) ;
beta = [ al_Q - al_i , al_Q ] ;

% Structure for storing initial conditions
ic.t0 = t0 ;
ic.problem_size = size(y0,1) ; ic.y0 = y0 ;
ic.Q = Q ; ic.m_Q = m_Q ; ic.m_i = m_i ; ic.beta = beta ;
ic.lam_rat_i = lam_rat_i ;
gamma_val = zeros(Q,m_Q) ;
for i = 1 : Q-1
    k = 0 : m_i(i)-1 ;
    gamma_val(i,k+1) = gamma(k+beta(i)+1) ;
end
k = 0 : m_Q-1 ; gamma_val(Q,:) = factorial(k) ;
ic.gamma_val = gamma_val ;

% Structure for storing information on the problem
Probl.ic = ic ;
Probl.fdefun = f_fun ;
Probl.Jfdefun = J_fun ;
Probl.problem_size = size(y0,1) ;
Probl.param = param ;
Probl.Q = Q ;
Probl.lam_Q = lam_Q ;
Probl.lam_rat_i = lam_rat_i ;

% Check number of initial conditions
if size(y0,2) < ic.m_Q
    error('MATLAB:NotEnoughInputs', ...
        ['A not sufficient number of initial conditions are assigned.' ...
        'Y0 must have as many columns as the number of derivatives at the ' ...
        'origin involved by initial conditions (%d for this problem).'], ...
        ic.m_Q);
end

% Check compatibility size of the problem with size of the vector field
f_temp = f_vectorfield(t0,y0(:,1),Probl) ;
if Probl.problem_size ~= size(f_temp,1)
    error('MATLAB:SizeNotCompatible', ...
        ['Size %d of the problem as obtained from initial conditions ' ...
        '(i.e. the number of rows of Y0) not compatible with the ' ...
        'size %d of the output of the vector field F_FUN. '], ...
        Probl.problem_size,size(f_temp,1));
end

% Number of points in which to evaluate the solution or the weights
r = 64 ; 
N = ceil((T-t0)/h) ;
Nr = ceil((N+1)/r)*r ;
Qr = ceil(log2((Nr)/r)) - 1 ;
NNr = 2^(Qr+1)*r ;

% Preallocation of some variables
y = zeros(Probl.problem_size,N+1) ;
fy = zeros(Probl.problem_size,N+1) ;
zn = zeros(Probl.problem_size,NNr+1,Q) ;

% Evaluation of weights of the method
nvett = 0 : NNr+1 ; 
bn = zeros(Q,NNr+1) ; 
for i = 1 : Q
    nbeta = nvett.^beta(i) ; 
    bn(i,:) = ( nbeta(2:end) - nbeta(1:end-1) ) * h^beta(i) / gamma(beta(i)+1) ; 
end
C = 0 ;
for i = 1 : Q-1
    C = C + lam_rat_i(i)*bn(i,1) ;
end
METH.bn = bn ; METH.C = C ;
METH.itmax = itmax ; METH.tol = tol ; 

% Initializing solution and proces of computation
t = (0 : N)*h ;
y(:,1) = y0(:,1) ;
fy(:,1) = f_vectorfield(t0,y0(:,1),Probl) ;
[y, fy] = Triangolo(1, r-1, t, y, fy, zn, N, METH, Probl ) ;

% Main process of computation by means of the FFT algorithm
ff = zeros(1,2.^(Qr+2)) ; ff(1:2) = [0 2] ; card_ff = 2 ;
nx0 = 0 ; ny0 = 0 ;
for qr = 0 : Qr
    L = 2^qr ; 
    [y, fy] = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, ...
        zn, N, METH, Probl ) ;
    ff(1:2*card_ff) = [ff(1:card_ff) ff(1:card_ff)] ; 
    card_ff = 2*card_ff ; 
    ff(card_ff) = 4*L ; 
end

% Evaluation solution in TFINAL when TFINAL is not in the mesh
if T < t(N+1)
    c = (T - t(N))/h ;
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
    zn, N , METH, Probl)

nxi = nx0 ; nxf = nx0 + L*r - 1 ;
nyi = ny0 ; nyf = ny0 + L*r - 1 ;
is = 1 ;
s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;

i_triangolo = 0 ; stop = 0 ;
while ~stop
    
    stop = nxi+r-1 == nx0+L*r-1 | (nxi+r-1>=Nr-1) ;
    
    zn = Quadrato(nxi, nxf, nyi, nyf, y, fy, zn, METH, Probl) ;
    
    [y, fy] = Triangolo(nxi, nxi+r-1, t, y, fy, zn, N, METH, Probl) ;
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
function zn = Quadrato(nxi, nxf, nyi, nyf, y, fy, zn, METH, Probl)

coef_beg = nxi-nyf ; coef_end = nxf-nyi+1 ;
funz_beg = nyi+1 ; funz_end = nyf+1 ;

% Evaluation convolution segments for each term
Q = Probl.Q ;
for i = 1 : Q
    vett_coef = METH.bn(i,coef_beg:coef_end) ;
    if nyi == 0 % Evaluation of the lowest square
        if i < Q
            vett_funz = [zeros(Probl.problem_size,1) , ...
                y(:,funz_beg+1:funz_end) , ...
                zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
        else
            vett_funz = [zeros(Probl.problem_size,1) , ...
                fy(:,funz_beg+1:funz_end) , ...
                zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
        end
    else % Evaluation of any square but not the lowest
        if i < Q
            vett_funz = [ y(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
        else
            vett_funz = [ fy(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
        end
    end
    zzn = real(FastConv(vett_coef,vett_funz)) ;
    zn(:,nxi+1:nxf+1,i) = zn(:,nxi+1:nxf+1,i) + zzn(:,nxf-nyf+1:end) ;
end

end


% =========================================================================
% =========================================================================
function [y, fy] = Triangolo(nxi, nxf, t, y, fy, zn, N, METH, Probl)

Q = Probl.Q ;

for n = nxi : min(N,nxf)
    
    St = StartingTerm_Multi(t(n+1),Probl.ic) ;
    
    Phi_n = St ;
    for i = 1 : Q-1
        temp = zn(:,n+1,i) ; 
        for j = nxi : n-1
            temp = temp + METH.bn(i,n-j+1)*y(:,j+1) ;
        end
        Phi_n = Phi_n - Probl.lam_rat_i(i)*temp ;
    end
    temp = zn(:,n+1,Q) ;
    for j = nxi : n-1
        temp = temp + METH.bn(Q,n-j+1)*fy(:,j+1) ;
    end
    Phi_n = Phi_n + temp/Probl.lam_Q ;
    
    yn0 = y(:,n) ; fn0 = f_vectorfield(t(n+1),yn0,Probl) ;
    Jfn0 = Jf_vectorfield(t(n+1),yn0,Probl) ;
    Gn0 = (1+METH.C)*yn0 - METH.bn(Q,1)/Probl.lam_Q*fn0 - Phi_n ;
    stop = 0 ; it = 0 ;
    
    while ~stop
        
        JGn0 = (1+METH.C)*eye(Probl.problem_size,Probl.problem_size) - ...
            METH.bn(Q,1)/Probl.lam_Q*Jfn0 ;
        yn1 = yn0 - JGn0\Gn0 ;
        fn1 = f_vectorfield(t(n+1),yn1,Probl) ;
        Gn1 = (1+METH.C)*yn1 - METH.bn(Q,1)/Probl.lam_Q*fn1 - Phi_n ;
        it = it + 1 ;
        
        stop = norm(yn1-yn0,inf) < METH.tol | norm(Gn1,inf) <  METH.tol;
        if it > METH.itmax && ~stop
            warning('MATLAB:NonConvegence',...
                strcat('Internal iterations do not convergence to the tolerance ', ...
                '%e in %d iterations. Try with a larger tolarence of a bigger',...
                'numeber of max iterations'),METH.tol,METH.itmax) ;
            stop = 1 ;
        end
        
        yn0 = yn1 ; Gn0 = Gn1 ;
        if ~stop
            Jfn0 = Jf_vectorfield(t(n+1),yn0,Probl) ;
        end
        
    end
    
    y(:,n+1) = yn1 ;
    fy(:,n+1) = fn1 ;
    
end

end


% =========================================================================
% =========================================================================
function z = FastConv(x, y)

Lx = length(x) ; Ly = size(y,2) ; problem_size = size(y,1) ;
if Lx ~= Ly
    disp('Warning: dimensions in FastConv do not agree') ;
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
function ys = StartingTerm_Multi(t,ic)

Q = ic.Q ; 
ys = zeros(ic.problem_size,1) ;
for k = 0 : ic.m_Q-1
    ys = ys + (t-ic.t0)^(k)/ic.gamma_val(Q,k+1)*ic.y0(:,k+1) ;
end
for i = 1 : ic.Q-1
    temp = zeros(ic.problem_size,1) ;
    for k = 0 : ic.m_i(i)-1
        temp = temp + (t-ic.t0)^(k+ic.beta(i))/ic.gamma_val(i,k+1)*ic.y0(:,k+1) ;
    end
    ys = ys + ic.lam_rat_i(i)*temp ;
end

end



