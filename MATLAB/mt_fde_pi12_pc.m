function [t, y] = MT_FDE_PI12_PC(al,lam,f_fun,t0,T,y0,h,param,mu,mu_tol)
% -------------------------------------------------------------------------
% MT_FDE_PI2_PC   Solves on the interval [t0,T] an initial value problem
%                 for the multiterm fractional differential equation
%
%    lam_Q D^(al_Q) y(t) + ... + lam_1 D^(al_1) y(t) = f(t,y(t))
%    y(0) = y0(1), y'(0) = y0(2), ... y^m(0) = y0(m)
%
%                 with m the smallest integer greater than max(al_1,...,al_Q). 
%                 The problem is solved by means of product-integration
%                 rules of rectangular and trapezoidal type implemented in 
%                 a predictor-corrector framework.
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
% t0, T     : starting and final time of the integration interval.
% y0        : initial conditions to be assigned to uniquely define the
%             solution. The set of initial conditions y0 must be a matrix
%             with a number of rows equal to the size of the problem (hence
%             equal to the number of rows of the output of f_fun) and a
%             number of columns equal to the smallest integer greater than
%             max(al_1,...,al_Q).  
% h         : step-size for integration and must be real and positive
% param     : vector of possible parameters for the evaluation of the
%             source or nonlinear term f_fun(t,y,param). This an optional
%             parameter: it can be omitted or it is possible to use the
%             empty vector [] if the function f_fun(t,y) does not require 
%             any parameter. 
% mu        : number of multiple corrector iterations (optional). The
%             following values for MU are admissible (default value mu=1):
%        mu = 0 : corrector not evaluated; the solution is provided just by 
%                 the predictor method (the first order rectangular rule);
%        mu > 0 : the corrector is evaluated by the selected number mu of
%                 times (the classical PECE method is obtained for mu=1); 
%      mu = Inf : the corrector is evaluated for a certain number of times
%                 until convergence of the iterations is reached (for
%                 convergence the difference between two consecutive 
%                 iterates is tested).
% mu_tol     : tolerance for testing convergence of corrector iterations
%              when mu = Inf (optional). If not specified, the default
%              value mu_tol = 1.E-6 is used.  
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
% [t, y] = MT_FDE_PI12_PC(al,lam,f_fun,t0,T,y0,h)
% [t, y] = MT_FDE_PI12_PC(al,lam,f_fun,t0,T,y0,h,param)
% [t, y] = MT_FDE_PI12_PC(al,lam,f_fun,t0,T,y0,h,param,mu)
% [t, y] = MT_FDE_PI12_PC(al,lam,f_fun,t0,T,y0,h,param,mu,mu_tol)
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
%  Version of November, 23 2017
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
if nargin < 10
    mu_tol = 1.0e-6 ;
    if nargin < 9
        mu = 1 ;
        if nargin < 8
            param = [] ;
            if nargin < 7
                error('MATLAB:NumberOfParameters',...
                    'The number of input parameters is not sufficient') ;
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
zn_pred = zeros(Probl.problem_size,NNr+1,Q) ;
if mu > 0
    zn_corr = zeros(Probl.problem_size,NNr+1,Q) ;
else
    zn_corr = 0 ;
end

% Evaluation of weights of the method
nvett = 0 : NNr+1 ; 
bn = zeros(Q,NNr+1) ; an = zeros(Q,NNr+1) ; a0 = zeros(Q,NNr+1) ; 
for i = 1 : Q
    nbeta = nvett.^beta(i) ; nbeta1 = nbeta.*nvett ;
    bn(i,:) = ( nbeta(2:end) - nbeta(1:end-1) ) * h^beta(i) / gamma(beta(i)+1) ; 
    an(i,:) = [ 1 , (nbeta1(1:end-2) - 2*nbeta1(2:end-1) + nbeta1(3:end)) ] ...
        * h^beta(i) / gamma(beta(i)+2) ;
    a0(i,:) = [ 0 , nbeta1(1:end-2) - nbeta(2:end-1).*(nvett(2:end-1)-beta(i)-1)] ...
        * h^beta(i) / gamma(beta(i)+2) ;
end
C = 0 ;
for i = 1 : Q-1
    C = C + lam_rat_i(i)*an(i,1) ;
end
METH.an = an ; METH.a0 = a0 ; METH.bn = bn ; METH.C = C ;
METH.mu = mu ; METH.mu_tol = mu_tol ;  

% Initializing solution and proces of computation
t = (0 : N)*h ;
y(:,1) = y0(:,1) ;
fy(:,1) = f_vectorfield(t0,y0(:,1),Probl) ;
[y, fy] = Triangolo(1, r-1, t, y, fy, zn_pred, zn_corr, N, METH, Probl ) ;

% Main process of computation by means of the FFT algorithm
ff = zeros(1,2.^(Qr+2)) ; ff(1:2) = [0 2] ; card_ff = 2 ;
nx0 = 0 ; ny0 = 0 ;
for qr = 0 : Qr
    L = 2^qr ; 
    [y, fy] = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, ...
        zn_pred, zn_corr, N, METH, Probl ) ;
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
    zn_pred, zn_corr, N , METH, Probl)

nxi = nx0 ; nxf = nx0 + L*r - 1 ;
nyi = ny0 ; nyf = ny0 + L*r - 1 ;
is = 1 ;
s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;

i_triangolo = 0 ; stop = 0 ;
while ~stop
    
    stop = nxi+r-1 == nx0+L*r-1 | (nxi+r-1>=Nr-1) ;
    
    [zn_pred, zn_corr] = Quadrato(nxi, nxf, nyi, nyf, y, fy, zn_pred, zn_corr, METH, Probl) ;
    
    [y, fy] = Triangolo(nxi, nxi+r-1, t, y, fy, zn_pred, zn_corr, N, METH, Probl) ;
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
function [zn_pred, zn_corr] = Quadrato(nxi, nxf, nyi, nyf, y, fy, zn_pred, zn_corr, METH, Probl)

coef_beg = nxi-nyf ; coef_end = nxf-nyi+1 ;
funz_beg = nyi+1 ; funz_end = nyf+1 ;
Q = Probl.Q ;

% Evaluation convolution segments for the predictor
for i = 1 : Q
    vett_coef = METH.bn(i,coef_beg:coef_end) ;
    if i < Q
        vett_funz = [y(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
    else
        vett_funz = [fy(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
    end
    zzn_pred = real(FastConv(vett_coef,vett_funz)) ;
    zn_pred(:,nxi+1:nxf+1,i) = zn_pred(:,nxi+1:nxf+1,i) + zzn_pred(:,nxf-nyf+1-1:end-1) ;
end

% Evaluation convolution segments for the corrector
if METH.mu > 0
    for i = 1 : Q
        vett_coef = METH.an(i,coef_beg:coef_end) ;
        if nyi == 0 % Evaluation of the lowest square
            if i < Q
                vett_funz = [zeros(Probl.problem_size,1) , y(:,funz_beg+1:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
            else
                vett_funz = [zeros(Probl.problem_size,1) , fy(:,funz_beg+1:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
            end
        else % Evaluation of any square but not the lowest
            if i < Q
                vett_funz = [ y(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
            else
                vett_funz = [ fy(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
            end
        end
        zzn_corr = real(FastConv(vett_coef,vett_funz)) ;
        zn_corr(:,nxi+1:nxf+1,i) = zn_corr(:,nxi+1:nxf+1,i) + zzn_corr(:,nxf-nyf+1:end) ;
    end
else
    zn_corr = 0 ;
end

end


% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
function [y, fy] = Triangolo(nxi, nxf, t, y, fy, zn_pred, zn_corr, N, METH, Probl)

Q = Probl.Q ;

for n = nxi : min(N,nxf)
    
    St = StartingTerm_Multi(t(n+1),Probl.ic) ;
    
    % Evaluation of the predictor
    Phi_n = St ;
    if nxi == 1 % Case of the first triangle
        j_beg = 0 ;
    else % Case of any triangle but not the first
        j_beg = nxi ;
    end
    for i = 1 : Q-1
        temp = zn_pred(:,n+1,i) ;
        for j = j_beg : n-1
            temp = temp + METH.bn(i,n-j)*y(:,j+1) ;
        end
        Phi_n = Phi_n - Probl.lam_rat_i(i)*temp ;
    end
    temp = zn_pred(:,n+1,Q) ;
    for j = j_beg : n-1
        temp = temp + METH.bn(Q,n-j)*fy(:,j+1) ;
    end
    Phi_n = Phi_n + temp/Probl.lam_Q ;
    y_pred = Phi_n ;
    f_pred = f_vectorfield(t(n+1),y_pred,Probl) ;     
    
    % Evaluation of the corrector
    if METH.mu == 0
        y(:,n+1) = y_pred ;
        fy(:,n+1) = f_pred ;
    else
        j_beg = nxi ;
        Phi_n = St ;
        for i = 1 : Q-1
            temp = METH.a0(i,n+1)*y(:,1) + zn_corr(:,n+1,i) ;
            for j = j_beg : n-1
                temp = temp + METH.an(i,n-j+1)*y(:,j+1) ;
            end
            Phi_n = Phi_n - Probl.lam_rat_i(i)*temp ;
        end
        temp = METH.a0(Q,n+1)*fy(:,1) + zn_corr(:,n+1,Q) ;
        for j = j_beg : n-1
            temp = temp + METH.an(Q,n-j+1)*fy(:,j+1) ;
        end
        Phi_n = Phi_n + temp/Probl.lam_Q ;
        
        yn0 = y_pred ; fn0 = f_pred ;
        stop = 0 ; mu_it = 0 ;
        while ~stop
            yn1 = Phi_n - METH.C*yn0 + METH.an(Q,1)/Probl.lam_Q*fn0 ;
            mu_it = mu_it + 1 ;
            if METH.mu == Inf
                stop = norm(yn1-yn0,inf) < METH.mu_tol ;
                if mu_it > 100 && ~stop
                    warning('MATLAB:NonConvegence',...
                        strcat('It has been requested to run corrector iterations until convergence but ', ...
                        'the process does not converge to the tolerance %e in 100 iteration'),METH.mu_tol) ;
                    stop = 1 ;
                end
            else
                stop = mu_it == METH.mu ;
            end
            fn1 = f_vectorfield(t(n+1),yn1,Probl) ;             
            yn0 = yn1 ; fn0 = fn1 ;
        end
        y(:,n+1) = yn1 ;
        fy(:,n+1) = fn1 ;
    end
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

