function [t, y] = fde_pi12_pcIM(alpha,f_fun,J_fun,t0,T,y0,h,param,mu,mu_tol)
% FDE_PI12_PC Solves on the interval [t0,T] an initial value problem for a 
%             for a system D^alpha y = f(t,y) of differential equations of
%             fractional order alpha. Multi-order systems (MOS) are also
%             possible, i.e. systems in which each equation has a different
%             fractional order. 
%             The system or MOS of FDEs is solved by means of a couple of 
%             product-integration rules of rectangular and trapezoidal type 
%             implemented in a predictor-corrector framework.
%
% Further information on this code are available in [1]; please, cite this
% code as [1] (see later on). 
%
% -------------------------------------------------------------------------
% Description of input parameters
% -------------------------------------------------------------------------
% alpha     : vector of the fractional orders alpha of each equation of the
%             multi-order system of FDEs. 
% f_fun     : function handle that evaluates the right side f(t,y(t)) of
%             the MOS of FDEs. f_fun(t,y) must return a vector function  
%             with the same number of entries as alpha. If necessary,  
%             f_fun(t,y,param) can have a vector of parameters. 
% t0, T     : starting and final time of the interval of integration.
% y0        : vector (or matrix) of initial conditions with a number of 
%             rows equal to the size of the problem (hence to the number of
%             entries of alpha or to the number of rows of the output of
%             f_fun) and a number of columns equal to the smallest integer
%             greater than max(alpha)
% h         : step-size for integration. It must be real and positive
% param     : vector of possible parameters for the evaluation of the
%             right side f(t,y(t)) and of its Jacobian. This an optional 
%             parameter: it can be omitted or it is possible to use the an
%             empty vector [] if the function f_fun(t,y) and its Jacobian 
%             do not have parameters.
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
% [t, y] = FDE_PI2_PC(alpha,f_fun,t0,T,y0,h)
% [t, y] = FDE_PI2_PC(alpha,f_fun,t0,T,y0,h,param)
% [t, y] = FDE_PI2_PC(alpha,f_fun,t0,T,y0,h,param,mu)
% [t, y] = FDE_PI2_PC(alpha,f_fun,t0,T,y0,h,param,mu,mu_tol)
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
%  Version of November, 30 2017
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
        mu = 2 ;
        if nargin < 8
            param = [] ;
        end
    end
end

% Check order of the FDE
alpha = alpha(:) ;
if any(alpha <= 0)
    i_alpha_neg = find(alpha<=0) ;
    error('MATLAB:NegativeOrder',...
        ['The orders ALPHA of the FDEs must be all positive. The value ' ...
        'ALPHA(%d) = %f can not be accepted.'], i_alpha_neg(1), alpha(i_alpha_neg(1)));
end

% Check the step--size of the method
if h <= 0
    error('MATLAB:NegativeStepSize',...
        ['The step-size H for the method must be positive. The value ' ...
        'H = %e can not be accepted.'], h);
end

% Check compatibility size of the problem with number of fractional orders
alpha_length = length(alpha) ; 
problem_size = size(y0,1) ;
if alpha_length > 1
    if problem_size ~= alpha_length
        error('MATLAB:SizeNotCompatible', ...
            ['Size %d of the problem as obtained from initial conditions ' ...
            '(i.e. the number of rows of Y0) not compatible with the ' ...
            'number %d of fractional orders for multi-order systems. ' ...
            ], problem_size, alpha_length);
    end
else
    alpha = alpha*ones(problem_size,1) ; 
    alpha_length = problem_size ;
end

% Storage of initial conditions
ic.t0 = t0 ;
ic.y0 = y0 ;
ic.m_alpha = ceil(alpha) ;
for i = 1 : alpha_length
    for j = 0 : ic.m_alpha(i)-1
        ic.m_alpha_factorial(i,j+1) = factorial(j) ;
    end
end

% Storage of information on the problem
Probl.ic = ic ;
Probl.f_fun = f_fun ;
Probl.J_fun = J_fun ;
Probl.problem_size = problem_size ;
Probl.param = param ;
Probl.alpha = alpha ;
Probl.alpha_length = alpha_length ;

% Check number of initial conditions
if size(y0,2) < max(ic.m_alpha)
    error('MATLAB:NotEnoughInitialConditions', ...
        ['A not sufficient number of assigned initial conditions. ' ...
        'Order ALPHA = %f requires %d initial conditions.'], ...
        max(alpha),max(ic.m_alpha));
end

% Check compatibility size of the problem with size of the vector field
f_temp = f_vectorfield(t0,y0(:,1),Probl) ;
if Probl.problem_size ~= size(f_temp,1)
    error('MATLAB:SizeNotCompatible', ...
        ['Size %d of the problem as obtained from initial conditions ' ...
        '(i.e. the number of rows of Y0) not compatible with the ' ...
        'size %d of the output of the vector field F_FUN. ' ...
        ], Probl.problem_size,size(f_temp,1));
end

% Check compatibility size of the problem with number of fractional order
if Probl.alpha_length > 1
    if Probl.problem_size ~= Probl.alpha_length
        error('MATLAB:SizeNotCompatible', ...
            ['Size %d of the problem as obtained from initial conditions ' ...
            '(i.e. the number of rows of Y0) not compatible with the ' ...
            'number %d of fractional order for noncommensurate systems. ' ...
            ], Probl.problem_size,Probl.alpha_length);
    end
end

% Number of points in which to evaluate weights and solution
r = 16 ;
N = ceil((T-t0)/h) ;
Nr = ceil((N+1)/r)*r ;
Qr = ceil(log2(Nr/r)) - 1 ;
NNr = 2^(Qr+1)*r ;

% Preallocation of some variables
y = zeros(Probl.problem_size,N+1) ;
fy = zeros(Probl.problem_size,N+1) ;
zn_pred = zeros(Probl.problem_size,NNr+1) ;
if mu > 0
    zn_corr = zeros(Probl.problem_size,NNr+1) ;
else
    zn_corr = 0 ;
end

% Evaluation of coefficients of the PECE method
nvett = 0 : NNr+1 ; 
bn = zeros(Probl.alpha_length,NNr+1) ; an = bn ; a0 = bn ;
for i_alpha = 1 : Probl.alpha_length
    find_alpha = find(alpha(i_alpha)==alpha(1:i_alpha-1)) ;
    if ~isempty(find_alpha)
        bn(i_alpha,:) = bn(find_alpha(1),:) ;
        an(i_alpha,:) = an(find_alpha(1),:) ;
        a0(i_alpha,:) = a0(find_alpha(1),:) ;
    else
        nalpha = nvett.^alpha(i_alpha) ;
        nalpha1 = nalpha.*nvett ;
        bn(i_alpha,:) = nalpha(2:end) - nalpha(1:end-1) ;
        an(i_alpha,:) = [ 1 , (nalpha1(1:end-2) - 2*nalpha1(2:end-1) + nalpha1(3:end)) ] ;
        a0(i_alpha,:) = [ 0 , nalpha1(1:end-2)-nalpha(2:end-1).*(nvett(2:end-1)-alpha(i_alpha)-1)] ;
    end
end
METH.bn = bn ; METH.an = an ; METH.a0 = a0 ;
METH.halpha1 = h.^alpha./gamma(alpha+1) ;
METH.halpha2 = h.^alpha./gamma(alpha+2) ;
METH.mu = mu ; METH.mu_tol = mu_tol ;

% Evaluation of FFT of coefficients of the PECE method
METH.r = r ;
if Qr >= 0
    index_fft = zeros(2,Qr+1) ;
    for l = 1 : Qr+1 % log2(NNr/r)
        if l == 1
            index_fft(1,l) = 1 ; index_fft(2,l) = r*2 ;
        else
            index_fft(1,l) = index_fft(2,l-1)+1 ; index_fft(2,l) = index_fft(2,l-1)+2^l*r;
        end
    end
    
    bn_fft = zeros(Probl.alpha_length,index_fft(2,Qr+1)) ; an_fft = bn_fft ;
    for l = 1 : Qr+1 % log2(NNr/r)
        coef_end = 2^l*r ;
        for i_alpha = 1 : Probl.alpha_length
            find_alpha = find(alpha(i_alpha)==alpha(1:i_alpha-1)) ;
            if ~isempty(find_alpha)
                bn_fft(i_alpha,index_fft(1,l):index_fft(2,l)) = bn_fft(find_alpha(1),index_fft(1,l):index_fft(2,l)) ;
                an_fft(i_alpha,index_fft(1,l):index_fft(2,l)) = an_fft(find_alpha(1),index_fft(1,l):index_fft(2,l)) ;
            else
                bn_fft(i_alpha,index_fft(1,l):index_fft(2,l)) = fft(METH.bn(i_alpha,1:coef_end),coef_end) ;
                an_fft(i_alpha,index_fft(1,l):index_fft(2,l)) = fft(METH.an(i_alpha,1:coef_end),coef_end) ;
            end
        end
    end
    METH.bn_fft = bn_fft ; METH.an_fft = an_fft ; METH.index_fft = index_fft ;
end

% Initializing solution and proces of computation
t = t0 + (0 : N)*h ;
y(:,1) = y0(:,1) ;
fy(:,1) = f_temp ;
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

% Evaluation solution in T when T is not in the mesh
if T < t(N+1)
    c = (T - t(N))/h ;
    t(N+1) = T ;
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
    
    stop = nxi+r-1 == nx0+L*r-1 | (nxi+r-1>=Nr-1) ; % Ci si ferma quando il triangolo attuale finisce alla fine del quadrato
    
    [zn_pred, zn_corr] = Quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, N, METH, Probl) ;
    
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
function [zn_pred, zn_corr] = Quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, N, METH, Probl)

%coef_beg = nxi-nyf ; 
coef_end = nxf-nyi+1 ;
i_fft = log2(coef_end/METH.r) ;
funz_beg = nyi+1 ; funz_end = nyf+1 ;
Nnxf = min(N,nxf) ; 

% Evaluation convolution segment for the predictor
vett_funz = fy(:,funz_beg:funz_end) ;
vett_funz_fft = fft(vett_funz,coef_end,2) ;
zzn_pred = zeros(Probl.problem_size,coef_end) ;
for i = 1 : Probl.problem_size
    i_alpha = min(Probl.alpha_length,i) ;
    Z = METH.bn_fft(i_alpha,METH.index_fft(1,i_fft):METH.index_fft(2,i_fft)).*vett_funz_fft(i,:) ;
    zzn_pred(i,:) = real(ifft(Z,coef_end)) ;
end
zzn_pred = zzn_pred(:,nxf-nyf+1-1:end-1) ;
zn_pred(:,nxi+1:Nnxf+1) = zn_pred(:,nxi+1:Nnxf+1) + zzn_pred(:,1:Nnxf-nxi+1) ;

% Evaluation convolution segment for the corrector
if METH.mu > 0
    if nyi == 0 % Evaluation of the lowest square
        vett_funz = [zeros(Probl.problem_size,1) , fy(:,funz_beg+1:funz_end) ] ;
        vett_funz_fft = fft(vett_funz,coef_end,2) ;
    end
    zzn_corr = zeros(Probl.problem_size,coef_end) ;
    for i = 1 : Probl.problem_size
        i_alpha = min(Probl.alpha_length,i) ;
        Z = METH.an_fft(i_alpha,METH.index_fft(1,i_fft):METH.index_fft(2,i_fft)).*vett_funz_fft(i,:) ;
        zzn_corr(i,:) = real(ifft(Z,coef_end)) ;
    end
    zzn_corr = zzn_corr(:,nxf-nyf+1:end) ;
    zn_corr(:,nxi+1:Nnxf+1) = zn_corr(:,nxi+1:Nnxf+1) + zzn_corr(:,1:Nnxf-nxi+1) ;
else
    zn_corr = 0 ;
end

end



% =========================================================================
% =========================================================================
function [y, fy] = Triangolo(nxi, nxf, t, y, fy, zn_pred, zn_corr, N, METH, Probl)

for n = nxi : min(N,nxf)
    
    % Evaluation of the predictor
    Phi = zeros(Probl.problem_size,1) ;
    if nxi == 1 % Case of the first triangle
        j_beg = 0 ;
    else % Case of any triangle but not the first
        j_beg = nxi ;
    end
    for j = j_beg : n-1
        Phi = Phi + METH.bn(1:Probl.alpha_length,n-j).*fy(:,j+1) ;
    end
    St = StartingTerm(t(n+1),Probl.ic) ;
    Phi_n = St + METH.halpha1.*(zn_pred(:,n+1) + Phi) ;
        
    yn0 = y(:,n) ; fn0 = f_vectorfield(t(n+1),yn0,Probl) ;
    Jfn0 = Jf_vectorfield(t(n+1),yn0,Probl) ;
    Gn0 = yn0 - METH.halpha1.*fn0 - Phi_n ;
    JGn0 = eye(Probl.problem_size,Probl.problem_size) - diag(METH.halpha1)*Jfn0 ;
    y_pred = yn0 - JGn0\Gn0 ;
    f_pred = f_vectorfield(t(n+1),y_pred,Probl) ;
    
    % Evaluation of the corrector
    if METH.mu == 0
        y(:,n+1) = y_pred ;
        fy(:,n+1) = f_pred ;
    else
        j_beg = nxi ;
        Phi = zeros(Probl.problem_size,1) ;
        for j = j_beg : n-1
            Phi = Phi + METH.an(1:Probl.alpha_length,n-j+1).*fy(:,j+1) ;
        end
        Phi_n = St + ...
            METH.halpha2.*(METH.a0(1:Probl.alpha_length,n+1).*fy(:,1) + zn_corr(:,n+1) + Phi) ;
        yn0 = y_pred ; fn0 = f_pred ;   
        Jfn0 = Jf_vectorfield(t(n+1),yn0,Probl) ;
        Gn0 = yn0 - METH.halpha2.*fn0 - Phi_n ;
        
        stop = 0 ; mu_it = 0 ;
        while ~stop
            
        JGn0 = eye(Probl.problem_size,Probl.problem_size) - diag(METH.halpha2)*Jfn0 ;
        yn1 = yn0 - JGn0\Gn0 ;
        fn1 = f_vectorfield(t(n+1),yn1,Probl) ;
        Gn1 = yn1 - METH.halpha2.*fn1 - Phi_n ;
            
            mu_it = mu_it + 1 ;
            if METH.mu == Inf
                stop = norm(yn1-yn0,inf) < METH.mu_tol | norm(Gn1,inf) <  METH.mu_tol;
                if mu_it > 100 && ~stop
                    warning('MATLAB:fde12:NonConvegence',...
                        strcat('It has been requested to run corrector iterations until convergence but ', ...
                        'the process does not converge to the tolerance %e in 100 iteration'),METH.mu_tol) ;
                    stop = 1 ;
                end
            else
                stop = mu_it == METH.mu ;
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

end


% =========================================================================
% =========================================================================
function f = f_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.f_fun,t,y) ;
else
    f = feval(Probl.f_fun,t,y,Probl.param) ;
end

end

% =========================================================================
% =========================================================================
function f = Jf_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.J_fun,t,y) ;
else
    f = feval(Probl.J_fun,t,y,Probl.param) ;
end

end

% =========================================================================
% =========================================================================
function ys = StartingTerm(t,ic)

ys = zeros(size(ic.y0,1),1) ;
for k = 1 : max(ic.m_alpha)
    if length(ic.m_alpha) == 1
        ys = ys + (t-ic.t0)^(k-1)/ic.m_alpha_factorial(k)*ic.y0(:,k) ;
    else
        i_alpha = find(k<=ic.m_alpha) ;
        ys(i_alpha,1) = ys(i_alpha,1) + (t-ic.t0)^(k-1)*ic.y0(i_alpha,k)./ic.m_alpha_factorial(i_alpha,k) ;
    end
end

end