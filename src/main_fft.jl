# If the 5th argument is integer then distpach (of FDEsolver) works for main_fft.That it not nice!
function FDEsolver(F, tSpan, y0, β, ::Int, ::Nothing, par...; h = 0.01, nc = 1, StopIt = "Standard", tol = 10e-6, itmax = 100)

# some renames
alpha = β
y0 = defineY0(y0, β)


# Check compatibility size of the problem with number of fractional orders
alpha_length = length(alpha)
problem_size = size(y0,1)

if alpha_length > 1
        # println("Oh! Good! ODE system!")
else
alpha = alpha*ones(problem_size,1)
alpha_length = problem_size
end

# Storage of initial conditions
initial_conditions = @SLVector (:t0 ,:y0,:m_alpha, :m_alpha_factorial)
ic = initial_conditions(tSpan[1], y0, Int64.(map(ceil, alpha)),
     zeros(alpha_length,Int64.(ceil(maximum(alpha)))))
for i in 1 : alpha_length
    for j in 0 : ic.m_alpha[i]-1
        ic.m_alpha_factorial[i,j+1] = factorial(j) ;
    end
end

# Storage of information on the problem
Problem = @SLVector (:ic , :f_fun, :problem_size, :param, :alpha, :alpha_length)
Probl = Problem(ic , F, problem_size,par, alpha, alpha_length)

# Time discretization
N = Int64.(cld(tSpan[2] - tSpan[1], h))
t = tSpan[1] .+ collect(0:N) .* h

# Enter the initial values
# Y = defineY(N, y0, β)

# Check compatibility size of the problem with size of the vector field
f_temp = F(t, 1, β, y0, par...)

# Number of points in which to evaluate weights and solution
r = 16
Nr::Int64 = ceil((N+1)/r)*r
Qr::Int64 = ceil(log2(Nr/r)) - 1
NNr::Int64 = 2^(Qr+1)*r

# Preallocation of some variables
y = zeros(Probl.problem_size, N+1)
fy = zeros(Probl.problem_size, N+1)
zn_pred = zeros(Probl.problem_size,NNr+1)
if nc > 0
    zn_corr = zeros(Probl.problem_size,NNr+1)
else
    zn_corr = 0
end


# Evaluation of coefficients of the PECE method
nvett = 0 : NNr+1 ;
bn = zeros(Probl.alpha_length,NNr+1)
an = zeros(Probl.alpha_length,NNr+1)
a0 = zeros(Probl.alpha_length,NNr+1)

for i_alpha in 1 : Probl.alpha_length
    find_alpha = findall(alpha[i_alpha]==alpha[1:i_alpha-1])
    if ~isempty(find_alpha) # it is for speeding up the computations; we can use multilpe distpach
        bn[i_alpha,:] = bn[find_alpha[1],:]
        an[i_alpha,:] = an[find_alpha[1],:]
        a0[i_alpha,:] = a0[find_alpha[1],:]
    else
        nalpha = nvett.^alpha[i_alpha]
        nalpha1 = nalpha.*nvett
        bn[i_alpha,:] = nalpha[2:end] - nalpha[1:end-1]
        an[i_alpha,:] = [ 1 ; (nalpha1[1:end-2] - 2*nalpha1[2:end-1] + nalpha1[3:end]) ]
        a0[i_alpha,:] = [ 0 ; (nalpha1[1:end-2] - nalpha[2:end-1].*(nvett[2:end-1] .- alpha[i_alpha] .- 1))]
    end
end
Method = @SLVector (:bn, :an, :a0, :hα1, :hα2, :μ, :μTol, :r, :StopIt,:itmax)
METH = Method(bn, an, a0, h.^alpha./Γ(alpha.+1), h.^alpha./Γ(alpha.+2), nc, tol, r, StopIt, itmax)

# Evaluation of FFT of coefficients of the PECE method
if Qr >= 0
    index_fft = Int64.(zeros(2,Qr+1)) # I have tried index_fft::Int64 = zeros(2,Qr+1) and I got an error for converting Type!
    for l in 1 : Qr+1 # log2(NNr/r)
        if l == 1
            index_fft[1,l] = 1 ; index_fft[2,l] = r*2
        else
            index_fft[1,l] = index_fft[2,l-1] + 1
            index_fft[2,l] = index_fft[2,l-1]+2^l*r
        end
    end

    bn_fft =ComplexF64.(zeros(Probl.alpha_length,index_fft[2,Qr+1]))
    an_fft =ComplexF64.(zeros(Probl.alpha_length,index_fft[2,Qr+1]))
    for l in 1:Qr+1 # log2(NNr/r)
        coef_end = 2^l*r
        for i_alpha = 1 : Probl.alpha_length
            find_alpha = findall(alpha[i_alpha]==alpha[1:i_alpha-1])
            if ~isempty(find_alpha)
                bn_fft[i_alpha,index_fft[1,l]:index_fft[2,l]] = bn_fft[find_alpha[1],index_fft[1,l]:index_fft[2,l]]
                an_fft[i_alpha,index_fft[1,l]:index_fft[2,l]] = an_fft[find_alpha[1],index_fft[1,l]:index_fft[2,l]]
            else
                bn_fft[i_alpha,index_fft[1,l]:index_fft[2,l]] = fft(METH.bn[i_alpha,1:coef_end])
                an_fft[i_alpha,index_fft[1,l]:index_fft[2,l]] = fft(METH.an[i_alpha,1:coef_end])
            end
        end
    end
    Method_fft = @SLVector (:bn_fft,:an_fft,:index_fft)
    METH_fft = Method_fft(bn_fft, an_fft, Int64.(index_fft))
end

# Initializing solution and proces of computation
y[:,1] = y0[:,1]
fy[:,1] .= f_temp
y, fy = Triangolo(1, r-1, t, y, fy, zn_pred, zn_corr, N, METH, METH_fft, Probl) ;

# Main process of computation by means of the FFT algorithm
ff = zeros(1,2 .^ (Qr+2)) ; ff[1:2] = [0, 2] ; card_ff = 2
nx0 = 0 ; ny0 = 0
for qr = 0 : Qr
    L = 2^qr
    y, fy = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy,
        zn_pred, zn_corr, N, METH, METH_fft, Probl)
    ff[1:2*card_ff] = [ff[1:card_ff]; ff[1:card_ff]]
    card_ff = 2*card_ff
    ff[card_ff] = 4*L
end

 # Evaluation solution in T when T is not in the mesh
if tSpan[end] < t[N+1]
    c = [tSpan[end] - t(N)]/h
    t[N+1] = tSpan[end]
    y[:,N+1] = (1-c)*y[:,N] + c*y[:,N+1]
end

t = t[1:N+1] ; y = y[:,1:N+1]

return t, y
end
