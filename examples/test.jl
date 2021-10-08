exit()
push!(LOAD_PATH, "src/")
using FdeSolver
using Plots

## inputs
tSpan = [0, 5]     # [intial time, final time]
y0 = 1             # intial value
β = 0.6            # order of the derivative

alpha = 0.6
lambda = -10
par = lambda
F(t, n, β, y,par)= par * y[n]

## Numerical solution
@time t, y = FDEsolver(F, tSpan, y0, β, 1, nothing, par)
@time t, y = FDEsolver(F, tSpan, y0, β, nothing, nothing, par)
@time t, y = FDEsolver(F, tSpan, y0, β, 1, nothing, par, StopIt = "Convergence", tol = 10e-10, itmax = 100)
# If the 5th argument is integer then distpach (of FDEsolver) works for main_fft.That it not nice!
plot(t,y[1,:]) # This is not nice that we have to say y[1,:] instead of y, in a 1D example.
# Juno.@enter FDEsolver(F, tSpan, y0, β, 1,  nothing, par)

############ Let's try FFT for a system
tSpan = [0, 50]
h = 2^(-5)
y0 = [ 1.2 , 2.8]
β = [0.8,0.7]
A = 1 ; B = 3
par = [ A , B ]
# Definition of the System
function F(t, n, β, y, par)

    F1 = par[1] .- (par[2]+1) .* y[1,n] .+ y[1,n] .^ 2 .* y[2,n]
    F2 = par[2] .* y[1,n] .- y[1,n] .^ 2 .* y[2,n]

    return [F1, F2]

end

@time t, Yapp = FDEsolver(F, tSpan, y0, β, 1, nothing, par)
# Juno.@enter FDEsolver(F, tSpan, y0, β, 1,  nothing, par)

#######################

## inputs
tSpan = [0, 50] # time span

h = 0.1 # time step

N = 10 # number of species

β = ones(N) # order of derivatives

X0 = 2*rand(N) # initial abundances

## System definition

# parametrisation
par = [2,
      2*rand(N),
      rand(N),
      4*rand(N,N),
       N]

function F(t, n, β, x, par)
    l = par[1] # Hill coefficient
    b = par[2] # growth rates
    k = par[3] # death rates
    K = par[4] # inhibition matrix
    N = par[5] # number of species

# ODE
    Fun = zeros(N)
    for i in 1:N
    # inhibition functions
    f = prod(K[1:end .!= i, i] .^ l ./
             (K[1:end .!= i, i] .^ l .+ x[1:end .!= i, n] .^l))
    # System of equations
    Fun[i] = x[i, n] .* (b[i] .* f .- k[i] .* x[i, n])
    end

    return Fun

end

t, Yapp = FDEsolver(F, tSpan, X0, β, 1, nothing, par, nc =1)
plot(t, transpose(Yapp))
