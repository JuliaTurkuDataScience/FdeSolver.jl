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
F(t, y, par)= par * y
JF(t, y, par) = y

## Numerical solution
t,y = FDEsolver(F, nothing, tSpan, y0, β, par)
t1,y1 = FDEsolver(F, JF, tSpan, y0, β, par)
plot(t,y)
plot!(t1,y1, ls = :dash)
# Juno.@enter FDEsolver(F, tSpan, y0, β, 1,  nothing, par)

############ Let's try FFT for a system
tSpan = [0, 50]
h = 2^(-5)
y0 = [ 1.2 , 2.8]
β = [0.8,0.7]
A = 1 ; B = 3
par = [ A , B ]
# Definition of the System
function F(t, y, par)

    F1 = par[1] .- (par[2]+1) .* y[1] .+ y[1] .^ 2 .* y[2]
    F2 = par[2] .* y[1] .- y[1] .^ 2 .* y[2]

    return [F1, F2]

end

@time t, Yapp = FDEsolver(F, nothing, tSpan, y0, β, par)
plot(t, Yapp)

#######################
# Benchmark with MATLAB
## inputs
tSpan = [0, 700] # time span

h = 0.01 # time step

N = 15 # number of species

β = ones(N) # order of derivatives

X0 = 1/15 .*ones(N) # initial abundances

## System definition

# parametrisation
par = [2,
      ones(N),
      ones(N),
      4*rand(N,N),
       N]

function F(t, x, par)
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
             (K[1:end .!= i, i] .^ l .+ x[1:end .!= i] .^l))
    # System of equations
    Fun[i] = x[i] .* (b[i] .* f .- k[i] .* x[i])
    end

    return Fun

end

using BenchmarkTools
@benchmark FDEsolver(F, nothing, tSpan, X0, β, par, nc =1)
