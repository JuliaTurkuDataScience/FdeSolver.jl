using FdeSolver
using Plots

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
    f = prod(K[i,1:end .!= i] .^ l ./
             (K[i,1:end .!= i] .^ l .+ x[n, 1:end .!= i] .^l))
    # System of equations
    Fun[i] = x[n, i] .* (b[i] .* f .- k[i] .* x[n, i])
    end

    return Fun

end

# numerical solution
t, Xapp = FDEsolver(F, tSpan, X0, β, nothing, par, h = h, nc = 3, tol = 10^(-8))

# plot
plot(t, Xapp, linewidth = 5,
     title = "Dynamic of microbial intercation model",
     xaxis = "Time (t)")
     yaxis!("Log Abundance", :log10, minorgrid = true)
