using FdeSolver
using Plots

## inputs
tSpan = [0, 5000] # time span

H=[2^(-5), 2^(-6),  2^(-7), 2^(-8)]

N = 20 # number of species

β = ones(N) # order of derivatives

X0 = 2*rand(N) # initial abundances

## System definition

# parametrisation
par = [2,
      2*rand(N),
      rand(N),
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
    f = prod(K[i,1:end .!= i] .^ l ./
             (K[i,1:end .!= i] .^ l .+ x[ 1:end .!= i] .^l))
    # System of equations
    Fun[i] = x[ i] .* (b[i] .* f .- k[i] .* x[ i])
    end

    return Fun

end

h=zeros()
Bench=zeros(length(H))

for i in 1:length(H)
        h = H[i]
    Bench[i]= mean(@benchmark FDEsolver(F, tSpan, X0, β, par, h=h)).time*10^-9
end

Mdata=DataFrame(CSV.File("BenchLongTerm.csv",header=0));
plot(H[:],Bench[:],linewidth=5,title="Becnhmark longTerm",yaxis="Time(Sc)",xaxis="Step size", label="Julia")
plot!(H[:],Mdata[!,1],linewidth=5, ls = :dot,label="MATLAB")
