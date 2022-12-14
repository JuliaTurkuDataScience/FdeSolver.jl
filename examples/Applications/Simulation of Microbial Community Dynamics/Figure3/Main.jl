using FdeSolver
using Plots
using SpecialFunctions
using Statistics

## inputs
tSpan = [0, 450]   # time span
h = 0.01           # time step
N = 3         # number of species
# order of derivatives
μ1=.9*[1;1;1] # No memory

X0 = [.99 ; .01; .01]   # initial abundances

## System definition

# parametrisation
par1 = [2,
       [1;.95; 1.05],
       ones(N),
      .1*ones(N,N),
       N,
       1]
par2 = copy(par1) #deepcopy is another option
par2[6]= 2

function F(t, x, par)
    l = par[1] # Hill coefficient
    b = par[2] # growth rates
    k = par[3] # death rates
    K = par[4] # inhibition matrix
    N = par[5] # number of species
    pulse = par[6] # pulse type

    # pulse pertubation
    B=copy(b)
    mm=20; m=ceil(mod(t/(mm*4),mm))
    bB1=0.2; bB2=4.5

    if pulse==1

        if t>60.0 && t<100.0
            B[1]=bB1
        elseif  t>200 && t<330.0
            B[1]=bB2
        end

    elseif pulse==2

        if t>=mm*(4*m-3) && t<mm*(4*m-2)
            B[1]=bB1
        elseif t>=mm*(4*m-1) && t<mm*(4*m)
            B[1]=bB2
        end

    end

# ODE
    Fun = zeros(N)
    for i in 1:N
    # inhibition functions
    f = prod(K[i,1:end .!= i] .^ l ./
             (K[i,1:end .!= i] .^ l .+ x[ 1:end .!= i] .^l))
    # System of equations
    Fun[i] = x[ i] .* (B[i] .* f .- k[i] .* x[ i])
    end

    return Fun

end

t, x1 = FDEsolver(F, tSpan, X0, μ1, par1, h = h)

t, x2 = FDEsolver(F, tSpan, X0, μ1, par2, h = h)


X1=x1./sum(x1,dims=2)
X2=x2./sum(x2,dims=2)
##plotting
using ColorTypes
using CairoMakie

myColor=[:grey,:black,:blue, :red, :green]
p1=plot(t, X1, w = 5, palette = myColor[3:5],legend=:right,
            labels=["X_B/X_total" "X_R/X_total" "X_G/X_total"])
            p11=vspan(p1,[60, 100],color = :grey, alpha = 0.2, labels= :false)
            p111=vspan(p11,[200, 330],color = :black, alpha = 0.3, labels= :false)

savefig(p111, "fig3ploscb1.svg")
#background
mm=20; T2=4501
m=Int16(ceil(mod(T2/(mm*4),mm)))
v=zeros(m,2)
v1=zeros(m,2)
for i=1:m
        v[i,:]=[mm*(4*i-3) ; mm*(4*i-2)]
        v1[i,:]=[mm*(4*i-1) ; mm*(4*i)]
end

p2=vspan([v'], color = :grey, alpha = 0.3, labels= :false)
          p22=vspan(p22,[v1'],color = :black, alpha = 0.3, labels= :false)
          p222=plot!(t, X2, w = 5, palette = myColor[3:5], legend=:false
                    , labels=["X_B/X_total" "X_R/X_total" "X_G/X_total"])

savefig(p222, "fig3ploscb2.svg")
