using FdeSolver
using Plots, LaTeXStrings
using SpecialFunctions

## inputs
tSpan = [0, 300]   # time span
h = 0.05           # time step
N = 3         # number of species
# order of derivatives
μ1=[1;1;1] # No memory
μ2=.9*[1;1;1] # With memory

X0 = [.99;.01;.01]   # initial abundances

## System definition

# parametrisation
par1 = [2,
       [1;.95;1.05],
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
    if pulse==1
        bG=2
    elseif pulse==2
        bG=2.2
    end
    bB=0.5

    if t>20.0 && t<60.0
        B[1]=bB
        B[3]=bG
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
_, x2 = FDEsolver(F, tSpan, X0, μ2, par1, h = h)
_, x12 = FDEsolver(F, tSpan, X0, μ1, par2, h = h)
_, x22 = FDEsolver(F, tSpan, X0, μ2, par2, h = h)


#perturbations
function Fun_pert1(t)
    bb,br,bg=  [1,.95, 1.05]
    if t > 20 && t <60
            bg=2
            bb=0.5
    end
    return [bb,br,bg]
end
fbb(x)=20<x<60 ? 0.5 : 1
fbr(x)=1.05
fbg1(x)=20<x<60 ? 2 : 0.95
fbg2(x)=20<x<60 ? 2.2 : 0.95

## Plotting
myColor=[:snow2 :snow3 :blue :red :gray30]

pltP1=vspan([20, 60],color = myColor[1], labels="Low pulse perturbation")
    plot!(0:300,[fbb,fbr,fbg1],palette= [myColor[5],myColor[3],myColor[4]],xaxis = "Time",
     labels=[L"b_B" L"b_R" L"b_{G}"],
    title = "(a)" , titleloc = :left, titlefont = font(10))
    yaxis!("Growth rate of species", :log10, minorgrid = true)
pltP2=vspan([20, 60],color = myColor[2], labels= "Strong pulse perturbation")
        plot!(0:300,[fbb,fbr,fbg2],palette= [myColor[5],myColor[3],myColor[4]],xaxis = "Time",
        title = "(c)" , titleloc = :left, titlefont = font(10), labels=:false)
        yaxis!("Growth rate of species", :log10, minorgrid = true)


p1=vspan([20, 60],color = myColor[1], labels= :false)
          plot!(t, x1, palette = [myColor[5],myColor[3],myColor[4]],
          linestyle=:dash, labels=[L"X_B \;, \qquad 1" L"X_R \;, \qquad 1" L"X_{G}\;, \qquad 1"],
          legendtitle="Species, Orders",  legendposition=:right,legendtitlefont=font(10))

          yaxis!("Species abundance", :log10, minorgrid = true,
          title = "(b) " , titleloc = :left, titlefont = font(10))

        plot!(t, x2, palette=[myColor[5],myColor[3],myColor[4]],
        labels=[L"X_B\;, \qquad 0.9" L"X_R\;, \qquad 0.9" L"X_{G}\;, \qquad 0.9"])

p2=vspan([20, 60],color = myColor[2], labels= :false)
        plot!(t, x12, palette = [myColor[5],myColor[3],myColor[4]],
        linestyle=:dash,    xaxis = "Time",labels= :false)
        yaxis!("Species abundance", :log10, minorgrid = true,
        title = "(d) ", titleloc = :left, titlefont = font(10))

        plot!(t, x22, palette=[myColor[5],myColor[3],myColor[4]],labels= :false)

# l = @layout [b{.4w} grid(1,1);b{.4w} grid(1,1)]
P=plot(pltP1, p1, pltP2, p2, layout = [2,2], size=(820,600))


savefig(P, "myplot.svg")
