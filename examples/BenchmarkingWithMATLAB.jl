push!(LOAD_PATH, "src/")
using FdeSolver
using Plots
using BenchmarkTools
using MittagLeffler#for the exact solution
using TimerOutputs
using CSV, DataFrames
using LinearAlgebra
## inputs
H=[2^(-2), 2^(-3), 2^(-4), 2^(-5),  2^(-6),  2^(-7), 2^(-8)]
tSpan = [0, 5]     # [intial time, final time]
y0 = 1             # intial value
β = 0.6            # order of the derivative

λ = -10
par = λ
F(t, y, par)= par * y
JF(t, y, par) = par
h=zeros()
Bench1=(zeros(length(H),2,2))

for i in 1:length(H)
        h = H[i]
    Bench1[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, h=h)).time*10^-9
    Bench1[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, J = JF, h=h)).time*10^-9

    if i>=4
    Exact(t)=mittleff(β,λ *t .^ β)
    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h)
    Bench1[i,1,2]=norm(y-map(Exact,t),2)
    Bench1[i,2,2]=norm(y1-map(Exact,t),2)
    end

end
Mdata=DataFrame(CSV.File("Bench1.csv",header=0));
plot(H,Bench1[:,:,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example1",yaxis="Time(Sc)",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H,Mdata[!,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H[4:end],Bench1[4:end,:,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H[4:end],Mdata[4:end,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H[4:end],Mdata[4:end,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")

##
using SpecialFunctions
H=[2^(-2), 2^(-3), 2^(-4), 2^(-5),  2^(-6),  2^(-7), 2^(-8)]
tSpan = [0, 1]     # [intial time, final time]
y0 = 0             # intial value
β = 0.5            # order of the derivative

# Equation
par = β
F(t, y, par) = (40320 ./ gamma(9 - par) .* t .^ (8 - par) .- 3 .* gamma(5 + par / 2)
           ./ gamma(5 - par / 2) .* t .^ (4 - par / 2) .+ 9/4 * gamma(par + 1) .+
           (3 / 2 .* t .^ (par / 2) .- t .^ 4) .^ 3 .- y .^ (3 / 2))
# Jacobian
JF(t, y, par) = -(3 / 2) .* y .^ (1 / 2)

h=zeros()
Bench2=(zeros(length(H),2,2))

for i in 1:length(H)
        h = H[i]
    Bench2[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, h=h)).time*10^-9
    Bench2[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, J = JF, h=h)).time*10^-9

    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h)
    Exact=t.^8 - 3 * t .^ (4 + β / 2) + 9 / 4 * t.^β
    Bench2[i,1,2]=norm(y-Exact)
    Bench2[i,2,2]=norm(y1-Exact)


end
Mdata=DataFrame(CSV.File("Bench2.csv",header=0));
plot(H,Bench2[:,:,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example2",yaxis="Time(Sc)",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H,Mdata[!,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H,Bench2[:,:,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H,Mdata[!,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")
