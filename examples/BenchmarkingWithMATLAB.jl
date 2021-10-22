push!(LOAD_PATH, "src/")
using FdeSolver
using Plots
using BenchmarkTools
using MittagLeffler#for the exact solution
using TimerOutputs
using CSV, DataFrames
using LinearAlgebra
using SpecialFunctions
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
Exact1(t)=mittleff(β,λ *t .^ β)
for i in 1:length(H)
        h = H[i]
    Bench1[i,1,1]= mean(@benchmark FDEsolver($F, $tSpan, $y0, $β, $par, h=$h)).time*10^-9
    Bench1[i,2,1]= mean(@benchmark FDEsolver($F, $tSpan, $y0, $β, $par, J = $JF, h=$h)).time*10^-9

    if i>=4

    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h)
    Bench1[i,1,2]=norm(y-map(Exact1,t),2)
    Bench1[i,2,2]=norm(y1-map(Exact1,t),2)
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


Exact2(t)=t.^8 - 3 * t .^ (4 + β / 2) + 9 / 4 * t.^β

for i in 2:length(H)
        h = H[i]
    Bench2[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, h=h)).time*10^-9
    Bench2[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, J = JF, h=h)).time*10^-9

    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h)
    Exact=map(Exact2,t)
    Bench2[i,1,2]=norm(y-Exact)
    Bench2[i,2,2]=norm(y1-Exact)


end
Mdata=DataFrame(CSV.File("Bench2.csv",header=0));
plot(H[2:end],Bench2[2:end,:,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example2",yaxis="Time(Sc)",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H[2:end],Mdata[2:end,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H[2:end],Mdata[2:end,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H[2:end],Bench2[2:end,:,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H[2:end],Mdata[2:end,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H[2:end],Mdata[2:end,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")


##
H=[2^(-2), 2^(-3), 2^(-4), 2^(-5),  2^(-6),  2^(-7), 2^(-8)]

tSpan = [0, 2]     # [intial time, final time]
y0 = 0             # intial value
β = 0.5            # order of the derivative
par =β
function F(t,y,β)
    if t > 1
        dy= 1/gamma(2-β)* t ^(1-β) - 2/gamma(3-β) *(t-1) ^(2-β)
    else
        dy= 1/gamma(2-β)*t^(1-β)
    end
    return dy
end
# Jacobian
JF(t, y, β) = 0

h=zeros()
Bench3=(zeros(length(H),2,2))

function Exact3(t)
    if t > 1
        y= t .-(t .-1).^2
        else
        y= t
    end
    return y
end


for i in 1:length(H)
        h = H[i]
    Bench3[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, h=h,)).time*10^-9
    Bench3[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, J = JF, h=h)).time*10^-9

    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h)

    Exact = map(Exact3,t)
    Bench3[i,1,2]=norm(y-Exact)
    Bench3[i,2,2]=norm(y1-Exact)


end
Mdata=DataFrame(CSV.File("Bench3.csv",header=0));
plot(H,Bench3[:,:,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example3",yaxis="Time(Sc)",xaxis="Step size", label=["J_PI_PC" "J_PI_IM"])
plot!(H,Mdata[!,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H,Bench3[:,:,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H,Mdata[!,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")


##
H=[2^(-2),2^(-3), 2^(-4), 2^(-5),  2^(-6),  2^(-7), 2^(-8)]

# Parameters
tSpan = [0, 5]         # Time Span
y0 = [1, 0.5, 0.3]     # Initial values
β = [0.5, 0.2, 0.6]    # Order of derivation

# Definition of the System
function F(t, y)

    F1 = 1 / sqrt(pi) * (((y[2] - 0.5) * (y[3] - 0.3))^(1 / 6) + sqrt(t))
    F2 = gamma(2.2) * (y[1] - 1)
    F3 = gamma(2.8) / gamma(2.2) * (y[2] - 0.5)

    return [F1, F2, F3]

end
#
# function JF(t,y)
#     # System equation
#     J11 = 0
#     J12 = (y[2]-0.5).^(-5/6).*(y[3]-0.3).^(1/6)/6/sqrt(pi)
#     J13 = (y[2]-0.5).^(1/6).*(y[3]-0.3).^(-5/6)/6/sqrt(pi)
#     J21 =  gamma(2.2)
#     J22 =  0
#     J23 =  0
#     J31 =  0
#     J32 =  gamma(2.8)/gamma(2.2)
#     J33 =  0
#
#     J = [J11 J12 J13
#          J21 J22 J23
#          J31 J32 J33]
#     return J
# end

h=zeros()
Bench4=(zeros(length(H),2,2))
Exact4(t)=[t .+ 1; t.^1.2 .+ 0.5; t.^1.8 .+ 0.3]

for i in 1:length(H)
        h = H[i]
    Bench4[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, h=h)).time*10^-9
    # Bench4[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, J = JF, h=h)).time*10^-9
    # Using Jacobian is not a good idea for this example, we face negative value (or zero) power negative value.

    t,y = FDEsolver(F, tSpan, y0, β, h=h)
    # t1,y1 = FDEsolver(F, tSpan, y0, β, J=JF, h=h)

    Exact = reshape(Exact4(t),length(t),3)
    Bench4[i,1,2]=norm(y-Exact)
    # Bench4[i,2,2]=norm(y1-Exact)


end
Mdata=DataFrame(CSV.File("Bench4.csv",header=0));
plot(H,Bench4[:,1,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example4",yaxis="Time(Sc)",xaxis="Step size", label="PI_PC")
plot!(H,Mdata[!,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H,Bench4[:,1,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label="PI_PC")
plot!(H,Mdata[!,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")

##

H=[2^(-2),2^(-3), 2^(-4), 2^(-5),  2^(-6),  2^(-7), 2^(-8)]

# Parameters
tSpan = [0, 1]         # Time Span
y0 = [0; 0]     # Initial values
β = 1.5    # Order of derivation

# Definition of the System
F(t, y , β) = t.^(β)*y.+ 4*sqrt(t/π).-t.^(2+β)

JF(t,y,β) = t.^(β)
par=β

h=zeros()
Bench5=(zeros(length(H),2,2))
Exact5(t)=t.^2

for i in 1:length(H)
        h = H[i]
    Bench5[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, h=h)).time*10^-9
    Bench5[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, J = JF, h=h)).time*10^-9

    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h)

    Exact = map(Exact5,t)
    Bench5[i,1,2]=norm(y-Exact)
    Bench5[i,2,2]=norm(y1-Exact)


end
Mdata=DataFrame(CSV.File("Bench5.csv",header=0));
plot(H,Bench5[:,:,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example5",yaxis="Time(Sc)",xaxis="Step size", label=["J_PI_PC" "J_PI_IM"])
plot!(H,Mdata[!,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H,Bench5[:,:,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H,Mdata[!,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")


##

H=[2^(-2),2^(-3), 2^(-4), 2^(-5),  2^(-6),  2^(-7), 2^(-8)]

tSpan = [0, 10]     # [intial time, final time]
y0 = [1; 1]             # intial value ([of order 0; of order 1])
β = 2            # order of the derivative
par = [16.0, 4.0] # [spring constant for a mass on a spring, inertial mass]
h = 0.01

function F(t, x, par)

      K = par[1]
      m = par[2]

      - K ./ m .* x

end


JF(t,y,par) = - par[1] ./ par[2]

h=zeros()
Bench6=(zeros(length(H),2,2))

Exact6(t)=y0[1] .* cos(sqrt(par[1] / par[2]) .* t) .+ y0[2] ./ sqrt(par[1] / par[2]) .* sin(sqrt(par[1] / par[2]) .* t)

for i in 1:length(H)
        h = H[i]
    Bench6[i,1,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, h=h)).time*10^-9
    Bench6[i,2,1]= mean(@benchmark FDEsolver(F, tSpan, y0, β, par, J = JF, h=h, StopIt = "Convergence")).time*10^-9

    t,y = FDEsolver(F, tSpan, y0, β, par, h=h)
    t1,y1 = FDEsolver(F, tSpan, y0, β, par, J=JF, h=h, StopIt = "Convergence")

    Exact = map(Exact6,t)
    Bench6[i,1,2]=norm(y-Exact)
    Bench6[i,2,2]=norm(y1-Exact)


end
Mdata=DataFrame(CSV.File("Bench6.csv",header=0));
plot(H,Bench6[:,:,1],linewidth=5,xscale = :log,yscale = :log,title="Becnhmark for Example6",yaxis="Time(Sc)",xaxis="Step size", label=["J_PI_PC" "J_PI_IM"])
plot!(H,Mdata[!,1],linewidth=5,xscale = :log,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,2],linewidth=5,xscale = :log,yscale = :log, label="M_PI_IM")
plot(H,Bench6[:,:,2],yscale = :log,xscale = :log,linewidth=5,title="Square norm of the erorrs",yaxis="Errors",xaxis="Step size", label=["PI_PC" "PI_IM"])
plot!(H,Mdata[!,3],linewidth=5,xscale = :log,ls = :dash,yscale = :log, label="M_PI_PC")
plot!(H,Mdata[!,4],linewidth=5,xscale = :log,yscale = :log, ls = :dot,label="M_PI_IM")
