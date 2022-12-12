using BenchmarkTools, FdeSolver, FractionalDiffEq, Plots, LinearAlgebra, SpecialFunctions, CSV, DataFrames

## insert data
#it should be based on the directory of CSV files on your computer
push!(LOAD_PATH, " /data_matlab")
# cd("../data_matlab/")
M_Ex1 = CSV.read("RndEx1.csv", DataFrame, header = 1)
M_Ex2 = CSV.read("RndEx2.csv", DataFrame, header = 1)
M_Ex3 = CSV.read("RndEx3.csv", DataFrame, header = 1)
M_Ex4 = CSV.read("RndEx4.csv", DataFrame, header = 1)
## inputs Ex1

# Equation
F1(t, y, par) = (40320 ./ gamma(9 - par) .* t .^ (8 - par) .- 3 .* gamma(5 + par / 2)
           ./ gamma(5 - par / 2) .* t .^ (4 - par / 2) .+ 9/4 * gamma(par + 1) .+
           (3 / 2 .* t .^ (par / 2) .- t .^ 4) .^ 3 .- y .^ (3 / 2))
# Jacobian
JF1(t, y, par) = -(3 / 2) .* y .^ (1 / 2)

## inputs Ex2

# Equation
F2(t, y, par)= par * y
# Jacobian
JF2(t, y, par) = par

## inputs Ex3
#Equation
function F3(t, x, par)

      K, m = par

      - K ./ m .* x

end
function JF3(t, x, par)

      K, m = par

      - K ./ m

end
## inputs Ex4
#FODE model
function SIR(t, y, par)

    # parameters
    β = par[1]    # infection rate
    γ = par[2]    # recovery rate

    S = y[1]   # Susceptible
    I = y[2]   # Infectious
    R = y[3]   # Recovered

    # System equation
    dSdt = - β .* S .* I
    dIdt = β .* S .* I .- γ .* I
    dRdt = γ .* I

    return [dSdt, dIdt, dRdt]

end


# Jacobian of ODE system
function JSIR(t, y, par)

    # parameters
    β = par[1]     # infection rate
    γ = par[2]     # recovery rate

    S = y[1]    # Susceptible
    I = y[2]    # Infectious
    R = y[3]    # Recovered

    # System equation
    J11 = - β * I
    J12 = - β * S
    J13 =  0
    J21 =  β * I
    J22 =  β * S - γ
    J23 =  0
    J31 =  0
    J32 =  γ
    J33 =  0

    J = [J11 J12 J13
         J21 J22 J23
         J31 J32 J33]

    return J

end
function SIR!(dy, y, p, t)

    # System equation
    dy[1] = - p[1] * y[1] * y[2]
    dy[2] = p[1] * y[1] * y[2] - p[2] * y[2]
    dy[3] = p[2] * y[2]
        return nothing
end

## Check the results after the second run
E1 = Float64[]
T1 = Float64[]
E2 = Float64[]
T2 = Float64[]
E3 = Float64[]
T3 = Float64[]

for i in range(1, length(M_Ex1.HR))

    println("i: $i")
    tSpan = [0, 1]     # [intial time, final time]
    h =M_Ex1.HR[i] # step size
    nc=M_Ex1.NcR[i]
    tol=M_Ex1.TolR[i]
    β=M_Ex1.AlphaR[i] # order of the derivative
    par = β #parameter
        if  β>1
            y0 = [0 0]          # intial value
            else
            y0=0             # intial value
        end
    t1= @benchmark FDEsolver(F1, $(tSpan), $(y0), $(β), $(par) , h=$(h), nc=$(nc), tol=$(tol)) seconds=1
    t2= @benchmark FDEsolver(F1, $(tSpan), $(y0), $(β), $(par), JF = JF1, h=$(h), tol=$(tol)) seconds=1

    tt1, y1 = FDEsolver(F1, tSpan, y0, β, par, h=h, nc=nc , tol=tol)
    tt2, y2= FDEsolver(F1, tSpan, y0, β, par, JF = JF1, h=h, tol=tol)

    Exact(t) = t.^8 - 3 * t .^ (4 + β / 2) + 9 / 4 * t.^β
    ery1=norm(y1 - map(Exact, tt1), 2)
    ery2=norm(y2 - map(Exact, tt2), 2)

    push!(E1, ery1)
    push!(E2, ery2)

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
end

for i in range(1, length(M_Ex2.HR))

    println("i: $i")
    h =M_Ex2.HR[i] # step size
    nc=M_Ex2.NcR[i]
    tol=M_Ex2.TolR[i]
    β=M_Ex2.AlphaR[i] # order of the derivative
    par=M_Ex2.LambdaR[i]
    T=M_Ex2.TR[i]; tSpan=[0,T]
    y0 =1
    t1= @benchmark FDEsolver(F2, $(tSpan), $(y0), $(β), $(par) , h=$(h), nc=$(nc), tol=$(tol)) seconds=1
    t2= @benchmark FDEsolver(F2, $(tSpan), $(y0), $(β), $(par), JF = JF2, h=$(h), tol=$(tol)) seconds=1

    tt1, y1 = FDEsolver(F2, tSpan, y0, β, par, h=h, nc=nc , tol=tol)
    tt2, y2= FDEsolver(F2, tSpan, y0, β, par, JF = JF2, h=h, tol=tol)

    #exact solution: mittag-leffler
    λ=par
    Exact(t) = mittleff(β, λ * t .^ β)
    ery1=norm(y1 - map(Exact, tt1), 2)
    ery2=norm(y2 - map(Exact, tt2), 2)

    push!(E1, ery1)
    push!(E2, ery2)

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
end


for i in range(1, length(M_Ex3.HR))

    println("i: $i")
    h =M_Ex3.HR[i] # step size
    nc=M_Ex3.NcR[i]
    tol=M_Ex3.TolR[i]
    β=2 # order of the derivative
    par=[16.0, 4.0]
    T=M_Ex3.TR[i]; tSpan=[0,T]
    y0 = [1 1]

    t1= @benchmark FDEsolver(F3, $(tSpan), $(y0), $(β), $(par) , h=$(h), nc=$(nc), tol=$(tol)) seconds=1
    t2= @benchmark FDEsolver(F3, $(tSpan), $(y0), $(β), $(par), JF = JF3, h=$(h), tol=$(tol)) seconds=1

    t, y1 = FDEsolver(F3, tSpan, y0, β, par, h=h, nc=nc , tol=tol)
    _, y2= FDEsolver(F3, tSpan, y0, β, par, JF = JF3, h=h, tol=tol)

    #exact solution:
    Yex = y0[1] .* map(cos, sqrt(par[1] / par[2]) .* t) .+ y0[2] ./ sqrt(par[1] / par[2]) .* map(sin, sqrt(par[1] / par[2]) .* t)

    ery1=norm(y1 .- Yex, 2)
    ery2=norm(y2 .- Yex, 2)

    push!(E1, ery1)
    push!(E2, ery2)

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
end


for i in range(1, length(M_Ex4.HR))

    println("i: $i")
    h=2.0^(- M_Ex4.HR[i])# step size
    nc=M_Ex4.NcR[i]
    tol=M_Ex4.TolR[i]
    β=[M_Ex4.alp1R[i],M_Ex4.alp2R[i],M_Ex4.alp3R[i]] # order of the derivative
    par=[M_Ex4.BetaR[i], M_Ex4.GamaR[i]]
    T=M_Ex4.TR[i]; tSpan=[0,T]
    y0 = [1-M_Ex4.II0[i], M_Ex4.II0[i], 0]

    Tspan = (0, T)
    y00 = [1-M_Ex4.II0[i]; M_Ex4.II0[i]; 0]
    prob = FODESystem(SIR!, β, y00, Tspan, par)

    t1= @benchmark FDEsolver(SIR, $(tSpan), $(y0), $(β), $(par) , h=$(h), nc=$(nc), tol=$(tol)) seconds=1
    t2= @benchmark FDEsolver(SIR, $(tSpan), $(y0), $(β), $(par), JF = JSIR, h=$(h), tol=$(tol)) seconds=1
    t3 = @benchmark solve($(prob), $(h), PECE()) seconds=1

    t, y1 = FDEsolver(SIR, tSpan, y0, β, par, h=h, nc=nc , tol=tol)
    _, y2= FDEsolver(SIR, tSpan, y0, β, par, JF = JSIR, h=h, tol=tol)
    y3 = solve(prob, h, PECE())

    #fine solution:
    t, Yex= FDEsolver(SIR, tSpan, y0, β, par, JF = JSIR, h=2^-10, tol=1e-12) # Solution with a fine step size

    ery1=norm(y1 .- Yex[1:2^(10-M_Ex4.HR[i]):end,:],2)
    ery2=norm(y2 .- Yex[1:2^(10-M_Ex4.HR[i]):end,:],2)
    ery3=norm(y3.u' .- Yex[1:2^(10-M_Ex4.HR[i]):end,:],2)

    push!(E1, ery1)
    push!(E2, ery2)
    push!(E3, ery3)

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
    push!(T3, mean(t3).time / 10^9)
end


# plotting
scatter(T1, E1, xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :bottomleft, fc=:transparent,framestyle=:box)
scatter!(T2, E2,linewidth = 3, markersize = 5,label = "J-NR", shape = :rect)
scatter!(T3, E3,linewidth = 3,  markersize = 5, label = "J-PECE", shape = :circle)

Mt1=vcat(M_Ex1.Bench1_1,M_Ex2.Bench2_1,M_Ex3.Bench3_1,M_Ex4.Bench4_1)
Merr1=vcat(M_Ex1.Bench1_5,M_Ex2.Bench2_5,M_Ex3.Bench3_5,M_Ex4.Bench4_5)

Mt2=vcat(M_Ex1.Bench1_2,M_Ex2.Bench2_2,M_Ex3.Bench3_2,M_Ex4.Bench4_2)
Merr2=vcat(M_Ex1.Bench1_6,M_Ex2.Bench2_6,M_Ex3.Bench3_6,M_Ex4.Bench4_6)

Mt3=vcat(M_Ex1.Bench1_3,M_Ex2.Bench2_3,M_Ex3.Bench3_3,M_Ex4.Bench4_3)
Merr3=vcat(M_Ex1.Bench1_7,M_Ex2.Bench2_7,M_Ex3.Bench3_7,M_Ex4.Bench4_7)

Mt4=vcat(M_Ex1.Bench1_4,M_Ex2.Bench2_4,M_Ex3.Bench3_4,M_Ex4.Bench4_4)
Merr4=vcat(M_Ex1.Bench1_8,M_Ex2.Bench2_8,M_Ex3.Bench3_8,M_Ex4.Bench4_8)

scatter!(Mt1, Merr1, linewidth = 3, markersize = 5,label = "M-PI-EX",shape = :rtriangle)
scatter!(Mt3, Merr3, linewidth = 3, markersize = 5,label = "M-PI-IM1", shape = :diamond)
scatter!(Mt2, Merr2, linewidth = 3, markersize = 5,label = "M-PI-PC", shape = :circle)
scatter!(Mt4, Merr4, linewidth = 3, markersize = 5,label = "M-PI-IM2", shape = :rect)



#save data
using Tables
CSV.write("ErrRnd1.csv",  Tables.table(E1))
CSV.write("ErrRnd2.csv",  Tables.table(E2))
CSV.write("ErrRnd3.csv",  Tables.table(E3))
CSV.write("tRnd1.csv",  Tables.table(T1))
CSV.write("tRnd2.csv",  Tables.table(T2))
CSV.write("tRnd3.csv",  Tables.table(T3))
#
# J_E1 = CSV.read("ErrRnd1.csv", DataFrame, header = 1)
# J_E2 = CSV.read("ErrRnd2.csv", DataFrame, header = 1)
# J_T1 = CSV.read("tRnd1.csv", DataFrame, header = 1)
# J_T2 = CSV.read("tRnd2.csv", DataFrame, header = 1)
# J_E3 = CSV.read("ErrRnd3.csv", DataFrame, header = 1)
# J_T3 = CSV.read("tRnd3.csv", DataFrame, header = 1)
#
# # plotting
# scatter(J_T1[:,1], J_E1[:,1], xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
#      label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
#      thickness_scaling = 1,legend_position= :bottomleft, fc=:transparent,framestyle=:box)
# scatter!(J_T2[:,1], J_E2[:,1],linewidth = 3, markersize = 5,label = "J-NR", shape = :rect)
# scatter!(J_T3[:,1], J_E3[:,1],linewidth = 3,  markersize = 5, label = "J-PECE", shape = :circle)
#
# Mt1=vcat(M_Ex1.Bench1_1,M_Ex2.Bench2_1,M_Ex3.Bench3_1,M_Ex4.Bench4_1)
# Merr1=vcat(M_Ex1.Bench1_5,M_Ex2.Bench2_5,M_Ex3.Bench3_5,M_Ex4.Bench4_5)
#
# Mt2=vcat(M_Ex1.Bench1_2,M_Ex2.Bench2_2,M_Ex3.Bench3_2,M_Ex4.Bench4_2)
# Merr2=vcat(M_Ex1.Bench1_6,M_Ex2.Bench2_6,M_Ex3.Bench3_6,M_Ex4.Bench4_6)
#
# Mt3=vcat(M_Ex1.Bench1_3,M_Ex2.Bench2_3,M_Ex3.Bench3_3,M_Ex4.Bench4_3)
# Merr3=vcat(M_Ex1.Bench1_7,M_Ex2.Bench2_7,M_Ex3.Bench3_7,M_Ex4.Bench4_7)
#
# Mt4=vcat(M_Ex1.Bench1_4,M_Ex2.Bench2_4,M_Ex3.Bench3_4,M_Ex4.Bench4_4)
# Merr4=vcat(M_Ex1.Bench1_8,M_Ex2.Bench2_8,M_Ex3.Bench3_8,M_Ex4.Bench4_8)
#
# scatter!(Mt1, Merr1, linewidth = 3, markersize = 5,label = "M-PI-EX",shape = :rtriangle)
# scatter!(Mt3, Merr3, linewidth = 3, markersize = 5,label = "M-PI-IM1", shape = :diamond)
# scatter!(Mt2, Merr2, linewidth = 3, markersize = 5,label = "M-PI-PC", shape = :circle)
# plttt=scatter!(Mt4, Merr4, linewidth = 3, markersize = 5,label = "M-PI-IM2", shape = :rect)
#
# savefig(plttt,"plttt.png")
#
