using BenchmarkTools, FdeSolver, FractionalDiffEq, Plots, LinearAlgebra, SpecialFunctions, CSV, DataFrames

## insert data
#it should be based on the directory of CSV files on your computer
push!(LOAD_PATH, "./data_matlab")
# cd("../data_matlab/")
Mdata = Matrix(CSV.read("BenchSIR.csv", DataFrame, header = 0)) #Benchmark from Matlab
M_ExactSIR = Matrix(CSV.read("M_ExactSIR.csv", DataFrame, header = 0)) #Exact from Matlab
## inputs
I0 = 0.1             # intial value of infected
tSpan = [0, 100]       # [intial time, final time]
y0 = [1 - I0, I0, 0]   # initial values [S0,I0,R0]
α = [.9, .6, .7]          # order of derivatives

## FODE model
par = [0.4, 0.04] # parameters [β, recovery rate]

function F(t, y, par)

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

## Jacobian of ODE system
function JF(t, y, par)

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

##scifracx
function SIR!(dy, y, p, t)

    # System equation
    dy[1] = - .4 * y[1] * y[2]
    dy[2] = .4 * y[1] * y[2] - .04 * y[2]
    dy[3] = .04 * y[2]
        return dy
end

Tspan = (0, 100)

y00 = [1 - I0; I0; 0]
prob = FODESystem(SIR!, α, y00, Tspan)

t, Yex= FDEsolver(F, tSpan, y0, α, par, JF = JF, h=2^-10, tol=1e-12) # Solution with a fine step size
Norm2Ex=norm(Yex .- M_ExactSIR',2) # norm2 of solutions with fine stepsize in Matlab and Julia
# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[]
E3 = Float64[];T3 = Float64[];E4 = Float64[];T4 = Float64[]
E5 = Float64[];T5 = Float64[];E6 = Float64[];T6 = Float64[]
E7 = Float64[];T7 = Float64[];E8 = Float64[];T8 = Float64[]
E9 = Float64[];T9 = Float64[];h = Float64[]


for n in range(2, length=6)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computting the time
    t1= @benchmark FDEsolver(F, $(tSpan), $(y0), $(α), $(par) , h=$(h)) seconds=1
    t2= @benchmark FDEsolver(F, $(tSpan), $(y0), $(α), $(par), JF = JF, h=$(h), tol=1e-8) seconds=1
    t3= @benchmark solve($(prob), $(h), $(GL())) seconds=1
    t4= @benchmark solve($(prob), $(h), $(PIEX())) seconds=1
    t5= @benchmark solve($(prob), $(h), $(NonLinearAlg())) seconds=1
    t6 = @benchmark solve($(prob), $(h), FLMMBDF()) seconds=1
    t7 = @benchmark solve($(prob), $(h), FLMMNewtonGregory()) seconds=1
    t8 = @benchmark solve($(prob), $(h), FLMMTrap()) seconds=1
    t9 = @benchmark solve($(prob), $(h), PECE()) seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
    push!(T3, mean(t3).time / 10^9)
    push!(T4, mean(t4).time / 10^9)
    push!(T5, mean(t5).time / 10^9)
    push!(T6, mean(t6).time / 10^9)
    push!(T7, mean(t7).time / 10^9)
    push!(T8, mean(t8).time / 10^9)
    push!(T9, mean(t9).time / 10^9)
    #computting the error
    _, y1 = FDEsolver(F, tSpan, y0, α, par , h=h)
    _, y2 = FDEsolver(F, tSpan, y0, α, par, JF = JF, h=h, tol=1e-8)
    y3 =  solve(prob, h, GL())
    y4 =  solve(prob, h, PIEX())
    y5 =  solve(prob, h, NonLinearAlg())
    y6 = solve(prob, h, FLMMBDF())
    y7 = solve(prob, h, FLMMNewtonGregory())
    y8 = solve(prob, h, FLMMTrap())
    y9 = solve(prob, h, PECE())

    ery1=norm(y1 .- Yex[1:2^(10-n):end,:],2)
    ery2=norm(y2 .- Yex[1:2^(10-n):end,:],2)
    ery3=norm(y3.u' .- Yex[1:2^(10-n):end,:],2)
    ery4=norm(y4.u' .- Yex[1:2^(10-n):end,:],2)
    ery5=norm(y5.u' .- Yex[1:2^(10-n):end,:],2)
    ery6=norm(y6.u' .- Yex[1:2^(10-n):end,:],2)
    ery7=norm(y7.u' .- Yex[1:2^(10-n):end,:],2)
    ery8=norm(y8.u' .- Yex[1:2^(10-n):end,:],2)
    ery9=norm(y9.u' .- Yex[1:2^(10-n):end,:],2)

    push!(E1, ery1)
    push!(E2, ery2)
    push!(E3, ery3)
    push!(E4, ery4)
    push!(E5, ery5)
    push!(E6, ery6)
    push!(E7, ery7)
    push!(E8, ery8)
    push!(E9, ery9)

end

## plotting
# plot Matlab and FdeSolver outputs
plot(T1, E1, xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :right,fc=:transparent,framestyle=:box)
plot!(T2, E2,linewidth = 3, markersize = 5,label = "J-NR", shape = :rect)
plot!(Mdata[2:end, 1], Mdata[2:end, 5], linewidth = 3, markersize = 5,label = "M-PI-EX",shape = :rtriangle)
plot!(Mdata[:, 3], Mdata[:, 7], linewidth = 3, markersize = 5,label = "M-PI-IM1", shape = :diamond)
plot!(Mdata[:, 2], Mdata[:, 6], linewidth = 3, markersize = 5,label = "M-PI-PC", shape = :circle)
pSIR=plot!(Mdata[:, 4], Mdata[:, 8], linewidth = 3, markersize = 5,label = "M-PI-IM2", shape = :rect)


# savefig(pSIR,"pSIR.svg")

# plot Scifracx outputs
plot!(T5, E5, xscale = :log, yscale = :log, linewidth = 3, label = "J-NonLinearAlg", shape = :star5)
plot!(T9, E9,linewidth = 3,  markersize = 5, label = "J-PECE", shape = :circle)
plot!(T4, E4,linewidth = 3, markersize = 5, label = "J-PI-EX", shape = :hexagon)
plot!(T6, E6,linewidth = 3,markersize = 5,  label = "J-FLMMBDF", shape = :utriangle)
plot!(T7, E7,linewidth = 3,markersize = 5,  label = "J-FLMMNewtonG", shape = :dtriangle)
# pSIR1=plot!(T8, E8,linewidth = 2, markersize = 5, label = "J-FLMMTrap", shape = :pentagon, legend=:outerright)
# savefig(pSIR1,"pSIR1.svg")
pSIR2=plot!(T8, E8,linewidth = 3, markersize = 5, label = "J-FLMMTrap", shape = :pentagon, legend=false)

savefig(pSIR2,"pSIR2.svg")

#save data
using Tables
CSV.write("SIR_E1.csv",  Tables.table(E1))
CSV.write("SIR_E2.csv",  Tables.table(E2))
# CSV.write("SIR_E3.csv",  Tables.table(E3))
CSV.write("SIR_E4.csv",  Tables.table(E4))
CSV.write("SIR_E5.csv",  Tables.table(E5))
CSV.write("SIR_E6.csv",  Tables.table(E6))
CSV.write("SIR_E7.csv",  Tables.table(E7))
CSV.write("SIR_E8.csv",  Tables.table(E8))
CSV.write("SIR_E9.csv",  Tables.table(E9))
CSV.write("SIR_T1.csv",  Tables.table(T1))
CSV.write("SIR_T2.csv",  Tables.table(T2))
CSV.write("SIR_T3.csv",  Tables.table(T3))
CSV.write("SIR_T4.csv",  Tables.table(T4))
CSV.write("SIR_T5.csv",  Tables.table(T5))
CSV.write("SIR_T6.csv",  Tables.table(T6))
CSV.write("SIR_T7.csv",  Tables.table(T7))
CSV.write("SIR_T8.csv",  Tables.table(T8))
CSV.write("SIR_T9.csv",  Tables.table(T9))
