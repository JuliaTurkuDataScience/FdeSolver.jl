using BenchmarkTools, FdeSolver, FractionalDiffEq, Plots, LinearAlgebra, SpecialFunctions, CSV, DataFrames

## insert data
#it should be based on the directory of CSV files on your computer
push!(LOAD_PATH, "./data_matlab")
# cd("../data_matlab/")
Mdata = Matrix(CSV.read("BenchHarmonic.csv", DataFrame, header = 0)) #Benchmark from Matlab

## inputs
tSpan = [0, 10]     # [intial time, final time]
y0 = [1 1]             # intial value ([of order 0      of order 1])
α = 2            # order of the derivative
par = [16.0, 4.0] # [spring constant for a mass on a spring, inertial mass]

## Equation
function F(t, x, par)

      K, m = par

      - K ./ m .* x

end
function JF(t, x, par)

      K, m = par

      - K ./ m

end

# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[]
h = Float64[]


for n in range(2, length=6)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computting the time
    t1= @benchmark FDEsolver($(F), $(tSpan), $(y0), $(α), $(par) , h=$(h)) seconds=1
    t2= @benchmark FDEsolver($(F), $(tSpan), $(y0), $(α), $(par), JF = $(JF), h=$(h)) seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
    #computting the error
    t, y1 = FDEsolver(F, tSpan, y0, α, par , h=h)
    _, y2 = FDEsolver(F, tSpan, y0, α, par, JF = JF, h=h)
    #exact solution
    Yex = y0[1] .* map(cos, sqrt(par[1] / par[2]) .* t) .+ y0[2] ./ sqrt(par[1] / par[2]) .* map(sin, sqrt(par[1] / par[2]) .* t)

    ery1=norm(y1 .- Yex,2)
    ery2=norm(y2 .- Yex,2)

    push!(E1, ery1)
    push!(E2, ery2)

end

## plotting
# plot Matlab and FdeSolver outputs
plot(T1, E1, xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :right, fc=:transparent,framestyle=:box)
plot!(T2, E2,linewidth = 3, markersize = 5,label = "J-NR", shape = :rect)
plot!(Mdata[:, 1], Mdata[:, 5], linewidth = 3, markersize = 5,label = "M-PI-EX",shape = :rtriangle)
plot!(Mdata[:, 3], Mdata[:, 7], linewidth = 3, markersize = 5,label = "M-PI-IM1", shape = :diamond)
plot!(Mdata[:, 2], Mdata[:, 6], linewidth = 3, markersize = 5,label = "M-PI-PC", shape = :circle)
pHar=plot!(Mdata[:, 4], Mdata[:, 8], linewidth = 3, markersize = 5,label = "M-PI-IM2", shape = :rect, legend=:false)


savefig(pHar,"Harmonic.svg")

#save data
using Tables
CSV.write("Harmonic_E1.csv",  Tables.table(E1))
CSV.write("Harmonic_E2.csv",  Tables.table(E2))
CSV.write("Harmonic_T1.csv",  Tables.table(T1))
CSV.write("Harmonic_T2.csv",  Tables.table(T2))
