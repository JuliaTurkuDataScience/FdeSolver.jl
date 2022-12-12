using BenchmarkTools, FdeSolver, FractionalDiffEq, Plots, LinearAlgebra, SpecialFunctions, CSV, DataFrames

## insert data
#it should be based on the directory of CSV files on your computer
push!(LOAD_PATH, "./data_matlab")
# cd("../data_matlab/")
Mdata = CSV.read("BenchNonStiff.csv", DataFrame, header = 0) #it should be based on the directory of CSV files on your computer
## inputs
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

Exact(t) = t.^8 - 3 * t .^ (4 + β / 2) + 9 / 4 * t.^β
## Check the results after the second run
E1 = Float64[]
T1 = Float64[]
E2 = Float64[]
T2 = Float64[]
h = Float64[]

for n in range(3, 8)
    println("n: $n")
    h = 2.0^-n # step size
    t1= @benchmark FDEsolver(F, $(tSpan), $(y0), $(β), $(par) , h=$(h)) seconds=1
    t2= @benchmark FDEsolver(F, $(tSpan), $(y0), $(β), $(par), JF = JF, h=$(h)) seconds=1

    tt1, y1 = FDEsolver(F, tSpan, y0, β, par, h=h)
    tt2, y2= FDEsolver(F, tSpan, y0, β, par, JF = JF, h=h)

    ery1=norm(y1 - map(Exact, tt1), 2)
    ery2=norm(y2 - map(Exact, tt2), 2)

    push!(E1, ery1)
    push!(E2, ery2)

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
end

plot(T1, E1, xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :bottomleft, fc=:transparent,framestyle=:box)
plot!(T2, E2,linewidth = 3, markersize = 5,label = "J-NR", shape = :rect)
plot!(Mdata[:, 1], Mdata[:, 5], linewidth = 3, markersize = 5,label = "M-PI-EX",shape = :rtriangle)
plot!(Mdata[:, 3], Mdata[:, 7], linewidth = 3, markersize = 5,label = "M-PI-IM1", shape = :diamond)
plot!(Mdata[:, 2], Mdata[:, 6], linewidth = 3, markersize = 5,label = "M-PI-PC", shape = :circle)
pNonStiff=plot!(Mdata[:, 4], Mdata[:, 8], linewidth = 3, markersize = 5,label = "M-PI-IM2", shape = :rect
            ,legend=:false)


savefig(pNonStiff,"NonStiff.svg")

#save data
using Tables
CSV.write("NonStiff_E1.csv",  Tables.table(E1))
CSV.write("NonStiff_E2.csv",  Tables.table(E2))
CSV.write("NonStiff_T1.csv",  Tables.table(T1))
CSV.write("NonStiff_T2.csv",  Tables.table(T2))
