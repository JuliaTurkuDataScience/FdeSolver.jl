using BenchmarkTools, FdeSolver, FractionalDiffEq, Plots, LinearAlgebra, SpecialFunctions, CSV, DataFrames

## insert data
#it should be based on the directory of CSV files on your computer
push!(LOAD_PATH, "./data_matlab")
# cd("./data_matlab/")
Mdata = Matrix(CSV.read("BenchStiff.csv", DataFrame, header = 0)) #Benchmark from Matlab

## inputs
tSpan = [0, 5]     # [intial time, final time]
y0 = 1             # intial value
β = 0.8            # order of the derivative

λ = -10
par = λ
F(t, y, par)= par * y
JF(t, y, par) = par

#exact solution: mittag-leffler
Exact(t) = mittleff(β, λ * t .^ β)
# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[]
h = Float64[]


for n in range(3, 8)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computing the time
    t1= @benchmark FDEsolver(F, $(tSpan), $(y0), $(β), $(par) , h=$(h), nc=4) seconds=1
    t2= @benchmark FDEsolver(F, $(tSpan), $(y0), $(β), $(par), JF = JF, h=$(h)) seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
    #computing the error
    tt1, y1 = FDEsolver(F, tSpan, y0, β, par , h=h, nc=4)
    _, y2 = FDEsolver(F, tSpan, y0, β, par, JF = JF, h=h)

    ery1=norm(y1 - map(Exact, tt1), 2)
    ery2=norm(y2 - map(Exact, tt1), 2)

    push!(E1, ery1)
    push!(E2, ery2)

end

## plotting
# plot Matlab and FdeSolver outputs
plot(T1, E1, xscale = :log, yscale = :log, linewidth = 2, markersize = 5,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :right, c=:blue,fc=:transparent,framestyle=:box, mc=:white)
plot!(T2, E2,linewidth = 2, markersize = 5,label = "J-NR", shape = :rect, color = :blue, mc=:white)
plot!(Mdata[:, 1], Mdata[:, 5], linewidth = 2, markersize = 5,label = "M-PI-EX",shape = :rtriangle, color = :red, mc=:white)
plot!(Mdata[:, 2], Mdata[:, 6], linewidth = 2, markersize = 5,label = "M-PI-PC", shape = :circle, color = :red, mc=:white)
plot!(Mdata[:, 3], Mdata[:, 7], linewidth = 2, markersize = 5,label = "M-PI-IM1", shape = :diamond, color = :red, mc=:white)
pStiff1=plot!(Mdata[:, 4], Mdata[:, 8], linewidth = 2, markersize = 5,label = "M-PI-IM2", shape = :rect, color = :red, mc=:white)


savefig(pStiff1,"Stiff1.svg")


plot(T2, E2, linewidth = 2, markersize = 5,
     label = "J-NR", shape = :rect, xlabel="Execution time (sc)", ylabel="Error: 2-norm",
     thickness_scaling = 1,legend_position= :right, c=:blue,fc=:transparent,framestyle=:box, mc=:white)
pStiff2=plot!(Mdata[:, 4], Mdata[:, 8], linewidth = 2, markersize = 5,label = "M-PI-IM2", shape = :rect, color = :red, mc=:white)


savefig(pStiff2,"Stiff2.svg")

plot(T1[2:end], E1[2:end], xscale = :log, yscale = :log, linewidth = 3,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :topright,fc=:transparent,framestyle=:box)
plot!(T2[1:end], E2[1:end],linewidth = 3,label = "J-NR", shape = :rect)
plot!(Mdata[2:end, 1], Mdata[2:end, 5], linewidth = 3, label = "M-PI-EX",shape = :rtriangle)
plot!(Mdata[1:end, 3], Mdata[1:end, 7], linewidth = 3, label = "M-PI-IM1", shape = :diamond)
plot!(Mdata[2:end, 2], Mdata[2:end, 6], linewidth =3, label = "M-PI-PC", shape = :circle)
pStiff3=plot!(Mdata[1:end, 4], Mdata[1:end, 8], linewidth = 3, markersize = 5,label = "M-PI-IM2",shape = :rect, legend=:false)

savefig(pStiff3,"Stiff3.svg")


#save data
using Tables
CSV.write("Stiff_E1.csv",  Tables.table(E1))
CSV.write("Stiff_E2.csv",  Tables.table(E2))
CSV.write("Stiff_T1.csv",  Tables.table(T1))
CSV.write("Stiff_T2.csv",  Tables.table(T2))
