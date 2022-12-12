using BenchmarkTools, FdeSolver, FractionalDiffEq, Plots, LinearAlgebra, SpecialFunctions, CSV, DataFrames

## insert data
#it should be based on the directory of CSV files on your computer
push!(LOAD_PATH, "./data_matlab")
# cd("../data_matlab/")
Mdata = Matrix(CSV.read("BenchLV.csv", DataFrame, header = 0)) #Benchmark from Matlab
M_ExactLV3 = Matrix(CSV.read("M_Exact_LV3.csv", DataFrame, header = 0)) #Exact from Matlab
## inputs
tSpan = [0, 60]       # [intial time, final time]
y0 = [1,1,1]   # initial values [X1(0),X2(0),X3(0)]
α = [1, .9, .7]          # order of derivatives

## FODE model
par = [3,3,3,5,3,3,3]

function F(t, x, par)

    # parameters
    a1, a2, a3, a4, a5, a6, a7 = par

    # System equation
    dx1 = x[1]*(a1-a2*x[2]-x[3])
    dx2 = x[2]*(1-a3+a4*x[1])
    dx3 = x[3]*(1-a5+a6*x[1]+a7*x[2])

    return [dx1, dx2, dx3]

end

## Jacobian of ODE system
function JF(t, x, par)

    # parameters
    a1, a2, a3, a4, a5, a6, a7 = par

    # System equation
    J11 = a1-2*a2*x[1]-x[2]-x[3]
    J12 = -x[1]
    J13 = -x[1]
    J21 = a4*x[2]
    J22 = 1-a3+a4*x[1]
    J23 = 0
    J31 = a6*x[3]
    J32 = a7*x[3]
    J33 = a6*x[1]-a5+a7*x[2]+1

    J = [J11 J12 J13
         J21 J22 J23
         J31 J32 J33]

    return J

end

##scifracx
function LV!(dx, x, p, t)

    # parameters
    a1, a2, a3, a4, a5, a6, a7 = [3,3,3,5,3,3,3] # this package is not ready for calling the parameters!
    # System equation
    dx[1] = x[1]*(a1-a2*x[2]-x[3])
    dx[2] = x[2]*(1-a3+a4*x[1])
    dx[3] = x[3]*(1-a5+a6*x[1]+a7*x[2])
end

# Tspan = (0, 60)

# y00 = [1;1;1]
# prob = FODESystem(LV!, α, y00, Tspan)

t, Yex= FDEsolver(F, tSpan, y0, α, par, JF = JF, h=2^-10, tol=1e-12) # Solution with a fine step size
Norm2Ex=norm(Yex .- M_ExactLV3',2) # norm2 of solutions with fine stepsize in Matlab and Julia
# Benchmarking
E1 = Float64[];T1 = Float64[];E2 = Float64[];T2 = Float64[]
# E3 = Float64[];T3 = Float64[];
h = Float64[]


for n in range(4,8)
    println("n: $n")# to print out the current step of runing
    h = 2.0^-n #stepsize of computting
        #computting the time
    t1= @benchmark FDEsolver(F, $(tSpan), $(y0), $(α), $(par) , h=$(h), nc=4, tol =1e-8) seconds=1
    t2= @benchmark FDEsolver(F, $(tSpan), $(y0), $(α), $(par), JF = JF, h=$(h), tol=1e-8) seconds=1
    # t3 = @benchmark solve($(prob), $(h), PECE()) seconds=1

    # convert from nano seconds to seconds
    push!(T1, mean(t1).time / 10^9)
    push!(T2, mean(t2).time / 10^9)
    # push!(T3, mean(t3).time / 10^9)
    #computting the error
    _, y1 = FDEsolver(F, tSpan, y0, α, par , h=h, nc=4, tol =1e-8)
    _, y2 = FDEsolver(F, tSpan, y0, α, par, JF = JF, h=h, tol=1e-8)
    # y3 = solve(prob, h, PECE())

    ery1=norm(y1 .- Yex[1:2^(10-n):end,:],2)
    ery2=norm(y2 .- Yex[1:2^(10-n):end,:],2)
    # ery3=norm(y3.u' .- Yex[1:2^(10-n):end,:],2)

    push!(E1, ery1)
    push!(E2, ery2)
    # push!(E3, ery3)

end

## plotting
# plot Matlab and FdeSolver outputs
plot(T1, E1, xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
     label = "J-PC", shape = :circle, xlabel="Execution time (sc, Log)", ylabel="Error: 2-norm (Log)",
     thickness_scaling = 1,legend_position= :bottomleft, fc=:transparent,framestyle=:box)
plot!(T2, E2,linewidth = 3, markersize = 5,label = "J-NR", shape = :rect)
plot!(Mdata[2:end, 1], Mdata[2:end, 5], linewidth = 3, markersize = 5,label = "M-PI-EX",shape = :rtriangle)
plot!(Mdata[:, 3], Mdata[:, 7], linewidth = 3, markersize = 5,label = "M-PI-IM1", shape = :diamond)
plot!(Mdata[:, 2], Mdata[:, 6], linewidth = 3, markersize = 5,label = "M-PI-PC", shape = :circle)
pLV3=plot!(Mdata[:, 4], Mdata[:, 8], linewidth = 3, markersize = 5,label = "M-PI-IM2", shape = :rect
            ,legend=:false)


savefig(pLV3,"LV3.svg")
# plot Scifracx outputs
plotd2=plot!(T3, E3,linewidth = 3,  markersize = 5, label = "J-PECE (FractionalDiffEq.jl)", shape = :diamond)
savefig(plotd2,"LV3_1.svg")


# plot the dynamics
dynamicLV=plot(t[1:100:end],Yex[1:100:end,:], linewidth = 3,
                xlabel="Time", ylabel="Abundance of species" ,
                thickness_scaling = 1 , framestyle=:box, labels=["X1" "X2" "X3"])
savefig(dynamicLV,"dynamicLV.svg")


#save data
using Tables
CSV.write("LV_E1.csv",  Tables.table(E1))
CSV.write("LV_E2.csv",  Tables.table(E2))
CSV.write("LV_T1.csv",  Tables.table(T1))
CSV.write("LV_T2.csv",  Tables.table(T2))

DynLV=[t[1:5:end] Yex[1:5:end,:]]
CSV.write("DynLV.csv",  Tables.table(DynLV))
