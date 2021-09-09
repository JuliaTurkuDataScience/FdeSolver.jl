module FdeSolver

using SpecialFunctions

greet() = print("Hey, let's solve some FDEs!")

include("main.jl")
include("defineY.jl")
include("taylor_expansion.jl")
include("a_n0_function.jl")
include("alpha_function.jl")
include("phi_function.jl")
include("gamma_function.jl")

export(FDEsolver)

end
