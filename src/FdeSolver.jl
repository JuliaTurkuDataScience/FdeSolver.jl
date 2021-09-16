# """
# $(DocStringExtensions.README)
# """

module FdeSolver

using SpecialFunctions

"""
    greet()
Motivates to get started with FdeSolver
"""
greet() = print("Hey, let's solve some FDEs!")

include("main.jl")
include("defineY.jl")
include("taylor_expansion.jl")
include("a_n0_function.jl")
include("alpha_function.jl")
include("phi_function.jl")
include("gamma_function.jl")

export(FDEsolver)
export(my_f)

end
