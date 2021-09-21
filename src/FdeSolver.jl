# """
# $(DocStringExtensions.README)
# """
module FdeSolver

using SpecialFunctions
using LinearAlgebra

greet() = print("Hey, let's solve some FDEs!")

include("main.jl")
include("main_Jacob.jl")
include("SupFuns.jl")

export(FDEsolver)

end
