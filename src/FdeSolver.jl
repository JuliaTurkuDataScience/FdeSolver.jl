module FdeSolver

using SpecialFunctions
using LinearAlgebra

"""
    greet()
Motivates to get started with FdeSolver
"""
greet() = print("Hey, let's solve some FDEs!")

include("main.jl")
include("main_Jacob.jl")
include("SupFuns.jl")

export(FDEsolver)

end
