
module FdeSolver_Jacob

using SpecialFunctions
using LinearAlgebra

greet() = print("Hey, let's solve some FDEs!")

include("main_Jacob.jl")
include("SupFuns.jl")


export(FDEsolver_Jacob)

end
