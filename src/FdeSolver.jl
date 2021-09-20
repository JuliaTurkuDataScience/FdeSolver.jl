# """
# $(DocStringExtensions.README)
# """

module FdeSolver

using SpecialFunctions
using MittagLeffler

greet() = print("Hey, let's solve some FDEs!")

include("main.jl")
include("SupFuns.jl")


export(FDEsolver)

end
