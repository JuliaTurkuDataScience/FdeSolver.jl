module FdeSolver

import SpecialFunctions: gamma
import LinearAlgebra: Diagonal, norm, I
import FFTW: fft, ifft

include("main.jl")
include("main_Jacob.jl")
include("SupFuns.jl")
include("SupFuns_Jacob.jl")
include("dummy_fun.jl")
include("structs.jl")

export(FDEsolver)

end
