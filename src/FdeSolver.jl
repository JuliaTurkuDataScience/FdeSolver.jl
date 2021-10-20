"""
    FDEsolver(F, tSpan, y0, β, ::Nothing, par...; h = 0.01, nc = 1, StopIt = "Standard", tol = 10e-10, itmax = 10)

Solves fractional differential equations with a predictor-corrector approach.

    FDEsolver(F, tSpan, y0, β, J, par...; h = 0.01, nc = 1, StopIt = "Standard", tol = 10e-6, itmax = 10)

Additionally takes the Jacobian of the system J to evaluate the solution.

# Arguments
- `F`: the right side of the system of differential equations. It must be expressed
   in the form of a function and return a vector function with the same number of
   entries of order of derivatives. This function can also include a vector of
   parameters: `par...` .
- `tSpan::Vector{<:Real}`: the time span along which computation is performed.
- `y0`: the initial values in the form of a row vector (`Vector{<:Real}`) for ``β <= 1``
   and a matrix (`Matrix{<:Real}`) for ``β > 1``, where each column corresponds to the
   initial values of one differential equation and each row to the order of derivation.
- `β::Vector{<:Real}`: the orders of derivation in the form of a row vector, where
   each element corresponds to the order of one differential equation. It can take
   decimal as well as integer values.
- `J::Matrix{<:Real}`: the Jacobian of F. If not provided, its position within the
   function definition has to be replaced with `nothing`.
- `par...`: additional parameters for the function F.
- `h::Number`: the step size for correction.
- `nc::Int64`: the desired number of corrections.
- `StopIt::String`: the method to stop correction. It can take either "Standard"
   (by default) or "Convergence". In the former case, the function will repeat
   correction as many times as specified in nc; in the latter case, correction will
   stop only when tolerance (tol) or the iteration max (itmax) is reached.
- `tol::Float64`: the tolerance.
- `ìtmax::Int64`: the maximal number of iterations.
"""
module FdeSolver

using SpecialFunctions
using LinearAlgebra
using FFTW

include("main.jl")
include("main_Jacob.jl")
include("SupFuns.jl")
include("SupFuns_Jacob.jl")
include("dummy_fun.jl")

export(FDEsolver)

end
