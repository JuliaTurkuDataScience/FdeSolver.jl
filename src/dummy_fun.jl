# this dummy function is meant to keep the second positional argument of _FDEsolver,
# the Jacobian matrix J, jidden from the user, so that it is not necessary to give
# nothing as the second argument when the Jacobian matrix is not provided.

"""
    FDEsolver(F::Function, tSpan::Vector{<:Real}, y0::Union{Real, Vector{<:Real}, Matrix{<:Real}}, β::Union{Real, Vector{<:Real}}, par...; h = 2^-6, nc = 2, JF = nothing, StopIt = "Standard", tol = 10e-6, itmax = 100)

Solves fractional differential equations with a predictor-corrector approach.

# Arguments
- `F::Function`: the right side of the system of differential equations. It must be expressed
   in the form of a function and return a vector function with the same number of
   entries of order of derivatives. This function can also include a vector of
   parameters: `par...` .
- `tSpan::Vector{<:Real}`: the time span along which computation is performed.
- `y0::Union{Real, Vector{<:Real}, Matrix{<:Real}}`: the initial values in the form of a
   row vector (`Vector{<:Real}`) for ``β <= 1`` and a matrix (`Matrix{<:Real}`) for ``β > 1``,
   where each column corresponds to the initial values of one differential equation and
   each row to the order of derivation.
- `β::Union{Real, Vector{<:Real}}`: the orders of derivation in the form of a row vector, where
   each element corresponds to the order of one differential equation. It can take
   decimal as well as integer values.
- `JF::Function`: the Jacobian of F. If not provided, the solver will evaluate the solution
   without the aid of the Jacobian matrix.
- `par...`: additional parameters for the function F.
- `h::Real`: the step size for correction.
- `nc::Int64`: the desired number of corrections.
- `StopIt::String`: the method to stop correction. It can take either "Standard"
   (by default) or "Convergence". In the former case, the function will repeat
   correction as many times as specified in nc; in the latter case, correction will
   stop only when tolerance (tol) or the iteration max (itmax) is reached.
- `tol::Float64`: the tolerance.
- `ìtmax::Int64`: the maximal number of iterations.
"""
function FDEsolver(F, tSpan, y0, β, par...; h = default_values[1], nc = default_values[2], JF = default_values[3], StopIt = default_values[4], tol = default_values[5], itmax = default_values[6])

    pos_args = PositionalArguments(F, tSpan, y0, β)
    opt_args = OptionalArguments(h, nc, StopIt, tol, itmax)

    _FDEsolver(pos_args, opt_args, JF, par...)

end
