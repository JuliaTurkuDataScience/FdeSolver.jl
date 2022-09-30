# this dummy function is meant to keep the second positional argument of _FDEsolver,
# the Jacobian matrix J, jidden from the user, so that it is not necessary to give
# nothing as the second argument when the Jacobian matrix is not provided.

"""
    FDEsolver(F::Function, tSpan::Vector{<:Real}, y0::Union{Real, Vector{<:Real}, Matrix{<:Real}}, β::Union{Real, Vector{<:Real}}, par...; h = 2^-6, nc = 2, JF = nothing, tol = 10e-6, itmax = 100)

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
- `nc::Int64`: the desired number of corrections for predictor-corrector method, when there is no Jacobian.
- `tol::Float64`: the tolerance of errors, the norm inf of each iteration (for NR method) or correction when nc>10 (for PC method).
- `ìtmax::Int64`: the maximal number of iterations.
"""
function FDEsolver(F, tSpan, y0, β, par...; h = default_values[1], nc = default_values[2], JF = default_values[3], tol = default_values[4], itmax = default_values[5])

    pos_args = PositionalArguments(F, tSpan, y0, β)
    opt_args = OptionalArguments(h, nc, tol, itmax)

    if JF==default_values[3]
        if itmax != default_values[5] && nc == default_values[2]

            @warn "Setting the maximum number of iterations (itmax) is relevant
            only if you use a Jacobian function, otherwise the method is
                predictor-corrector method that no need itmax but you can set
                number of corrections (nc)."

        end

    else # JF!=default_values[3]

        if nc != default_values[2] && itmax == default_values[5]

            @warn "Setting the number of corrections (nc) is relevant
            only if you do NOT use a Jacobian function, otherwise the method is
                predictor-corrector method that no need nc but you can set
                the maximum number of iterations (itmax)."

        end

    end

    _FDEsolver(pos_args, opt_args, JF, par...)

end
