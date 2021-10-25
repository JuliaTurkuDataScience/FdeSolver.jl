
const default_values = (2^-6, 2, nothing, "Standard", 10e-6, 100)

struct PositionalArguments

    F::Function
    tSpan::Vector{<:Real}
    y0::Union{Real, Vector{<:Real}, Matrix{<:Real}}
    β::Union{Real, Vector{<:Real}}

    function PositionalArguments(F, tSpan, y0, β)

        if length(tSpan) != 2

            error("tSpan should contain exactly 2 values: one for the initial time and the other for the final time")

        end

        if tSpan[2] <= tSpan[1]

            error("The final time should be greater than the initial time")

        end

        if typeof(β) <: Real

            if β <= 0

                error("β must be a positive integer or float")

            end

        else

            if !isempty(β[β .<= 0])

                error("β must be positive integers or floats")

            end

        end

        new(F, tSpan, y0, β)

    end

end

struct OptionalArguments

    h::Float64
    nc::Int64
    StopIt::String
    tol::Float64
    itmax::Int64

    function OptionalArguments(h, nc, StopIt, tol, itmax)

        if !(StopIt == "Standard" || StopIt == "Convergence")

            error("StopIt can take on either 'Standard' or 'Convergence'")

        end

        new(h, nc, StopIt, tol, itmax)

    end

end

## some structures for standard case
struct Problem
    ic
    f_fun
    problem_size::Int64
    param
    β # we should think about its type
    β_length::Int64
end

struct Method
    bn::Matrix{Float64}
    an::Matrix{Float64}
    a0::Matrix{Float64}
    hα1 # we should think about its type
    hα2 # we should think about its type
    μ::Int64
    μTol::Float64
    r::Int64
    StopIt::String
    itmax::Int64
end

struct Method_fft
    bn_fft::Matrix{ComplexF64}
    an_fft::Matrix{ComplexF64}
    index_fft::Matrix{Int64}
end

struct initial_conditions
    t0::Float64
    y0::Any
    m_β # we should think about its type
    m_β_factorial::Matrix{Int64}
end

## some structures for Jacobian case
struct JProblem
    ic
    f_fun
    problem_size::Int64
    param
    β # we should think about its type
    β_length::Int64
    JF
end

struct JMethod
    an::Matrix{Float64}
    a0::Matrix{Float64}
    hα1 # we should think about its type
    hα2 # we should think about its type
    μ::Int64
    μTol::Float64
    r::Int64
    StopIt::String
    itmax::Int64
end

struct JMethod_fft
    an_fft::Matrix{ComplexF64}
    index_fft::Matrix{Int64}
end
