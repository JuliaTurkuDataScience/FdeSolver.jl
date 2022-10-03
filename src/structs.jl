
const default_values = (2^-6, 1, nothing, 1e-6, 100)

SZ=Base.size # we have to define a new function size() because size is defined as a value!
    # catch exceptions for tSpan

struct PositionalArguments

    F::Function
    tSpan::Vector{<:Real}
    y0::Union{Real, Vector{<:Real}, Matrix{<:Real}}
    β::Union{Real, Vector{<:Real}}

    function PositionalArguments(F, tSpan, y0, β)


        if length(tSpan) != 2

            error("tSpan should contain exactly 2 values: first value is for the initial time and the second one for the final time.")

        end

        if tSpan[2] <= tSpan[1]

            error("The final time should be greater than the initial time.")

        end

        # catch exceptions for y0
        if typeof(y0) <: Real

            if length(β) != 1

                error("Either the order of derivatives (β) has too many values or the initial values (y0) has too few.")

            end
            if length(y0) < Int64(ceil(β))

                error("The number of elements in initial values (y0) should equal the smallest integer greater than (or equal to) order of derivative (ceil(β)). For instance, if β = [1.2], its next integer is 2 and y0 should have 2 elements for the initial values.")

            end

        elseif typeof(y0) <: Vector{<:Real}

            if SZ(y0,2) != Int64(ceil(maximum(β)))

                error("The number of elements in initial values (y0) should equal the smallest integer greater than or equal to order of derivative (ceil(max(β)). For instance, if β = [1.2] or β = [0.5, 1.2], its next integer is 2 and y0 should have 2 elements for the initial values.")

            end
            if (SZ(y0) != SZ(β) && length(β) != 1)

                error("The number of elements in initial values (y0) should match that in order of derivatives (β).")

            elseif (SZ(y0) != SZ(β) && length(β) == 1)

                error("Initial values of higher-order derivatives should not be given as a vector (y0 = [y0', y0'', y0''', ...]). Instead, they should be given as a matrix (y0 = [y0' y0'' y0''' ...]).")

                # if length(y0) != Int64(ceil(β))
                #
                #     error("The number of elements in initial values (y0) should equal the smallest integer greater than or equal to order of derivative (ceil(max(β)). For instance, if β = [1.2] or β = [0.5, 1.2], its next integer is 2 and y0 should have 2 elements for the initial values.")
                #
                # end

            end

        else # typeof(y0) <: Matrix{<:Real}

            if length(β) == 1

                if SZ(y0, 2) != Int64(ceil(β))

                    error("The number of elements in initial values (y0) should equal the smallest integer greater than or equal to order of derivative (ceil(max(β)). For instance, if β = [1.2] or β = [0.5, 1.2], its next integer is 2 and y0 should have 2 elements for the initial values.")

                end

            else # length(β) != 1

                if SZ(y0, 2) != Int64(ceil(maximum(β)))

                    error("The number of columns in y0 should equal the next integer of β. For instance, if β = 1.2, its next integer is 2 and y0 should have 2 columns of initial values.")

                end

                if SZ(y0, 2) != length(β)

                    error("The number of rows in y0 should match the length of β, that is, the number of equations to solve.")

                end

            end

        end

        # catch exceptions for β
        if typeof(β) <: Real

            if β <= 0

                error("β must be a positive integer or float.")

            end

        else # typeof(β) == Vector{<:Real}

            if !isempty(β[β .<= 0])

                error("β must be positive integers or floats.")

            end

        end

        if !(typeof(y0) <: Matrix{<:Real})

            Y0 = zeros(SZ(y0, 1), Int64(ceil(maximum(β))))
            Y0 .= y0

        else

            Y0 = y0

        end

        new(F, tSpan, Y0, β)

    end

end

struct OptionalArguments

    h::Float64
    nc::Int64
    tol::Float64
    itmax::Int64

    function OptionalArguments(h, nc, tol, itmax)

        # catch exceptions for h
        if h <= 0

            error("The step size h must be positive.")

        end

        # give warnings for tol, itmax and nc
        if nc != default_values[2] && itmax != default_values[5]

                @warn "Setting the maximum number of iterations (itmax) makes sense only if you use a Jacobian function
                    and setting the number of corrections (nc) make sense only if you do NOT use a Jacobian function."
        end

        new(h, nc, tol, itmax)

    end

end

## some structures for standard case
struct Problem

    ic
    f_fun
    problem_size::Int64
    param
    β::Union{Real, Vector{<:Real}, Matrix{<:Real}}
    β_length::Int64

    function Problem(ic, f_fun, problem_size, param, β, β_length)

        if (β_length > 1 && problem_size != β_length)

            error("The size of the problem obtained from the initial conditions (the number of rows of y0) is not compatible with the number of fractional orders.")

        end

        new(ic, f_fun, problem_size, param, β, β_length)

    end

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

    function initial_conditions(t0, y0, m_β, m_β_factorial)

        if SZ(y0, 2) < maximum(m_β)

            error("There are not enough initial conditions to solve a system of order β.")

        end

        new(t0, y0, m_β, m_β_factorial)

    end

end

## some structures for Jacobian case
struct JProblem

    ic
    f_fun
    problem_size::Int64
    param
    β::Union{Real, Vector{<:Real}, Matrix{<:Real}} # matrix for a system with orders greater than 1
    β_length::Int64
    JF

    function JProblem(ic, f_fun, problem_size, param, β, β_length, JF)

        if problem_size != β_length

            error("The size of the problem obtained from the initial conditions (the number of rows of y0) is not compatible with the number of fractional orders.")

        end

        new(ic, f_fun, problem_size, param, β, β_length, JF)

    end

end

struct JMethod

    an::Matrix{Float64}
    a0::Matrix{Float64}
    hα1 # we should think about its type
    hα2 # we should think about its type
    μTol::Float64
    r::Int64
    itmax::Int64

end

struct JMethod_fft

    an_fft::Matrix{ComplexF64}
    index_fft::Matrix{Int64}

end
