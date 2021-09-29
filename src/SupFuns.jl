## define initial values ##

function defineY(N, y0, β)

    Y = zeros(N + 1)
    Y[1] = y0
    Y

end

function defineY(N, y0::Vector{<:Real}, β)

    if size(y0) == size(β)

        Y = zeros(N + 1, length(β))
        Y[1, :] .= y0
        Y

    elseif size(y0) != size(β)

        Y = zeros(N + 1)
        Y[1] = y0[1]
        Y

    end

end

function defineY(N, y0::Matrix{<:Real}, β)

    Y = zeros(N + 1, size(y0)[2])
    Y[1, :] .= y0[1, :]
    Y

end

## indexY ##
function indexY(n, Y::Vector{<:Real}, Ynext)

    Y[n + 1] = Ynext
    Y

end

function indexY(n, Y::Matrix{<:Real}, Ynext)

    Y[n + 1, :] = Ynext
    Y

end

## taylor_expansion ##
taylor_expansion(t0, t, y0, β, m) = y0

function taylor_expansion(t0, t, y0::Vector{<:Real}, β, m)

    if size(y0) == size(β)

        return y0

    elseif size(y0) != size(β)

        return sum([(t - t0) ^ k / factorial(k) * y0[k + 1] for k in 0:m - 1])

    end

end

taylor_expansion(t0, t, y0::Matrix{<:Real}, β, m) = sum([(t - t0) ^ k / factorial(k) * y0[:, k + 1] for k in 0:m - 1])

## function ϕ ##
function Phi(y, F, β, t, n, par...)

    alpha = zeros(n - 1)
    fun = zeros(n - 1)

    alpha = map((x) -> α(x, β), n - 1:-1:1)
    fun = map((x) -> F(t, x, β, y, par...), 2:n)

    a_n0(n, β) .* F(t, 1, β, y, par...) .+ sum([alpha[i] .* fun[i] for i in 1:n - 1])

end

a_n0(n, β) = ((n - 1) .^ (β .+ 1) .- n .^ β .* (n .- β .- 1)) ./ Γ(β .+ 2)

## function α ##
function α(n, β)

    if n == 0

        return 1 ./ Γ(β .+ 2)

    else

        return ((n .- 1) .^ (β .+ 1) .- 2 .* n .^ (β .+ 1) .+ (n .+ 1) .^ (β .+ 1)) ./ Γ(β .+ 2)

    end

end

## Gamma function for vectors ##
Γ(b) = map(gamma, b)
