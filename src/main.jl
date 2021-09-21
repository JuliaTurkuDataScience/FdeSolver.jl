"""
    FDEsolver(F, tSpan, y0, β, par...; h = 0.01, nc = 2, tol = 10^(-6), itmax = 30)
Solves fractional differential equations with a predictor-corrector approach.
"""
function FDEsolver(F, tSpan, y0, β, par...; h = 0.01, nc = 2, tol = 10^(-6), itmax = 30)

    # Time discretization
    N::Int64 = cld(tSpan[2] - tSpan[1], h)
    t = tSpan[1] .+ collect(0:N) .* h

    # Enter the initial values
    Y = defineY(N, y0)

    # Calculate T with taylor expansion
    m::Int64 = ceil(β[1])

    for n in 1:N

        T0 = taylor_expansion(tSpan[1], t[n], y0)

        if n == 1

            # Y1
            Y[2, :] .= T0 .+ h .^ β .* F(t, n, β, Y, par...) ./ Γ(β .+ 1)

            for j in 1:nc

                # Y11
                Y11 = T0 .+ h .^ β .* β .* F(t, n, β, Y, par...) ./ Γ(β .+ 2) .+ h .^ β .* F(t, n + 1, β, Y, par...) ./ Γ(β .+ 2)
                σ = sqrt(sum((Y11 .- Y[2, :]) .^ 2))
                Y[2, :] .= Y11

                if (σ < tol || j >= itmax)
                    break
                end

            end

        else

            ϕ = Phi(Y, F, β, t, n, par...)

            # Yp
            Y[n + 1, :] .= T0 .+ h .^ β .* (ϕ .- α(0, β) .* F(t, n - 1, β, Y, par...) .+ 2 .* α(0, β) .* F(t, n, β, Y, par...))

            for j in 1:nc

                # Y2
                Y2 = T0 .+ h .^ β .* (ϕ .+ α(0, β) .* F(t, n + 1, β, Y, par...))
                σ = sqrt(sum((Y2 .- Y[2, :]) .^ 2))
                Y[n + 1, :] .= Y2

                if (σ < tol || j >= itmax)
                    break
                end

            end

        end

    end

    # Output
    t, Y

end
