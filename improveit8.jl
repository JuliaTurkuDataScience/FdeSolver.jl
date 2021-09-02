function improveit8(F, tSpan, y0 , β; h = 0.01, c = 4)::Tuple{Vector{Float64}, Matrix{Float64}}

    # Time discretization
    N::Int64 = round((tSpan[2] - tSpan[1]) / h)
    t = tSpan[1] .+ collect(0:N) .* h

    # Enter the initial values
    Y = zeros(N + 1, length(y0))

    Y[1, :] .= y0

        for n in 1:N

            if n == 1

                # Y1
                Y1 = zeros(1, length(y0))
                Y1[1, :] .= y0 .+ h .^ β .* F(t, n, β, Y) ./ Γ(β .+ 1)

                for j in 1:c

                    # Y11
                    Y[2, :] .= y0 .+ h .^ β .* β .* F(t, n, β, Y) ./ Γ(β .+ 2) .+ h .^ β .* F(t, n, β, Y1) ./ Γ(β .+ 2)

                end

            else

                ϕ = Phi(Y, F, β, t, n)

                # Yp
                Y[n + 1, :] .= y0 .+ h .^ β .* (ϕ .- α(0, β) .* F(t, n - 1, β, Y) .+ 2 .* α(0, β) .* F(t, n, β, Y))

                for j in 1:c

                    # Y2
                    Y[n + 1, :] .= y0 .+ h .^ β .* (ϕ .+ α(0, β) .* F(t, n + 1, β, Y))

                end

            end

        end

    # Output
    t, Y

end
