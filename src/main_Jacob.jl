function FDEsolver_Jacob(F, tSpan, y0, β, J, par...; h = 0.01, nc = 3, tol = 10^(-9), itmax = 10)

    # Time discretization
    N::Int64 = cld(tSpan[2] - tSpan[1], h)
    t = tSpan[1] .+ collect(0:N) .* h

    # Enter the initial values
    Y = defineY(N, y0)

    # Number of equations
    Neq = size(y0,1)

    # Calculate T with taylor expansion
    m::Int64 = ceil(β[1])

    for n in 1:N

        T0 = taylor_expansion(tSpan[1], t[n], y0)
        JF = zeros(Neq, Neq)
        YnΨnC0fn = zeros(Neq)

        if n == 1

            # Y1
            Y[2, :] .= T0 .+ h .^ β .* F(t, n, β, Y, par...) ./ Γ(β .+ 1)

            # inverse matrix including Jacobian function
        if Neq == 1

            JF = inv(I(Neq) .- (h .^ β ./ Γ(β .+ 2)) * J(t, n + 1, β, Y, par...))

        else

            JF = inv(I(Neq) .- Diagonal(h .^ β ./ Γ(β .+ 2)) * J(t, n + 1, β, Y, par...))

        end

            for j in 1:nc

            # Y11
            YnΨnC0fn = (Y[2, :] .- T0 .- h .^ β .* β .* F(t, n, β, Y, par...)
            ./ Γ(β .+ 2) .- h .^ β .* F(t, n + 1, β, Y, par...) ./ Γ(β .+ 2))
            Y11 = Y[2, :] - JF * YnΨnC0fn

            Y[2, :] = Y11

            # check the convergence criteria when number of corrections is Inf
            if nc == Inf
                σ = sqrt(sum((Y11 .- Y[2, :]) .^ 2))
                if (σ < tol || j >= itmax)
                    break
                end
            end

            end

        else

            ϕ = Phi(Y, F, β, t, n, par...)

            # Yp
            Y[n + 1, :] .= T0 .+ h .^ β .* (ϕ .- α(0, β) .*
                           F(t, n - 1, β, Y, par...) .+ 2 .* α(0, β) .*
                           F(t, n, β, Y, par...))
          # inverse matrix including Jacobian function
      if Neq == 1
          JF = inv(I(Neq) .- (α(0, β) .* h .^ β) * J(t, n+1, β, Y, par...))
      else
          JF = inv(I(Neq) .- Diagonal(α(0, β) .* h .^ β) * J(t, n+1, β, Y, par...))
      end
          # multiple corrections
          for j in 1:nc
          # Y2
          YnΨnC0fn .= (Y[n+1, :] .- (T0 .+ h .^ β .* ϕ) .-
             h .^ β.* α(0, β) .* F(t, n + 1, β, Y, par...))
          Y2 = Y[n+1, :] - JF * YnΨnC0fn

          Y[n + 1, :] = Y2

           # check the convergence criteria when number of corrections is Inf
            if  nc == Inf
                σ = sqrt(sum((Y2 .- Y[n+1, :]) .^ 2))
                if (σ < tol || j >= itmax)
                    break
                end
            end
            end

        end

    end

    # Output
    t, Y

end
