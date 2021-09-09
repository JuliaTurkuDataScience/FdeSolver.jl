using Revise
push!(LOAD_PATH, "./src")
using FdeSolver

using Plots
using SpecialFunctions

tSpan = [0, 5]
y0 = [1, 0.5, 0.3]
β = [0.5, 0.2, 0.6]

function F(t, n, β, y)

    F1 = 1 / sqrt(pi) * (((y[n, 2] - 0.5) * (y[n, 3] - 0.3))^(1 / 6) + t[n]^(1 / 2))
    F2 = gamma(2.2) * (y[n, 1] - 1)
    F3 = gamma(2.8) / gamma(2.2) * (y[n, 2] - 0.5)

    return [F1, F2, F3]

end

t, Yapp = FDEsolver(F, tSpan, y0, β)

plot(t, Yapp, linewidth = 5, title = "Solution of system 33",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")

plot!(t, t -> (t .+ 1), lw = 3, ls = :dash, label = "Exact solution")
plot!(t, t -> (t.^1.2 .+ 0.5), lw = 3, ls = :dash, label = "Exact solution")
plot!(t, t -> (t.^1.8 .+ 0.3), lw = 3, ls = :dash, label = "Exact solution")
