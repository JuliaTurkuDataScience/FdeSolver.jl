using FdeSolver
using Plots

tSpan = [0, 25]
y0 = [34, 6]
β = [0.98, 0.99]

function F(t, n, β, y)

    # α1 = 0.55
    # β1 = 0.028
    # γ = 0.84
    # δ = 0.026

    α1 = 0.55
    β1 = 0.028
    γ = 0.80
    δ = 0.024

    F1 = α1 .* y[n, 1] .- β1 .* y[n, 1] .* y[n, 2]
    F2 = -γ .* y[n, 2] .+ δ .* y[n, 1] .* y[n, 2]

    return [F1, F2]

end

t, Yapp = FDEsolver(F, tSpan, y0, β)

plot(t, Yapp, linewidth = 5, title = "Solution to LV model with 2 FDEs",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")
