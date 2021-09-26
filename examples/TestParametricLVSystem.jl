using FdeSolver
using Plots

tSpan = [0, 25]
y0 = [34, 6]
β = [0.98, 0.99]
par = [0.55, 0.028, 0.80, 0.024]

function F(t, n, β, y, par)

    F1 = par[1] .* y[n, 1] .- par[2] .* y[n, 1] .* y[n, 2]
    F2 = - par[3] .* y[n, 2] .+ par[4] .* y[n, 1] .* y[n, 2]

    [F1, F2]

end

t, Yapp = FDEsolver(F, tSpan, y0, β, par)

plot(t, Yapp, linewidth = 5, title = "Solution to LV model with 2 FDEs",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")
